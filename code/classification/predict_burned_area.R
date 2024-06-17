# Script to classify Planet VNIR imagery to burned area maps
# Nils Rietze nils.rietze@uzh.ch 16 November 2023

library(terra)
library(sf)
library(spatialEco)
library(tidyverse)
library(tidyterra)
library(caret)
library(randomForest)
library(cowplot)
library(scico)
library(parallel)
library(pbapply)

#  Configuration ----
set.seed(10)

predict_burn_raster <- TRUE  # predict and export burned area raster ?
use_top5_TD <- TRUE           # use top 5 predictors based on transformed divergence (=TRUE) or GINI(=FALSE)
predict_all_aois <- TRUE      # model all AOIs or just one?
plot_results <- TRUE          # Plot figures generated during prediction?
save_plots <- TRUE            # Save those figures?

# °°°°°°°°°°°°°
# Load ... ----

##...areas of interest ----
aois <- read_sf('./data/geodata/feature_layers/aoi_wv/aois_analysis.geojson') %>%
  mutate(., id = 1:nrow(.))

# Provide a part of the AOI name for processing
if (predict_all_aois){
  aoi_names <- aois$site
} else {
  aoi_names <- 'LargeScarCenter'
}

##... training polygons----
training_polygons <- vect('./data/geodata/feature_layers/training_polygons/training_polygons_burn_area.shp') %>% 
  project('EPSG:32655')

##... polygon_masks----
poly_mask <- vect('data/geodata/feature_layers/planet_masks.shp')

##... list of raster files----
pre_year <- 2019
post_year <- 2020

raster_files_pre <- list.files(sprintf('C:/data/8_planet/%s/cropped/',pre_year),
                           pattern = '*composite.tif$',
                           full.names = T
)
raster_files_post <- list.files(sprintf('C:/data/8_planet/%s/cropped/',post_year),
                           pattern = '*composite.tif$',
                           full.names = T
)

# Set up functions ----

# function to compute the GEMI
GetGEMI <- function(NIR, Red){
  term1 <- (2 * (NIR^2 - Red^2) + 1.5 * NIR + 0.5 * Red) /
    (NIR + Red + 0.5)
  gemi <- term1 * (1 - 0.25 * term1) - ((Red - 0.125) / (1 - Red))
  return(gemi)
}

# Function to extract raster values within training polygons
get_raster_values <- function(fn_raster, poly_mask, features, sel_bands){
  cat("\n Extracting raster values: ", fn_raster, "\n")
  
  # Load raster data
  rast_refl <- rast(fn_raster) %>% mask(.,poly_mask,inverse = T)
  
  # Load focal raster (NOT USED)
  fn_raster_focal <- gsub("composite", paste0("focal_5_sd"), fn_raster)
  rast_focal <- rast(fn_raster_focal) %>% mask(.,poly_mask,inverse = T)
  names(rast_focal) <- paste0(names(rast_focal),'_sd')
  
  # Select training polygons within this raster
  sel_features <- features %>% 
    terra::intersect(ext(rast_refl)) %>% 
    mutate(.,ID = 1:nrow(.)) %>% 
    mutate(.,Burned = factor(Burned,labels = c('unburned','burned')))
  
  # Extract raster data
  raster_vals <- terra::extract(rast_refl, sel_features)[,c('ID',sel_bands)]
  
  raster_focal_vals <- terra::extract(rast_focal, 
                                      sel_features)[,c('ID',
                                                       'gcc_sd',
                                                       'rcc_sd',
                                                       'bcc_sd',
                                                       'bai_sd',
                                                       'ndvi_sd',
                                                       'ndwi_sd')]
  
  # Change column_names
  # names(raster_focal_vals)[-1] <- paste0(names(raster_focal_vals)[-1], "_sd")
  
  # Combine raster value tables
  raster_vals_all <- cbind(raster_vals, 
                           select(raster_focal_vals,-ID) ) 
  
  # Join labels and raster values
  training_data <- sel_features %>%
    as.data.frame() %>%
    select(ID, Burned) %>%
    full_join(raster_vals_all) 
  
  # compute chromatic coordinates
  training_data <- training_data %>%
    mutate(ndvi = (NIR - Red) / (NIR + Red),
           ndwi = (Green - NIR) / (Green + NIR),
           bai = 1 / ( (0.1 - Red)^2 + (0.06 - NIR)^2 ),
           gemi = GetGEMI(NIR,Red),
           rcc = Red / (Red + Blue + Green),
           bcc = Blue / (Red + Blue + Green),
           gcc = Green / (Red + Blue + Green),
    )
  
  return(training_data)
  
}

predict_burned_area <- function(aoi_name,
                                raster_files_pre, 
                                raster_files_post, 
                                training_polygons,
                                poly_mask){
  ### 1. configure inputs ----
  
  cat(sprintf("Processing %s ... \n",aoi_name))
  
  # Get filenames in that aoi
  fn_raster_pre <- raster_files_pre[grepl(aoi_name,raster_files_pre)]
  fn_raster_post <- raster_files_post[grepl(aoi_name,raster_files_post)]
  
  # extract site name from filename
  site <- unlist(strsplit( basename(fn_raster_post),'_')) [2]
  
  if (use_top5_TD){
    pred_choice <- "top5TD"
  } else {
    pred_choice <- "top5GINI"
  }
  
  # define output filenames
  FNAME_OUT <- paste0("data/geodata/raster/burned_area/planet/",
                      site,"_burned_area_",pred_choice,".tif")
  
  FNAME_PREDS_OUT <- paste0("data/geodata/raster/burned_area/planet/",
                            site,"_predictors_",pred_choice,".tif")
  
  # check if file exists
  if (file.exists(FNAME_OUT)){
    cat("File exists. Skipping burned area prediction for this site. \n")
    return(NULL)
  }
  
  ### 2. gather raster values ----
  
  # Define bands to extract (SuperDove = ('Coastal Blue','Blue','Green I','Green','Yellow','Red','Red Edge','NIR'))
  sel_bands <- c('Blue','Green','Red','NIR')
  
  training_data_pre <- get_raster_values(fn_raster_pre, poly_mask, training_polygons, sel_bands)
  training_data_post <- get_raster_values(fn_raster_post, poly_mask, training_polygons, sel_bands)

  # Compute differences in normalized bands & indices and add to post-fire data
  cols_to_diff <- c("ndvi","ndwi", "bai", "rcc", "bcc", "gcc","gemi")
  training_data <- training_data_post %>%
    mutate(across(all_of(cols_to_diff), ~ . - training_data_pre[[cur_column()]], .names = "D{.col}"))
  
  ### 3. format training data ----
  npoly <- nrow(training_polygons)
  # print(training_data %>% group_by(ID, Burned) %>% tally(), n = npoly)
  
  #  sample 50 pixels with replacement
  training_sample <- training_data %>%
    group_by(ID,Burned) %>%
    na.omit() %>%
    sample_n(50,replace = TRUE)
  
  # Split into training and validation
  training_sample$id <- 1:nrow(training_sample)
  
  groups <- training_sample %>%
    group_by(ID, Burned) %>%
    summarize(n = n(), .groups = 'drop')
  train_indices <- createDataPartition(groups$Burned, p = 0.8, list = FALSE)
  train_ids <- groups$ID[train_indices]
  val_ids <- setdiff(groups$ID, train_ids)
  
  training <- training_sample  %>% filter(ID %in% train_ids)
  validation <- training_sample %>% filter(ID %in% val_ids)
  
### 4. Calculate separability ----
  calculate_separability <- function(variable_name) {
    separability(
      filter(training_sample, Burned == 'unburned') %>% ungroup() %>% select(!!variable_name ),
      filter(training_sample, Burned == 'burned') %>% ungroup() %>% select(!!variable_name )
    ) %>%
      select(-!!variable_name , -!!paste0(variable_name,".1"))
  }
  
  bands_of_interest <- c("Red", "Green", "Blue","NIR",
                         # "gcc_sd","rcc_sd","bcc_sd","bai_sd","ndvi_sd","ndwi_sd",
                         "rcc", "gcc", "bcc","ndvi","ndwi", "bai","gemi",
                         "Dndvi","Dndwi","Dbai","Dgemi","Drcc","Dbcc","Dgcc" )
  
  sep_stats <- map_dfr(bands_of_interest, calculate_separability)
  
  if (plot_results){
    # Plot separability
      plot_separability <- function(band="Red"){
        training_sample$band_to_plot <- pull(training_sample, !!band)
        y_pos <- max(training_sample$band_to_plot) * .9
        
        p <- ggplot(training_sample) +
          geom_violin(aes(x = Burned,
                          y = band_to_plot, 
                          fill = Burned)) +
          geom_text(data = sep_stats[band,], 
                    aes(x = 2, y = y_pos, 
                        label = paste("TD =", 
                                      formatC(TD, format = "f", digits = 2)),
                        hjust = 0,
                        fontface = "bold")) +
          labs(x = 'Burned class',y = band) +
          theme_cowplot() +
          theme(legend.position = "none",
                strip.background = element_rect(fill = "white"))
        
        return(p)
      }
    
    plot_list <- map(bands_of_interest,plot_separability)
    
    # dynamically adjust grid layout
    n <- length(plot_list)
    nr <- ifelse(n / 2 > 4,3,2)
    nc <- ifelse(ceiling(n%%nr) == 0, n / nr, ceiling(n / nr))
    
    p <- cowplot::plot_grid(plotlist = plot_list, nrow = nr,
                            labels = paste0(letters[1:n], ")"))
    
    if (save_plots){
      cowplot::save_plot(filename = paste0("figures/band_separability_",site,".png"),
                         p,
                         nrow = nr,
                         ncol = nc,
                         base_asp = 1.2,
                         bg = "white")
    }
  } # end of separability plot
  
  # define the 5 predictor variables
  if (use_top5_TD){
    # Get top 5 bands based on TD separability 
    top_rows <- rownames(sep_stats)[order(-sep_stats$TD)[1:5]]
    cat("Top 5 Row Names with highest TD Values:\n", top_rows, "\n")
    
  } else {
    # Get top 5 bands based on mean GINI
    rf_fit <- randomForest(formula(paste('Burned ~',paste(bands_of_interest,collapse = "+"))),
                       data = select(training, -id))
    top_rows <- varImp(rf_fit, scale = FALSE) %>% 
    top_n(5, Overall) %>%
    rownames()
    cat("Top 5 Row Names with highest GINI:\n", top_rows, "\n")
  }
  
  ### 5. Train random forest model on top 5 ----
  rf_fit <- randomForest(formula(paste('Burned ~',paste(top_rows,collapse = "+"))),
                         data = select(training, -id))
  
  ### 6. Accuracy tests & confusion matrix ----
  # Test accuracy on validation split
  test_preds <- predict(rf_fit, newdata = validation,type = "response")
  
  # Export confusion matrix
  fn_conf_matrix <- paste0("tables/rf_confusion_matrix_",site,"_",pred_choice,".txt")
  
  file_connection <- file(fn_conf_matrix)
  
  confusionMatrix(data = test_preds, validation$Burned) %>%
    print() %>%
    capture.output() %>%
    writeLines(file_connection) 
  
  cat("writing confusion matrix ... \n")
  
  write(paste(rownames(sep_stats[order(-sep_stats$TD),]),
              sep_stats[order(-sep_stats$TD),'TD'],
              collapse = '\n', sep = ':'),
      file = fn_conf_matrix, append = TRUE)
  
  ### 7. Prepare rasters for prediction ----
  
  # Only predict and produce raster if wanted
  if (predict_burn_raster){
    
    # Set up functions
    get_CC <- function(img, band) {
      if (band == "rcc") band_og <- "Red"
      if (band == "gcc") band_og <- "Green"
      if (band == "bcc") band_og <- "Blue"
      img_cc <- img[[band_og]] / (img[["Red"]] + img[["Green"]] + img[["Blue"]])
      return(img_cc)
    }
    
    compute_image <- function(fn_raster,band){
      img <- rast(fn_raster)
      
      if(band %in% c("rcc", "gcc", "bcc")) {
        r <- get_CC(img,band)
      } else if(band == "ndvi"){
        r <- (img[["NIR"]] - img[["Red"]]) / (img[["NIR"]] + img[["Red"]])
      } else if(band == "ndwi"){
        r <- (img[["Green"]] - img[["NIR"]]) / (img[["Green"]] + img[["NIR"]])
      } else if(band == "bai"){
        r <- 1 / ( (0.1 - img[["Red"]])^2 + (0.06 - img[["NIR"]])^2 )
      } else if(band == "gemi"){
        r <- GetGEMI(img[["NIR"]],img[["Red"]])
      }
      names(r) <- band
      
      return(r)
    }
    
    # Load predictors
    vnir_bands <- c("Red","Green","Blue","NIR")
    band_indices <- c("bai","ndvi","ndwi","gemi",
                      "rcc","bcc","gcc")
    band_diffs <- c("Dndvi","Dndwi","Dbai","Dgemi",
                    "Drcc","Dbcc","Dgcc")
    focal_indices <- c("Blue_sd","Green_sd","Red_sd",
                       "rcc_sd","gcc_sd","bcc_sd",
                       "bai_sd","ndvi_sd","ndwi_sd")
    
    if (any(vnir_bands %in% top_rows) ) {
      bands <- intersect(vnir_bands,top_rows)
      cat(paste("Retrieving", bands, "for prediction.\n") )
      
      # Load raster
      rast_vnir <- rast(fn_raster_post)[[bands]]
      
    } else {rast_vnir <- NULL}
    if (any(band_indices %in% top_rows) ) {
      bands <- intersect(band_indices,top_rows)
      cat(paste("Retrieving", bands, "for prediction.\n") )
      
      # Load raster
      rast_indices <- rast( map(bands, ~compute_image(fn_raster_post, .)) )
      names(rast_indices) <- bands
      
    } else {rast_indices <- NULL}
    if (any(band_diffs %in% top_rows) ){
      bands <- intersect(band_diffs,top_rows)
      cat(paste("Retrieving", bands, "for prediction.\n") )
      
      # Extract which band indices to subtract
      temp_bands <- gsub("D", "", bands)
      
      rast_temp_pre <- rast( map(temp_bands, ~compute_image(fn_raster_pre, .)) )
      rast_temp_post <- rast( map(temp_bands, ~compute_image(fn_raster_post, .)) )
      
      # Load raster
      rast_diffs <- rast_temp_post - rast_temp_pre
      names(rast_diffs) <- bands
      
    } else {rast_diffs <- NULL}
    if(any(focal_indices %in% top_rows) ){
      bands <- intersect(focal_indices,top_rows)
      cat(paste("Retrieving", bands, "for prediction.\n") )
      
      window <- 5
      f_stat <- 'sd'
      
      fn_raster_focal <- gsub("composite", paste0("focal_", window, "_", f_stat), fn_raster_post)
      
      # Load raster
      rast_focal <- rast(fn_raster_focal)
      names(rast_focal) <- paste0(names(rast_focal),'_sd')
      rast_focal <- rast_focal[bands]
    } else {rast_focal <- NULL}
    
    # merge all predictor rasters
    rast_list <- list(rast_indices, rast_vnir, rast_diffs, rast_focal) # first create list
    
    # remove any NULL entity from the list
    rast_list <- rast_list[!sapply(rast_list, is.null)]
    
    # Concatenate the raster objects
    predictors <- do.call(c, rast_list) %>% 
      subst(., -Inf, NA)
    
    # Predict entire image
    cat("Predicting raster...\n")
    preds <- terra::predict(predictors, rf_fit,na.rm = TRUE)
    
    # Export rasters
    cat("Writing raster...\n")
    writeRaster(preds,
                filename = FNAME_OUT,
                overwrite = T
    )
    
    writeRaster(predictors,
                filename = FNAME_PREDS_OUT,
                overwrite = T
    )
    
  }
  cat("done.")
  
  return(NULL)
} # end of function


# Execute function ----

if(predict_all_aois){
  # setup cluster
  cl <- makeCluster(7)
  
  clusterEvalQ(cl, c(library(terra),
                     library(tidyverse),
                     library(tidyterra),
                     library(caret),
                     library(sf),
                     library(spatialEco),
                     library(randomForest),
                     library(cowplot),
                     library(scico))
               )
  clusterExport(cl, c("predict_burned_area","get_raster_values"), 
                envir=environment())
  
  # Crop images
  pblapply(
    aoi_names,
    raster_files_pre = raster_files_pre,
    raster_files_post = raster_files_post,
    training_polygons = training_polygons,
    poly_mask = poly_mask,
    FUN = predict_burned_area,
    cl = 7)
  
  # Stop cluster
  stopCluster(cl)
} else {
  aoi_name <- "Berelech"
  predict_burned_area(aoi_name, raster_files_pre, raster_files_post, training_polygons, poly_mask)
}
