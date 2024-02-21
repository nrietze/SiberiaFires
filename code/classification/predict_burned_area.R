# Script to classify Planet VNIR imagery to burned area maps
# Nils Rietze nils.rietze@uzh.ch 16 November 2023

library(terra)
library(tidyverse)
library(tidyterra)
library(sf)
library(spatialEco)
library(caret)
library(randomForest)
library(cowplot)
library(scico)
library(parallel)
library(pbapply)

#  Configuration ----
set.seed(10)

# Should a burned area raster be produced?
predict_burn_raster <- TRUE

# Should the RF model use the 5 predictors with highest transformed divergence index?
use_top5_TD <- TRUE # if FALSE, the top 5 GINI predictors will be used

# Should an individual AOI be processed or all?
predict_all_aois <- FALSE

# Provide a part of the AOI name for processing
if (predict_all_aois){
  aoi_names <- c('Center','Control','Kosukhino','Berelech')
} else {
  aoi_names <- 'Center'
}

# Define bands to extract (SuperDove = ('Coastal Blue','Blue','Green I','Green','Yellow','Red','Red Edge','NIR'))
sel_bands <- c('Blue','Green','Red','NIR')

# Prepare config table
config <- list(aoi_names,
               sel_bands)

# Plot figures generated during prediction?
plot_results <- TRUE

# Save those figures?
save_plots <- TRUE

# °°°°°°°°°°°°°
# Load ... ----

##...areas of interest ----
aois <- read_sf('./data/geodata/feature_layers/aoi_wv/aois.shp') %>%
  mutate(., id = 1:nrow(.))

##... training polygons----
training_polygons <- vect('./data/geodata/feature_layers/training_polygons/training_polygons_burn_area.shp') %>% 
  project('EPSG:32655')

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

# Function to extract raster values within training polygons
get_raster_values <- function(fn_raster, features, sel_bands){
  cat("\nProcessing", fn_raster, "\n")
  
  window <- 5
  f_stat <- 'sd'
  
  fn_raster_focal <- gsub("composite", paste0("focal_", window, "_", f_stat), fn_raster)
  
  # Load raster data
  rast_refl <- rast(fn_raster)
  rast_focal <- rast(fn_raster_focal)
  names(rast_focal) <- paste0(names(rast_focal),'_sd')
  
  # Select training polygons within this raster
  sel_features <- features %>% 
    terra::intersect(ext(rast_refl)) %>% 
    mutate(.,ID = 1:nrow(.)) %>% 
    mutate(.,Burned = factor(Burned,labels = c('unburned','burned')))
  
  # Extract raster data
  raster_vals <- terra::extract(rast_refl, sel_features)[,c('ID',sel_bands)]
  
  raster_focal_vals <- terra::extract(rast_focal, sel_features)[,c('ID','gcc_sd','rcc_sd','bcc_sd','bai_sd','ndvi_sd')]
  
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
    mutate(ndvi = (Red - NIR) / (Red + NIR),
           bai = 1 / ( (0.1 - Red)^2 + (0.06 - NIR)^2 ),
           rcc = Red / (Red + Blue + Green),
           bcc = Blue / (Red + Blue + Green),
           gcc = Green / (Red + Blue + Green),
    )
  
  
  return(training_data)
  
}


predict_burned_area <- function(config, raster_files_pre, raster_files_post, training_polygons){
  ### 1. configure inputs ----
  aoi_name <- config[[1]] # get aoi name
  
  # Get filenames in that aoi
  fn_raster_pre <- raster_files_pre[grepl(aoi_name,raster_files_pre)]
  fn_raster_post <- raster_files_post[grepl(aoi_name,raster_files_post) & 
                                        grepl('PS2',raster_files_post)]
  
  # extract site name from filename
  site <- unlist(strsplit( basename(fn_raster_post),'_')) [2]
  
  cat(sprintf("Predict burned area for: %s.",site))
  
  ### 2. gather raster values ----
  training_data_pre <- get_raster_values(fn_raster_pre, training_polygons, sel_bands)
  training_data_post <- get_raster_values(fn_raster_post, training_polygons, sel_bands)

  # Compute differences in normalized bands & indices and add to post-fire data
  cols_to_diff <- c("ndvi", "bai", "rcc", "bcc", "gcc")
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
  training <- training_sample %>%
    group_by(Burned) %>%
    slice_sample(prop = 0.8)
  validation <- filter(training_sample, !(id %in% training$id)) 
  
  ### 4. Calculate separability ----
  calculate_separability <- function(variable_name) {
    separability(
      filter(training_sample, Burned == 'unburned') %>% ungroup() %>% select(!!variable_name ),
      filter(training_sample, Burned == 'burned') %>% ungroup() %>% select(!!variable_name )
    ) %>%
      select(-!!variable_name , -!!paste0(variable_name,".1"))
  }
  
  bands_of_interest <- c("Red", "Green", "Blue","NIR",
                         "gcc_sd","rcc_sd","bcc_sd","bai_sd","ndvi_sd",
                         "rcc", "gcc", "bcc","ndvi", "bai",
                         "Dndvi","Dbai","Drcc","Dbcc","Dgcc" )
  
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
    pred_choice <- "top5TD"
    
  } else {
    # Get top 5 bands based on mean GINI
    rf_fit <- randomForest(formula(paste('Burned ~',paste(bands_of_interest,collapse = "+"))),
                       data = select(training, -id))
    top_rows <- varImp(rf_fit, scale = FALSE) %>% 
    top_n(5, Overall) %>%
    rownames()
    cat("Top 5 Row Names with highest GINI:\n", top_rows, "\n")
    pred_choice <- "top5GINI"
  }
  # top_rows_aoi2 <- c("NIR","bai","Dndvi","Dbai","bai_sd")
  
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
  close(file_connection)
  
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
        r <- (img[["Red"]] - img[["NIR"]]) / (img[["Red"]] + img[["NIR"]])
      } else if(band == "bai"){
        r <- 1 / ( (0.1 - img[["Red"]])^2 + (0.06 - img[["NIR"]])^2 )
      }
      names(r) <- band
      
      return(r)
    }
    
    # Load predictors
    vnir_bands <- c("Red","Green","Blue","NIR")
    band_indices <- c("bai","ndvi","rcc","bcc","gcc")
    band_diffs <- c("Dndvi","Dbai","Drcc","Dbcc","Dgcc")
    focal_indices <- c("Blue_sd","Green_sd","Red_sd","rcc_sd","gcc_sd","bcc_sd","bai_sd","ndvi_sd")
    
    if (any(vnir_bands %in% top_rows) ) {
      bands <- intersect(vnir_bands,top_rows)
      cat(paste("Retrieving", bands, "for prediction.\n") )
      
      # Load raster
      rast_vnir <- rast(fn_raster_post)[[bands]]
      
    } 
    if (any(band_indices %in% top_rows) ) {
      bands <- intersect(band_indices,top_rows)
      cat(paste("Retrieving", bands, "for prediction.\n") )
      
      # Load raster
      rast_indices <- rast( map(bands, ~compute_image(fn_raster_post, .)) )
      
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
    
    # Predict entire image
    predictors <- c(rast_vnir,
                    rast_indices,
                    rast_diffs,
                    rast_focal)
    
    cat("Predicting raster...\n")
    preds <- terra::predict(predictors, rf_fit)
    
    # Export rasters
    cat("Writing raster...\n")
    writeRaster(preds,
                filename = paste0("data/geodata/raster/burned_area/planet/",
                                  site,"_burned_area_",pred_choice,".tif"),
                overwrite = T
    )
    
    writeRaster(predictors,
                filename = paste0("data/geodata/raster/burned_area/planet/",
                                  site,"_predictors_",pred_choice,".tif"),
                overwrite = T
    )
    
  }
  cat("done.")
  
  return(NULL)
} # end of function

predict_burned_area(config, raster_files_pre, raster_files_post, training_polygons)

