# Script to classify Planet VNIR imagery to water mask
# Nils Rietze nils.rietze@uzh.ch 16 November 2023

library(terra)
library(tidyverse)
library(tidyterra)
library(sf)
library(spatialEco)
library(caret)
library(randomForest)
library(ggplot2)
library(purrr)
library(cowplot)
library(scico)
library(dplyr)
library(parallel)
library(pbapply)

set.seed(10)

# Load ... ----

# ...areas of interest
aois <- read_sf('./data/geodata/feature_layers/aoi_wv/aois.shp') %>%
  mutate(., id = 1:nrow(.))

# ... training polygons
training_polygons <- vect('./data/geodata/feature_layers/training_polygons/training_polygons_burn_area.shp') %>% 
  project('EPSG:32655')

# ... list of raster files
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

# Config ----

# Define bands to extract (SuperDove = ('Coastal Blue','Blue','Green I','Green','Yellow','Red','Red Edge','NIR'))
sel_bands <- c('Blue','Green','Red','NIR')

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
    mutate(iswater = factor(if_else(is.na(LandCover), 'no_water','water'),labels = c('no_water','water') ) ) 
  
  # Extract raster data
  raster_vals <- terra::extract(rast_refl, sel_features)[,c('ID',sel_bands)]
  
  raster_focal_vals <- terra::extract(rast_focal, sel_features)[,c('ID','gcc_sd','rcc_sd','bcc_sd','ndwi_sd','ndvi_sd')]
  
  # Change column_names
  # names(raster_focal_vals)[-1] <- paste0(names(raster_focal_vals)[-1], "_sd")
  
  # Combine raster value tables
  raster_vals_all <- cbind(raster_vals, 
                           select(raster_focal_vals,-ID) ) 
  
  # Join labels and raster values
  training_data <- sel_features %>%
    as.data.frame() %>%
    select(ID, iswater) %>%
    full_join(raster_vals_all) 
  
  # compute chromatic coordinates
  training_data <- training_data %>%
    mutate(ndvi = (Red - NIR) / (Red + NIR),
           ndwi = (Green - NIR) / (Green + NIR),
           rcc = Red / (Red + Blue + Green),
           bcc = Blue / (Red + Blue + Green),
           gcc = Green / (Red + Blue + Green),
    )
  
  
  return(training_data)
  
}

# Raster value extraction ----

# Extract raster values in polygons
regex_name <- 'Kosukhino'
fn_raster_pre <- raster_files_pre[grepl(regex_name,raster_files_pre)]
fn_raster_post <- raster_files_post[grepl(regex_name,raster_files_post) & 
                                    grepl('PS2',raster_files_post)]

training_data_pre <- get_raster_values(fn_raster_pre, training_polygons, sel_bands)
training_data_post <- get_raster_values(fn_raster_post, training_polygons, sel_bands)

# Compute differences in normalized bands & indices and add to post-fire data
cols_to_diff <- c("ndvi", "ndwi", "rcc", "bcc", "gcc")

training_data <- training_data_post %>%
  mutate(across(all_of(cols_to_diff), ~ . - training_data_pre[[cur_column()]], .names = "D{.col}"))

# extract site name from filename
site <- unlist(strsplit( basename(fn_raster_post),'_')) [2]

npoly <- nrow(training_polygons)
print(training_data %>% group_by(ID, iswater) %>% tally(), n = npoly)

training_sample <- training_data %>%
  group_by(ID,iswater) %>%
  na.omit() %>%
  sample_n(50,replace = TRUE)

training_sample %>% group_by(ID, iswater) %>% tally()

# Calculate separability
calculate_separability <- function(variable_name) {
    separability(
      filter(training_sample, iswater == 'uniswater') %>% ungroup() %>% select(!!variable_name ),
      filter(training_sample, iswater == 'iswater') %>% ungroup() %>% select(!!variable_name )
    ) %>%
      select(-!!variable_name , -!!paste0(variable_name,".1"))
}

bands_of_interest <- c("Red", "Green", "Blue","NIR",
                       "gcc_sd","rcc_sd","bcc_sd","ndwi_sd","ndvi_sd",
                       "rcc", "gcc", "bcc","ndvi", "ndwi",
                       "Dndvi","Dndwi","Drcc","Dbcc","Dgcc" )

sep_stats <- map_dfr(bands_of_interest, calculate_separability)

# Plot separability
plot_separability <- function(band="Red"){
  training_sample$band_to_plot <- pull(training_sample, !!band)
  y_pos <- max(training_sample$band_to_plot) * .9
  
  p <- ggplot(training_sample) +
    geom_violin(aes(x = iswater,
                    y = band_to_plot, 
                    fill = iswater)) +
    geom_text(data = sep_stats[band,], 
              aes(x = 2, y = y_pos, 
                  label = paste("TD =", 
                                formatC(TD, format = "f", digits = 2)),
                  hjust = 0,
                  fontface = "bold")) +
    labs(x = 'iswater class',y = band) +
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

plot_grid(plotlist = plot_list, 
          nrow = nr,
          labels = paste0(letters[1:n], ")")) 
  # save_plot(paste0("figures/band_separability_",site,"_",year,".png"),
  #           .,
  #           nrow = nr,
  #           ncol = nc,
  #           base_asp = 1.2,
  #           bg = "white")

# Get top 5 bands based on TD separability 
top_rows <- rownames(sep_stats)[order(-sep_stats$TD)[1:5]]
cat("Top 5 Row Names with Highest TD Values:\n", top_rows, "\n")

# Split into training and validation
training_sample$id <- 1:nrow(training_sample)
training <- training_sample %>%
  group_by(iswater) %>%
  slice_sample(prop = 0.8)
validation <- filter(training_sample, !(id %in% training$id)) 

# Train random forest model on all variables
rf_fit <- randomForest(formula(paste('iswater ~',paste(bands_of_interest,collapse = "+"))),
                       data = select(training, -id))
top_5_gini <- varImp(rf_fit, scale = FALSE) %>% 
  top_n(5, Overall) %>%
  rownames()
# top_rows <- top_5_gini

# Train random forest model on top 5 
rf_fit <- randomForest(formula(paste('iswater ~',paste(top_rows,collapse = "+"))),
                       data = select(training, -id))

# Test accuracy on validation split
test_preds <- predict(rf_fit, newdata = validation,
                      type = "response")

# Export confusion matrix
fn_conf_matrix <- paste0("tables/rf_confusion_matrix_water_",site,"_",year,".txt")

file_connection <- file(fn_conf_matrix)

confusionMatrix(data = test_preds, validation$iswater) %>%
  print() %>%
  capture.output() %>%
  writeLines(file_connection)
close(file_connection)

# Prepare raster indices
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
  } else if(band == "ndwi"){
    r <- (img[["Green"]] - img[["NIR"]]) / (img[["Green"]] + img[["NIR"]])
  }
  names(r) <- band
  
  return(r)
}
 
# Load predictors
vnir_bands <- c("Red","Green","Blue","NIR")
band_indices <- c("ndwi","ndvi","rcc","bcc","gcc")
band_diffs <- c("Dndvi","Dndwi","Drcc","Dbcc","Dgcc")
focal_indices <- c("Blue_sd","Green_sd","Red_sd","rcc_sd","gcc_sd","bcc_sd","ndwi_sd","ndvi_sd")

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
cat("Writing raster...\n")
writeRaster(preds,
            filename = paste0("data/geodata/raster/water_area/planet/",
                              site,"_",year,"_water_area.tif"),
            overwrite = T
)

# Export predictors as raster
writeRaster(predictors,
            filename = paste0("data/geodata/raster/water_area/planet/",
                              site,"_",year,"_predictors.tif"),
            overwrite = T
)
