# Script to classify Planet VNIR imagery to burned area maps
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

# Load areas of interest
aois <- read_sf('./data/geodata/feature_layers/aoi_wv/aois.shp') %>%
  mutate(., id = 1:nrow(.))

# Load training polygons
training_polygons <- read_sf('./data/geodata/feature_layers/training_polygons/training_polygons_burn_area.shp') 

# Load list of raster files
raster_path <- 'C:/data/8_planet/2020/cropped/'

raster_files <- list.files(raster_path,
                           pattern = '*.tif$',
                           full.names = T
)

# Bands to extract (SuperDove = ('Coastal Blue','Blue','Green I','Green','Yellow','Red','Red Edge','NIR'))
sel_bands <- c(2,4,6,8) # B, G, R, NIR

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
    vect() %>% 
    project('EPSG:32655') %>% 
    terra::intersect(ext(rast_refl)) %>% 
    mutate(.,ID = 1:nrow(.)) %>% 
    mutate(.,Burned = factor(Burned,labels = c('unburned','burned')))
  
  # Extract raster data
  raster_vals <- terra::extract(rast_refl, sel_features)[,c(1,sel_bands+1)]
  
  raster_focal_vals <- terra::extract(rast_focal, sel_features)[,c('ID','rcc','bcc','bai','ndvi')]
  
  # Change column_names
  # names(raster_focal_vals)[-1] <- paste0(names(raster_focal_vals)[-1], "_sd")
  
  # Combine raster value tables
  raster_vals_all <- cbind(raster_vals, 
                           select(raster_focal_vals, -ID))
  
  training_data <- sel_features %>%
    as.data.frame() %>%
    select(ID, Burned) %>%
    full_join(raster_vals_all) 
  
}

training_sample <- training_data %>%
  group_by(ID,Burned) %>%
  na.omit() %>%
  sample_n(10)

training_sample %>% group_by(ID, Burned) %>% tally()

training_sample <- training_sample %>% 
  mutate(ndvi = (Red - NIR) / (Red + NIR),
         bai = 1 / ( (0.1 - Red)^2 + (0.06 - NIR)^2 ),
         rcc = Red / (Red + Blue + Green),
         bcc = Blue / (Red + Blue + Green),
         gcc = Green / (Red + Blue + Green),
  )

# Calculate separability
calculate_separability <- function(variable_name) {
    separability(
      filter(training_sample, Burned == 'unburned') %>% ungroup() %>% select(!!variable_name ),
      filter(training_sample, Burned == 'burned') %>% ungroup() %>% select(!!variable_name )
    ) %>%
      select(-!!variable_name , -!!paste0(variable_name,".1"))
}

bands_of_interest <- c("Red", "Green", "Blue",
                       "rcc_sd","bcc_sd","bai_sd","ndvi_sd",
                       "rcc", "gcc", "bcc","ndvi", "bai")

sep_stats <- map_dfr(bands_of_interest, calculate_separability)

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

plot_grid(plotlist = plot_list, 
          nrow = nr,
          labels = paste0(letters[1:n], ")")) %>%
  save_plot("figures/band_separability.png",
            .,
            nrow = nr,
            ncol = nc,
            base_asp = 1.2,
            bg = "white")

# Get top 5 bands based on TD separabilities 
top_rows <- rownames(sep_stats)[order(-sep_stats$TD)[1:5]]
cat("Top 5 Row Names with Highest TD Values:\n", top_rows, "\n")

# Split into training and validation
training_sample$id <- 1:nrow(training_sample)
training <- training_sample %>%
  group_by(Burned) %>%
  slice_sample(prop = 0.8)
validation <- filter(training_sample, !(id %in% training$id)) 


# Train random forest model
rf_fit <- randomForest(formula(paste('Burned ~',paste(top_rows,collapse = "+"))),
                       data = select(training, -id))

# Test accuracy on validation split
test_preds <- predict(rf_fit, newdata = validation,
                      type = "response")

# Export confusion matrix
file_connection <- file("tables/rf_confusion_matrix.txt")
confusionMatrix(data = test_preds, validation$Burned) %>%
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

compute_image <- function(band){
  if(band %in% c("rcc", "gcc", "bcc")) {
    r <- get_CC(rast_refl,band)
  } else if(band == "ndvi"){
    r <- (rast_refl[["Red"]] - rast_refl[["NIR"]]) / (rast_refl[["Red"]] + rast_refl[["NIR"]])
  } else if(band == "bai"){
    r <- 1 / ( (0.1 - rast_refl[["Red"]])^2 + (0.06 - rast_refl[["NIR"]])^2 )
  }
  names(r) <- band
  
  return(r)
}
  
# Predict entire image
image_list <- map(c("rcc","gcc","bcc", "ndvi","bai"),compute_image)
predictors <- c(rast_refl,
                rast_focal,
                image_list[[1]], image_list[[2]], image_list[[3]], 
                image_list[[4]],image_list[[5]])

cat("Predicting raster...\n")
preds <- terra::predict(predictors, rf_fit)
cat("Writing raster...\n")
writeRaster(preds,
            filename = paste0("data/geodata/raster/burned_area/planet/",
                              "berelech_preds_top5_TD.tif"),
            overwrite = T
)
