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
  
  # Load raster data
  rast_refl <- rast(fn_raster)
  
  # Select training polygons within this raster
  sel_features <- features %>% 
    vect() %>% 
    project('EPSG:32655') %>% 
    terra::intersect(ext(rast_refl)) %>% 
    mutate(.,ID = 1:nrow(.)) %>% 
    mutate(.,Burned = factor(Burned,labels = c('unburned','burned')))
  
  # Extract raster data
  raster_vals <- terra::extract(rast_refl, sel_features)[,c(1,sel_bands+1)]
  
  training_data <- sel_features %>%
    as.data.frame() %>%
    select(ID, id, Burned) %>%
    full_join(raster_vals) %>%
    select(-ID)
  
}

# First cbh (10k per year and class)
training_sample <- training_data %>%
  group_by(id,Burned) %>%
  na.omit() %>%
  sample_n(10)

training_sample %>% group_by(id, Burned) %>% tally()

training_sample <- training_sample %>% 
  mutate(ndvi = (Red - NIR) / (Red + NIR),
         bai = 1 / ( (0.1 - Red)^2 + (0.06 - NIR)^2 ),
         rcc = Red / (Red + Blue + Green),
         bcc = Blue / (Red + Blue + Green),
         gcc = Green / (Red + Blue + Green),
  )

# Calculate separability
sep_stats <- bind_rows(
  separability(
    filter(training_sample, Burned == 0) %>% ungroup() %>% select(Red),
    filter(training_sample, Burned == 1) %>% ungroup() %>% select(Red)
  ) %>%
    select(-Red, -Red.1),
  separability(
    filter(training_sample, Burned == 0) %>% ungroup() %>% select(Green),
    filter(training_sample, Burned == 1) %>% ungroup() %>% select(Green)
  ) %>%
    select(-Green, -Green.1),
  separability(
    filter(training_sample, Burned == 0) %>% ungroup() %>% select(Blue),
    filter(training_sample, Burned == 1) %>% ungroup() %>% select(Blue)
  ) %>%
    select(-Blue, -Blue.1),
  separability(
    filter(training_sample, Burned == 0) %>% ungroup() %>% select(rcc),
    filter(training_sample, Burned == 1) %>% ungroup() %>% select(rcc)
  ) %>%
    select(-rcc, -rcc.1),
  separability(
    filter(training_sample, Burned == 0) %>% ungroup() %>% select(gcc),
    filter(training_sample, Burned == 1) %>% ungroup() %>% select(gcc)
  ) %>%
    select(-gcc, -gcc.1),
  separability(
    filter(training_sample, Burned == 0) %>% ungroup() %>% select(bcc),
    filter(training_sample, Burned == 1) %>% ungroup() %>% select(bcc)
  ) %>%
    select(-bcc, -bcc.1),
  separability(
    filter(training_sample, Burned == 0) %>% ungroup() %>% select(ndvi),
    filter(training_sample, Burned == 1) %>% ungroup() %>% select(ndvi)
  ) %>%
    select(-ndvi, -ndvi.1),
  separability(
    filter(training_sample, Burned == 0) %>% ungroup() %>% select(bai),
    filter(training_sample, Burned == 1) %>% ungroup() %>% select(bai)
  ) %>%
    select(-bai, -bai.1)
)


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

plot_list <- map(c("Red", "Green", "Blue", "rcc", "gcc", "bcc","ndvi", "bai"),
                    plot_separability)

plot_grid(plotlist = plot_list, 
          nrow = 2,
          labels = paste0(letters[1:8], ")")) %>%
  save_plot("figures/band_separability.png",
            .,
            nrow = 2,
            ncol = 4,
            base_asp = 1.2,
            bg = "white")

# Split into training and validation
training_sample$id <- 1:nrow(training_sample)
training <- training_sample %>%
  group_by(Burned) %>%
  slice_sample(prop = 0.8)
validation <- filter(training_sample, !(id %in% training$id)) 


# Train random forest model
rf_fit <- randomForest(Burned ~ Red + Green + Blue + rcc + bcc + gcc + ndvi + bai,
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
                image_list[[1]], image_list[[2]], image_list[[3]], 
                image_list[[4]],image_list[[5]])

cat("Predicting raster...\n")
preds <- terra::predict(predictors, rf_fit)
cat("Writing raster...\n")
writeRaster(preds,
            filename = paste0("data/geodata/raster/burned_area/planet/",
                              "berelech_preds_withBAI.tif"),
            overwrite = T
)
