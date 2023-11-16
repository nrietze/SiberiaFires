library(terra)
library(tidyverse)
library(sf)
library(caret)
library(randomForest)
library(ggplot2)
library(cowplot)
library(dplyr)
library(purrr)
library(parallel)
library(pbapply)

# Load areas of interest
aois <- read_sf('./data/geodata/feature_layers/aoi_wv/aois.shp') %>%
  mutate(., id = 1:nrow(.))

# Load list of raster files


# Function to compute band ratio vegetation index
VI <- function(img, k, i) {
  bk <- img[[k]]
  bi <- img[[i]]
  vi <- (bk - bi) / (bk + bi)
  return(vi)
}

# Plot Planet SuperDove RGB 
RGB(r1) <- c(6,3,2)
plot(r1)


# NDVI, NIR = 8, Red = 5
ndvi <- VI(r1, 8, 5)

plot(ndvi, col = scico(10, palette = 'cork'), main = 'Planet NDVI')
