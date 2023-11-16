# Script to crop the satellite imagery to the areas of interest
# Nils Rietze nils.rietze@uzh.ch 16 November 2023

library(terra)
library(tidyverse)
library(sf)
library(parallel)
library(pbapply)

# Load list of raster files
DIR_PLT_20 <- 'C:/data/8_planet/2020/original/'

raster_files_20 <- list.files(DIR_PLT_20,
                               pattern = '*composite.tif$',
                               full.names = T
)

# Load areas of interest
aois <- read_sf('./data/geodata/feature_layers/aoi_wv/aois.shp') %>%
  mutate(., id = 1:nrow(.))

# Function to extract satellite imagery from AOIs
extract_values <- function(fn_raster, features, path_out){
  cat("\nProcessing", fn_raster, "\n")
  
  # Load raster data
  rast_refl <- rast(fn_raster)
  rast_udm <- rast(gsub("(.*)(\\.tif)", "\\1_udm2\\2", fn_raster))

  # mask reflectance with usable data mask provided by Planet
  rast_ok <- mask(rast_refl,rast_udm[[1]])
  
  # Check which aoi is covered in this raster
  features_utm <- vect(features) %>% 
    project('EPSG:32655')
  
  # find which aois intersect with the raster and extract intersection area
  intersection <- terra::intersect(features_utm, ext(rast_ok))
  intersect_areas <- expanse(intersection)
  
  # select the larger intersection area to use as crop extent
  aoi <- features[intersection[which.max(intersect_areas),]$id,] %>% 
    vect() %>% 
    project('EPSG:32655')
  
  # crop raster to aoi
  out_raster <- terra::crop(rast_ok, aoi,mask = TRUE)
  
  # rename new file
  new_fn_raster <- gsub("(.*\\/\\d{4}-\\d{2}-\\d{2}_)(.*)(_composite\\.tif)", 
                        paste0("\\1", aoi$site, "\\3"), fn_raster)
  new_fn_raster <- gsub("(.*\\/\\d{4})(\\/original\\/)(.*)", "\\1/cropped/\\3", new_fn_raster)
  
  # write cropped image to raster
  writeRaster(out_raster,new_fn_raster, 
              gdal=c("COMPRESS=NONE", "TFW=YES"),
              names = c('Coastal Blue','Blue','Green I','Green','Yellow','Red','Red Edge','NIR'),
              overwrite = TRUE)
  
  return(NULL)
}

# setup cluster
cl <- makeCluster(7)

clusterEvalQ(cl, c(library(terra),library(sf),library(tidyverse)))
clusterExport(cl, c("extract_values"), 
              envir=environment())

# Correct images
pblapply(
  raster_files_20,
  features = aois,
  FUN = extract_values,
  cl = 7)

# Stop cluster
stopCluster(cl)
