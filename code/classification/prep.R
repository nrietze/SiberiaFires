# Script to crop the satellite imagery to the areas of interest
# Nils Rietze nils.rietze@uzh.ch 16 November 2023

library(terra)
library(tidyverse)
library(sf)
library(parallel)
library(pbapply)
library(geostats)

# Load list of raster files
raster_path <- 'C:/data/8_planet/2020/original/'

raster_files_20 <- list.files(raster_path,
                               pattern = '*composite.tif$',
                               full.names = T
)

# Load areas of interest
aois <- read_sf('./data/geodata/feature_layers/aoi_wv/aois.shp') %>%
  mutate(., id = 1:nrow(.))

# Function to crop satellite imagery from AOIs
crop_rasters <- function(fn_raster, features, path_out){
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
clusterExport(cl, c("crop_rasters"), 
              envir=environment())

# Crop images
pblapply(
  raster_files_20,
  features = aois,
  FUN = crop_rasters,
  cl = 7)

# Stop cluster
stopCluster(cl)


# Compute focal rasters
window <- 5
target_stat <- 'sd'

# Function to compute focal raster statistics
compute_focal_raster <- function(fn_raster,band,window,f_stat){
  cat("\nProcessing", fn_raster, "\n")
  
  # Load raster
  rast_refl <- rast(fn_raster)
  
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
  
  # Compute band indices
  image_list <- map(c("rcc","gcc","bcc", "ndvi","bai"),compute_image)
  
  # add band indices to multiband raster
  rsrc <- rast(image_list)
  r_all <- rast(list(rast_refl,rsrc))
  
  # rename new file
  new_fn_raster <- gsub("composite", paste0("focal_", window, "_", f_stat), fn_raster)
  
  focal(r_all, window, 
        fun = f_stat, 
        na.rm = T, na.policy = "omit",
        filename = new_fn_raster,
        overwrite = T)
  
  return(NULL)
}

raster_path_cropped <- 'C:/data/8_planet/2020/cropped/'
raster_files_cropped <- list.files(raster_path_cropped,
                              pattern = '*composite.tif$',
                              full.names = T
)

# setup cluster
cl <- makeCluster(7)

clusterEvalQ(cl, c(library(terra),library(sf),library(tidyverse)))
clusterExport(cl, c("compute_focal_raster"), 
              envir=environment())

# Crop images
pblapply(
  raster_files_cropped,
  window = 5,
  f_stat = 'sd',
  FUN = compute_focal_raster,
  cl = 7)

# Stop cluster
stopCluster(cl)
