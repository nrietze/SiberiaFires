# Script to crop the satellite imagery to the areas of interest
# Nils Rietze nils.rietze@uzh.ch 16 November 2023

library(terra)
library(tidyterra)
library(tidyverse)
library(sf)
library(parallel)
library(pbapply)
library(geostats)
library(jsonlite)

# define processing year
year <- 2020

# Load list of raster files
raster_path <- paste0('C:/data/8_planet/',year,'/original')

raster_files <- list.files(raster_path,
                           recursive = T,
                           pattern = '^composite.tif$',
                           full.names = T
)

# Load areas of interest
aois <- vect('./data/geodata/feature_layers/aoi_wv/aois.shp') %>%
  project('EPSG:32655') %>% 
  mutate(., id = 1:nrow(.))

# Function to crop satellite imagery from AOIs
crop_rasters <- function(feature_idx,features, fns_raster){
  feature <- features[feature_idx]
  
  cat("\nProcessing", feature$site, "\n")
  
  # Function to find the raster using the largest intersection area
  find_raster <- function(fns_raster, feature){
    intersect_areas <- lapply(raster_files, feature = feature, FUN = function(fn_raster, feature){
      
      # Load raster metadata files & geometries from json
      dir_path <- dirname(fn_raster)
      metadata <- vect(paste0(dir_path, '/composite_metadata.json') ) %>% 
        project('EPSG:32655')
      
      # find which aois intersect with the raster and extract intersection area
      intersection <- terra::intersect(feature, ext(metadata))
      area <- expanse(intersection)
    })
    
    # Check if there is any raster covering this AOI by looking for an intersect area
    any_non_empty_numeric <- any(sapply(intersect_areas, function(x) is.numeric(x) && length(x) > 0))
    
    # Report the filename of that raster
    if (any_non_empty_numeric) {
      idx <- which.max(intersect_areas)
      
      fn_raster_ok <- fns_raster[idx]
      
      return(fn_raster_ok)
    } else {
      return(NULL)
    }
    
  } 
  
  # Get raster covering this AOI
  fn_raster_ok <- find_raster(fns_raster,feature)
  
  cat(fn_raster_ok)
  
  # Break the function and return nothing if there is no raster for this AOI
  if (is.null(fn_raster_ok)) {
    cat("No raster covering this AOI was found.")
    return(NULL)
    stop()
  }
  
  dir_path <- dirname(fn_raster_ok)
  metadata <- vect(paste0(dir_path, '/composite_metadata.json') )
  instrument <- metadata$instrument
  
  # define band names per instrument
  if (instrument == 'PSB.SD'){
    band_names <- c('Coastal Blue','Blue','Green I','Green','Yellow','Red','Red Edge','NIR')
  } else {
    band_names <- c('Blue','Green','Red','NIR')
  }
  
  # Load raster data
  rast_refl <- rast(fn_raster_ok)
  rast_udm <- rast(gsub("(.*)(\\.tif)", "\\1_udm2\\2", fn_raster_ok))

  # mask reflectance with usable data mask provided by Planet
  rast_ok <- mask(rast_refl,rast_udm[[1]],maskvalue = 0)
  
  # crop raster to aoi
  out_raster <- terra::crop(rast_ok, feature,mask = TRUE)
  
  # rename new file
  new_dir_raster <- gsub("original", "cropped", sub("/\\d+$", "", dirname(fn_raster_ok)))
  new_fn_raster <- paste0(new_dir_raster,'/',feature$site,'_',year,'_composite.tif')
  new_fn_json <- paste0(new_dir_raster,'/',feature$site,'_',year,'_composite_metadata.json')
  
  # write cropped image to raster
  writeRaster(out_raster,
              new_fn_raster, 
              gdal=c("COMPRESS=NONE", "TFW=YES"),
              names = band_names,
              overwrite = TRUE)
  
  # write metadata to json
  writeVector(metadata, 
              new_fn_json,
              filetype = 'geojson',
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
  1:nrow(aois),
  features = aois,
  fns_raster = raster_files,
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

raster_path_cropped <- sprintf('C:/data/8_planet/%s/cropped/',year)
raster_files_cropped <- list.files(raster_path_cropped,
                              pattern = sprintf('*%s_composite.tif$',year),
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
