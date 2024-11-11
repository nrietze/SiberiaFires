# Script to crop the satellite imagery to the areas of interest
# Nils Rietze nils.rietze@uzh.ch 
# 16 November 2023

library(terra)
library(tidyterra)
library(tidyverse)
library(sf)
library(parallel)
library(pbapply)
library(geostats)
library(jsonlite)

# define processing year
year <- 2019

# Load list of raster files
raster_path <- paste0('C:/data/8_planet/',year,'/original')

raster_files <- list.files(raster_path,
                           recursive = T,
                           pattern = '^composite.tif$',
                           full.names = T
)

# Load areas of interest
aois <- vect('./data/geodata/feature_layers/aoi_wv/aois_analysis.geojson') %>%
  project('EPSG:32655') %>% 
  mutate(., id = 1:nrow(.))

# Function to crop satellite imagery from AOIs
crop_rasters <- function(aoi_name,aois, fns_raster){
  aoi <- aois[aois$site == aoi_name]
  
  FNAME_TEST <- paste0('C:/data/8_planet/',year,'/cropped')
  FLIST_EXIST <- list.files(FNAME_TEST,
                            recursive = T,
                            pattern = '*composite.tif$',
                            full.names = T
                            )
  
  # Check if files already exist
  if (any(grepl(aoi_name, FLIST_EXIST))){
    cat("Found processed data for this site. Skipping pre-processing.")
    return(NULL)
  }
  
  cat("\nProcessing", aoi$site, "\n")
  
  # Function to find the raster using the largest intersection area
  find_raster <- function(fns_raster, aoi){
    intersect_areas <- lapply(raster_files, feature = aoi, FUN = function(fn_raster, feature){
      
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
  fn_raster_ok <- find_raster(fns_raster,aoi)
  
  # Break the function and return nothing if there is no raster for this AOI
  if (is.null(fn_raster_ok)) {
    stop("No raster covering this AOI was found.")
  }
  
  dir_path <- dirname(fn_raster_ok)
  metadata <- vect(paste0(dir_path, '/composite_metadata.json') )
  instrument <- gsub("\\.","-",metadata$instrument)
  
  # define band names per instrument
  if (instrument == 'PSB-SD'){
    if (nlyr(rast(fn_raster_ok)) == 4){
      band_names <- c('Blue','Green','Red','NIR')
    } else {
      band_names <- c('Coastal Blue','Blue','Green I','Green','Yellow','Red','Red Edge','NIR')  
    }
  } else {
    band_names <- c('Blue','Green','Red','NIR')
  }
  
  # Load raster data
  rast_refl <- rast(fn_raster_ok)
  rast_udm <- rast(gsub("(.*)(\\.tif)", "\\1_udm2\\2", fn_raster_ok))

  # mask reflectance with usable data mask provided by Planet
  rast_ok <- mask(rast_refl,rast_udm[[1]],maskvalue = 0)
  
  # crop raster to aoi
  out_raster <- terra::crop(rast_ok, aoi,mask = TRUE)
  
  # rename new file
  new_dir_raster <- gsub("original", "cropped", dirname(fn_raster_ok))
  new_fn_raster <- paste0(new_dir_raster,'_',aoi$site,'_',instrument,'_composite.tif')
  cat(new_fn_raster)
  
  # write cropped image to raster
  writeRaster(out_raster,
              filename = new_fn_raster, 
              gdal=c("COMPRESS=NONE", "TFW=YES"),
              names = band_names,
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
  aois$site,
  aois = aois,
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
      r <- (rast_refl[["NIR"]] - rast_refl[["Red"]]) / (rast_refl[["NIR"]] + rast_refl[["Red"]])
    } else if(band == "ndwi"){
      r <- (rast_refl[["Green"]] - rast_refl[["NIR"]]) / (rast_refl[["Green"]] + rast_refl[["NIR"]])
    } else if(band == "bai"){
      r <- 1 / ( (0.1 - rast_refl[["Red"]])^2 + (0.06 - rast_refl[["NIR"]])^2 )
    }
    names(r) <- band
    
    return(r)
  }
  
  # Compute band indices
  image_list <- map(c("rcc","gcc","bcc", "ndvi","ndwi","bai"),compute_image)
  
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
                              pattern = sprintf('*_composite.tif$'),
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
