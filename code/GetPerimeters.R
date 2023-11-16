library(sf)
library(terra)
# library(raster)
library(tidyverse)

# Prepare paths ----
PATH_CCI <- 'C:/data/3_fire_data/burned_area/fire_cci/'

PATH_FIRMS <- 

year <- 2020

# Load raster ----
PATH_TIF <- paste(PATH_CCI,year, sep = '/')

# List all burned area files for that year
flist <- list.files(PATH_TIF,pattern = '*JD.tif', full.names = T)

i <- 1
ba <- rast(flist[i])

# Convert to binary ----
# Mask with unburned = 0, burned > 0
bin_ba <- ba > 0

# Vectorize ----

# Buffer ----

# Export to table & database ----