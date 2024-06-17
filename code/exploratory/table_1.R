# Script to compute and gather area statistics of rasters
# Nils Rietze: nils.rietze@uzh.ch 
# 08 May 2024

library(terra)
library(tidyterra)
library(tidyverse)
library(gt)

# Load study regions ----
aois <- vect('./data/geodata/feature_layers/aoi_wv/aois_analysis.geojson') %>%
  project('EPSG:32655') %>% 
  mutate(., id = 1:nrow(.))

poly_mask <- vect('data/geodata/feature_layers/planet_masks.shp') 

GetBurnedArea <- function(aoi_name,product = "descals"){
  cat(sprintf("processing %s ... \n", aoi_name))
  
  # load burn perimeter
  bp <- vect(
    sprintf('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_%s.shp',aoi_name)
    )
    
  if (aoi_name %in% c("Berelech","LargeScarCenter")){
    aoi <- aois %>% filter(site == aoi_name)
    poly_mask <- terra::intersect(poly_mask,aoi)
    
    has_mask <- TRUE
    
  } else{has_mask <- FALSE}
  
  # load burned fractions
  fraction_burned <- rast(
    sprintf('data/geodata/raster/predictors/%s_predictors_30m.tif',aoi_name)
    ) %>% 
    select(burned_fraction) %>% 
    crop(bp)
  
  # load water mask
  wa <- rast(
    sprintf('data/geodata/raster/water_area/planet/%s_water_area_top5TD.tif',aoi_name)
    )
  wa <- ifel(wa == 'water',wa,NA)
  
  # load PlanetScope burned area
  ba <- rast(
    sprintf('data/geodata/raster/burned_area/planet/%s_burned_area_top5TD.tif',aoi_name)
  ) %>% 
    mask(wa,maskvalues = 2,updatevalue = NA) %>% 
    mask(bp)
  
  # load comparison Landsat burned areas
  if (product == "descals"){
    fn_comp <- 'data/geodata/raster/burned_area/ba_descals_landsat_2020_utm_shifted.tif'
    # fn_desc <- 'data/geodata/raster/burned_area/ba_descals_sentinel_2020_utm_shifted.tif'
    burn_val <- 30
    
  } else if (product == "GABAM"){
    fn_comp <- 'data/geodata/raster/burned_area/N75E145_burn_class_UTM_55N.tif'
    burn_val <- 255
  }
  
  ba_comp <- rast(fn_comp) %>% 
    crop(bp) %>% 
    mask(fraction_burned,maskvalues = NA) %>% 
    rename(class = 1) 
  
  ba_comp_bin <- ifel(ba_comp == burn_val, 1, 0) %>%  
    as.factor()
  levels(ba_comp_bin) <- data.frame(id=0:1, class= c('unburned','burned'))

  if (has_mask){
    ba_comp_bin <- ba_comp_bin %>% 
      mask(poly_mask,inverse = T)
    
    ba <- ba %>% 
      mask(poly_mask,inverse = T)
  }
  
  # get burned area in Landsat product
  A_comp <- expanse(ifel(ba_comp_bin == "burned",TRUE,NA),
                    unit="km",
                    transform = FALSE)
  
  # get burned area in PlanetScope product
  A_PS <- expanse(ifel(ba == "burned",TRUE,NA),
                  unit="km",
                  transform = FALSE)
  
  return(data.frame(site = aoi_name,A_comp = A_comp$area,A_PS = A_PS$area))
}

A_desc <- do.call(rbind,lapply(aois$site,GetBurnedArea))
A_GABAM <- do.call(rbind,lapply(aois$site,GetBurnedArea,product = "GABAM"))

A <- A_desc %>% 
  mutate(prop_desc = A_PS / A_comp) %>% 
  rename(A_desc = A_comp) %>% 
  as_tibble() %>% 
  left_join(A_GABAM %>% 
              mutate(prop_gabam = A_PS / A_comp) %>% 
              rename(A_gabam = A_comp) %>% 
              as_tibble()) %>% 
  relocate(A_PS,.before = A_desc)

# Format and export as formatted table ----
A %>% 
  mutate(site = str_replace_all(site, "(?<=[a-z])(?=[A-Z])", " ")) %>% 
  gt() %>% 
  cols_label(
    site = "Site",
    A_PS = html("PlanetScope,<br>km<sup>2</sup>"),
    A_desc = html("Landsat (Descals et al., 2022),<br>km<sup>2</sup>"),
    prop_desc = html("Ratio <br> PlanetScope:Landsat"),
    A_gabam = html("Landsat (Wei et al., 2021),<br>km<sup>2</sup>"),
    prop_gabam = html("Ratio <br> PlanetScope:Landsat")
  ) %>% 
  fmt_number(decimals = 1) %>% 
  tab_style(
    style = cell_text(v_align = "top", weight = 'bold'),
    locations = cells_column_labels()) %>% 
  gtsave(filename = "tables/Table1.html")
