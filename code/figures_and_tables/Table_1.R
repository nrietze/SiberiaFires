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
    # aoi <- aois %>% filter(site == aoi_name)
    # poly_mask <- terra::intersect(poly_mask,aoi)

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
    sprintf('data/geodata/raster/water_area/%s_Landsat_mask.tif',aoi_name)
  )
  
  # load PlanetScope burned area
  ba <- rast(
    sprintf('data/geodata/raster/burned_area/planet/%s_burned_area_top5TD.tif',aoi_name)
  ) 
  
  wa_3m <- resample(wa,ba)
  
  # mask water areas in burned area map
  ba <- mask(ba,wa_3m,maskvalues = 2,updatevalue = NA) %>% 
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
  
  ba_comp_bin <- ifel(ba_comp == burn_val, 1, 0) %>%  #convert values to 0 & 1
    mask(ba_comp,maskvalues = NA, updatevalue = 0) %>% # convert NA values to 0 (unburned), important for GABAM
    mask(crop(wa,bp),maskvalues = 2, updatevalue = NA) %>%  # mask out water areas
    mask(bp) %>%
    as.factor()
  
  levels(ba_comp_bin) <- data.frame(id=0:1, class= c('unburned','burned'))

  if (has_mask){
    ba_comp_bin <- ba_comp_bin %>% 
      mask(poly_mask,inverse = T)
    
    ba <- ba %>% 
      mask(poly_mask,inverse = T)
  }
  
  # get burned area in Landsat product
  A_comp <- ifel(ba_comp_bin == "burned",TRUE,NA) %>% 
    expanse(unit="km",transform = FALSE)
  
  # get burned area in PlanetScope product
  A_PS_bd <- ifel(ba == "burned",TRUE,NA) %>% 
    expanse(unit="km",transform = FALSE)
  
  # get unburned area in PlanetScope product
  A_PS_unbd <- ifel(ba == "unburned",TRUE,NA) %>% 
    expanse(unit="km",transform = FALSE)
  
  return(data.frame(site = aoi_name,A_comp = A_comp$area,
                    A_PS_unbd = A_PS_unbd$area,A_PS_bd = A_PS_bd$area))
}

A_desc <- do.call(rbind,lapply(aois$site,GetBurnedArea))
A_GABAM <- do.call(rbind,lapply(aois$site,GetBurnedArea,product = "GABAM"))

A <- A_desc %>% 
  mutate(prop_desc = A_PS_bd / A_comp,
         pct_unburned = A_PS_unbd / (A_PS_bd + A_PS_unbd) *100 ) %>% 
  rename(A_desc = A_comp) %>% 
  as_tibble() %>% 
  left_join(A_GABAM %>% 
              mutate(prop_gabam = A_PS_bd / A_comp) %>% 
              rename(A_gabam = A_comp) %>% 
              as_tibble()) %>% 
  relocate(c(A_PS_unbd,A_PS_bd,pct_unburned),.before = A_desc)

# Format and export as formatted table ----
A %>% 
  mutate(site = str_replace_all(site, "(?<=[a-z])(?=[A-Z])", " ")) %>% 
  gt() %>% 
  tab_spanner(
    label = html("PlanetScope <br>(this study)"),
    columns  = c(A_PS_bd,A_PS_unbd,pct_unburned),
    id = "PS"
  ) %>% tab_spanner(
    label = html("Landsat <br>(Descals et al., 2022)"),
    columns  = c(A_desc,prop_desc),
    id = "DESC"
  ) %>% tab_spanner(
    label = html("Landsat <br>(Wei et al., 2021)"),
    columns  = c(A_gabam,prop_gabam),
    id = "GABAM"
  ) %>% 
  cols_label(
    site = "Fire Scar",
    A_PS_unbd = html("unburned area<br>km<sup>2</sup>"),
    A_PS_bd = html("burned area<br>km<sup>2</sup>"),
    pct_unburned = html("Percentage unburned <br>%"),
    A_desc = html("burned area<br>km<sup>2</sup>"),
    prop_desc = html("Ratio <br> PS:LS8"),
    A_gabam = html("burned area<br>km<sup>2</sup>"),
    prop_gabam = html("Ratio <br> PS:LS8")
  ) %>% 
  fmt_number(decimals = 1) %>% 
  tab_style(
    style = cell_text(v_align = "top", weight = 'bold'),
    locations = cells_column_spanners()) %>%
  tab_style(
    style = cell_text(v_align = "top", weight = 'bold'),
    locations = cells_column_labels()) %>% 
  gtsave(filename = "tables/Table1.html")
