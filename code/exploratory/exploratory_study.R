rm(list = ls())

library(reticulate)
use_python("C:/Users/nils/AppData/Local/r-miniconda/envs/earthengine")

library(rgee)
ee_Initialize('ee-nrietze',gcs = FALSE)
ee_get_earthengine_path()

library(sf)
library(tidyverse)
library(terra)
library(raster)

library(scico)
scico_palette_show()

setwd('C:/Users/nils/OneDrive - Universität Zürich UZH/Dokumente 1/1_PhD/8_CHAPTER2/')

# LOAD DATA ----  
## Burned area from Descals et al. (2022) ----

# Sentinel 2 - Tile 17
burned_area_dsc_s2 <- raster('C:/data/3_fire_data/burned_area/descals_et_al_2022/raster/Arctic_BA_Sentinel2_2019-2020_tile-17_v1-1.tif') 

# Landsat 7 - 9 - Tile 17
burned_area_dsc_ls <- raster('C:/data/3_fire_data/burned_area/descals_et_al_2022/raster/landsat/Arctic_BA_Landsat78_2019-2020_tile-17_v1-1.tif') 
# burned_area_dsc_ls <- raster('C:/Users/nils/OneDrive - Universität Zürich UZH/Dokumente 1/1_PhD/8_CHAPTER2/data/geodata/raster/burned_area/ba_descals_landsat_2013-2020.tif') %>% 
#   crop(aoi)

# Set burned area product for analysis
burned_area_dsc <- burned_area_dsc_ls

# get extent of tile 17
aoi_tile_17 <- extent(burned_area_dsc)

# Set cavm extent as AOI
# aoi <- extent(cavm)

## Fire perimeteres from Talucci et al. (2022) ----
fire_perimeters <- read_sf('C:/data/3_fire_data/burned_area/siberia_talucci/data/SiberiaFires/SiberiaFires2001-2020_wgs.shp',
                           crs = 4326)

## CAVM raster ----
cavm_legend <- read.csv2('C:/data/6_vegetation/cavm/Raster CAVM legend.csv')

# subset in Siberia, original CRS
FNAME_CAVM_og <- './data/geodata/raster/cavm/cavm_descals_aoi_og_resolution.tif'
cavm_og <- raster(FNAME_CAVM_og)

# subset in Siberia buffered around fires (Lambert Azimuthal Equal Area, resampled to 30m)
FNAME_CAVM <- './data/geodata/raster/cavm/cavm_descals_subset_30m.tif'

# subset around tile 17 (UTM 55N, resampled to 30m)
FNAME_CAVM <- './data/geodata/raster/cavm/test_tile/raster_cavm_30m.tif'

# subset raster (Lambert Azimuthal Equal Area, resampled to 30m)
# FNAME_CAVM <- './data/geodata/raster/cavm/test_tile/raster_cavm_30m_az_eqar.tif'

# read CAVM raster & intersect with tile 17
cavm <- raster(FNAME_CAVM) %>% 
  intersect(burned_area_dsc_s2) %>% 
  as.factor()

# transform extent of og CAVM data to UTM 55N to crop it
ext_aoi_laea <- extent( projectExtent( cavm, crs(cavm_og) ) )
cavm_og <- crop(cavm_og, ext_aoi_laea)

# get final aoi
aoi <- extent(cavm)

# Save aoi as shapefile
p <- as(aoi, "SpatialPolygons")
crs(p) <- crs(cavm)

if (F){
  shapefile(p, './data/geodata/feature_layers/aoi_tile_17.shp')
}

# crop burned area to that aoi
burned_area_dsc <- crop(burned_area_dsc, aoi)

# cavm_classes <- unique(cavm)
cavm_classes <- c(1, 22, 23, 24, 33, 34, 42, 43, 91, 92, 99)

# load CAVM subbzone borders
cavm_zones <- read_sf('C:/data/6_vegetation/cavm/vegetation_zones.shp')

extent(burned_area_dsc) <- aoi
extent(cavm) <- aoi

# Function to derive burned area per CAVM class
get_burned_area <- function(year_value, burned_raster, land_cover_raster){
  stats_df <- zonal(burned_raster == year_value, land_cover_raster, "sum",na.rm = T) %>% 
    as.data.frame()
  # Convert number of pixels to ha
  stats_df$sum <- (stats_df$sum * 30 **2) / 1e4 
  return(stats_df$sum)
}

# integer values of each year (except unburned == 0)
year_values <- unique(values(burned_area_dsc)[values(burned_area_dsc) != 0]) 

# Create data frame and loop through each year and compute statistics
results_df <- data.frame(zone = cavm_classes,
                         data = sapply(year_values,get_burned_area,burned_area_dsc, cavm ))
colnames(results_df)[-1] <- sapply(year_values,function(val){paste0(20, val-10)})

merged_data <- merge(results_df, cavm_legend[,(1:2)], by.x = "zone", by.y = "Raster.code")

# map the colors to the data
merged_data <- merged_data %>% 
  mutate(color = cavm_colors[match(zone, cavm_classes)])

# Reshape the dataframe
merged_data_long <-  merged_data %>% 
  pivot_longer(cols = c("2019","2020"), 
               names_to = 'year',
               values_to = 'area') %>% 
  group_by(year) %>% 
  mutate(percent = area / sum(area) *100 )

merged_data_long$year <- as.numeric(merged_data_long$year)

cavm_colors <- rainbow(length(cavm_classes))

# Print the resulting data frame
print(merged_data[,1:4])

# PLOT BURNED AREA ----
OUT_PATH <- "./output/cavm_vs_burned_area/"

# Plot the second band (cavm class) with random colors
ggplot(merged_data, aes(x = Vegetation.Unit, y = sum, fill = Vegetation.Unit)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Vegetation unit", y = "Burned area (ha)", title = "Burned area per CAVM vegetation unit \n 2019") +
  theme_classic() + 
  scale_colour_brewer(palette = "Set1")

# Absolute burned area in tile 17
(p <- merged_data_long %>% 
    ggplot( aes(x=year, y=area, fill=Vegetation.Unit, text=Vegetation.Unit)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_viridis_d() +
    geom_text(data = subset(merged_data_long, area > 200),
              aes(label = Vegetation.Unit), position = position_stack(vjust = .5),
              size = 4, color = 'black') +  ggtitle("Burned area per CAVM vegetation unit (ha)") +
    theme_classic() +
    theme(legend.position="none"))
ggsave(paste0(OUT_PATH, "ba_cavm_absolute_tile_17_ls.png") )

# Relative burned area
(p <- merged_data_long %>%
    ggplot( aes(x=year, y=percent, fill=Vegetation.Unit, text=Vegetation.Unit)) +
    geom_bar(position="stack", stat="identity") +
    scale_fill_viridis_d() +
    geom_text(data = subset(merged_data_long, percent > 1),
              aes(label = Vegetation.Unit), position = position_stack(vjust = .5),
              size = 4, color = 'white') +
    ggtitle("Relative burned area per CAVM vegetation unit (%)") +
    theme_classic() +
    theme(legend.position="none"))
ggsave(paste0(OUT_PATH, "ba_cavm_relative_tile_17_ls.png") )

# Turn it interactive
p <- ggplotly(p, tooltip="text")
p

# LST PREFIRE in Earth engine ----

## Define areas & products ----
# area_burned <- ee$ImageCollection("MODIS/061/MCD64A1")
area_burned <- ee$ImageCollection("ESA/CCI/FireCCI/5_1")

# import fire perimeteres from Talucci et al. (2022)
fire_perimeters <- read_sf('C:/data/3_fire_data/burned_area/siberia_talucci/data/SiberiaFires/SiberiaFires2001-2020_wgs.shp',
                           crs = 4326)
unburned_perimeters <- read_sf('./data/geodata/feature_layers/talucci_perimeters_unburned_buffer_5km_join.shp',
                               crs = 4326)

# import burned area raster from Descals et al. (2022)
ee_ba_dsc <- ee$Image('projects/ee-nrietze/assets/burned_area_dsc_tile-17')

# Load CAVM data & legend
cavm <- ee$Image('projects/ee-nrietze/assets/raster_cavm_v1')
cavm_legend <- read.csv2('C:/data/6_vegetation/cavm/Raster CAVM legend.csv')

# Look at LST before fire ----
LST <- ee$ImageCollection("MODIS/061/MOD11A1") #daily
bands <- c('LST_Day_1km')

LST <- ee$ImageCollection("MODIS/061/MOD11A2") # 8-day composite
bands <- c('LST_Day_1km')

# Define the date range to filter
start_date <- ee$Date('2020-05-01')
end_date <- ee$Date('2020-09-30')

# Apply the conversion formula: °C = K - 273.15
convert_K_to_C <- function(image) {
  return(image$multiply(0.02)$subtract(273.15))
}

# Load tile 17 extent
aoi <- read_sf('./data/geodata/feature_layers/aoi_tile_17.shp', crs = 32655) %>% 
  st_transform(4326) 

aoi_ee <- sf_as_ee(aoi)

## Extract LST time series in Talucci perimeters ----

# filter out 2020 fires within AOI
fire_perimeters_20 <- fire_perimeters %>% 
  filter(FireYr == 2020) %>% 
  st_intersection(aoi)

unburned_perimeters_20 <- unburned_perimeters %>% 
  filter(FireYr == 2020) %>% 
  st_intersection(aoi)

fire_perimeters_ee <- sf_as_ee(fire_perimeters_20[,c('IDobj')])
unburned_perimeters_ee <- sf_as_ee(unburned_perimeters_20[,c('IDobj')])

## Extract LST data ----

GetImageryData <- function(imagery, features, bands, start_date, end_date){
  imagery_filtered <- imagery$
    select(bands)$
    filterDate(start_date, end_date)$
    filterBounds(features)
  
  image_dates <- ee_get_date_ic(imagery_filtered)$time_start
  
  imagery_filtered <- imagery_filtered$
    map(convert_K_to_C)
  
  # Extract LST from fire perimeters
  df_imagery_data <- ee_extract(
    x = imagery_filtered,
    y = features,
    scale = 1000,
    fun = ee$Reducer$mean(),
    sf = FALSE
  ) %>% 
    t() %>% # transpose
    as.data.frame() %>% 
    `colnames<-`(.[1, ]) %>% # assign fire perimeter id as column names
    .[-1, ] %>%  # remove the original id row
    mutate(date = ymd( sub("^X(\\d{4}_\\d{2}_\\d{2}).*$", "\\1", rownames(.)) ) )
  
  n_img <- nrow(df_imagery_data)
  rownames(df_imagery_data) <- 1:n_img
  
  # df_out <- data.frame(df_imagery_data, date = image_dates) # if imagery_dates are extracted separately and not from final dataframe index
  df_out <- df_imagery_data
  
  return(df_out)
}

df_fireperim_lst <- GetImageryData(LST, fire_perimeters_ee, bands, start_date, end_date)
df_unburned_lst <- GetImageryData(LST, unburned_perimeters_ee, bands, start_date, end_date)

df_fireperim_lst %>%
  pivot_longer(-date, names_to = "perim_id", values_to = "lst") %>%
  ggplot(aes(x = date, y = lst, color = perim_id)) +
  geom_line(alpha = 0.8, linewidth = 2) +
  xlab("Date") +
  ylab("Land surface temperature (°C)") +
  theme_minimal() 

# Combine the dataframes and add a 'site_type' column
df_fireperim_lst$site_type <- "burned"
df_unburned_lst$site_type <- "unburned"

# Combine the dataframes
combined_df <- rbind(df_fireperim_lst, df_unburned_lst)

# get day of year to cut off by fire start
combined_df$doy <- yday(combined_df$date)

# Extract the IDobj values from fire_perimeters_20
id_values <- fire_perimeters_20$IDobj

# Create a list of MinDay values corresponding to each column
nfeat <- nrow(fire_perimeters_20)
min_day_values <- fire_perimeters_20$MinDay[match(names(combined_df)[c(1:nfeat)], id_values)]

# Filter and extract the subset of data for each column using purrr::map2
subset_data_list <- map2(names(combined_df)[c(1:nfeat)], min_day_values, ~
                           combined_df %>%
                           filter(doy <= .y) %>%
                           dplyr::select(all_of(.x), doy, site_type ) # Adjust if needed
                         )

prefire_df <- bind_rows(subset_data_list) %>%
  pivot_longer(cols = -c(2,3), names_to = "perim_id", values_to = "lst")

# Boxplot of summer LST (with fire dates)
combined_df %>%
  pivot_longer(cols = starts_with("X"), names_to = "perim_id", values_to = "lst") %>%
  ggplot(aes(x = perim_id, y = lst, fill = site_type)) +
  geom_boxplot() +
  xlab("Perimeter ID") +
  ylab("Land Surface Temperature (°C)") +
  theme_minimal() +
  theme(legend.position = "top")

# Boxplot of summer LST (before fire dates)
prefire_df %>%
  ggplot(aes(x = perim_id, y = lst, fill = site_type)) +
  geom_boxplot() +
  geom_jitter(color="black", size=1, alpha=0.9) +
  xlab("Perimeter ID") +
  ylab("Land Surface Temperature (°C)") +
  theme_minimal() +
  theme(legend.position = "top")

### Plot LST over burn perimeters ----
image <- ee$Image( LST_filtered$toList(LST_filtered$size())$get(7) )

viz <- list(
  bands = bands,
  max = 35,
  min = -20,
  palette = rev(scico(30, palette = 'lajolla'))
)

Map$setCenter(145, 72, 4)
# Add all LST scenes to the map
Map$addLayers(eeObject = LST_filtered, 
              vis = viz,
              nmax = 2,
              name = LST_filtered$toBands()$bandNames()$getInfo()) +
Map$addLayer(fire_perimeters_ee$style(color = "006600", width = 2,
                                        fillColor = "00000000"), {}, "Fire perimeters")+
Map$addLegend(visParams = viz,name = 'LST (°C)')


data <- LST_filtered$map(function(image) {
  image$reduceRegions(
    collection = fire_perimeters_ee,
    reducer = ee$Reducer$mean(),
    scale = 1000
  )
  }
  )$flatten()
print(data$first()$getInfo())





# Select data from 2019 (==29) or 2020 (==30) and append CAVM band
ee_ba_single_year <- ee_ba_dsc$
  eq(30)

ee_ba_dsc_cavm <- ee_ba_single_year$updateMask(ee_ba_single_year)$
  addBands(cavm)

ee_stats <- ee_ba_dsc_cavm$reduceRegion(
  reducer = ee$Reducer$sum()$group(
    groupField = 1, # 1 because the CAVM classes are stored in the second band, indexed by a 1 (first is 0)
    groupName = "class"
  ),
  geometry = aoi,
  scale = 30,
  # maxPixels = 1e8
  bestEffort = T
)

# Print the resultant Dictionary.
stats <- ee_stats$getInfo()

# Extract the groups into separate lists
group_codes <- sapply(stats$groups, function(group) group$class)
group_stat <- sapply(stats$groups, function(group) group$sum)

# Create a data frame from the extracted lists
stats_df <- data.frame(class = group_codes, count = group_stat)

# Merge the data with vegetation data based on "code"
merged_data <- merge(stats_df, cavm_legend, by.x = "class", by.y = "Raster.code")

# Print the resulting data frame
print(merged_data[,1:4])

# Map the colors to the stats data frame
merged_data <- merged_data %>%
  mutate(color = cavm_colors[match(code, cavm_classes)])

# Plot the second band (cavm class) with random colors
ggplot(merged_data, aes(x = Vegetation.Unit, y = count, fill = Vegetation.Unit)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(x = "Vegetation Unit", y = "Count", title = "Burned Pixels per CAVM Vegetation Unit")

viz1 <- list(
  bands = 'b1',
  max = 1,
  min = 0,
  palette = list("FFFFFF", # unburned
                 "FF0000") # burned
)

viz2 <- list(
  bands = 'b1_1',
  max = 43,
  min = 0,
  palette = rev(scico(30, palette = 'hawaii'))
)

Map$setCenter(145, 72, 4)
Map$addLayer(ee_ba_dsc_cavm,viz2, name = 'CAVM') + 
Map$addLayer(ee_ba_dsc_cavm,viz1, name = 'Burned area') +
Map$addLayer(aoi,list(color = 'FF0000'),name = 'AOI')



ALAN_cum <- terra::rast("C:/Users/nils/OneDrive - Universität Zürich UZH/Dokumente 1/9_various/cengiz_paper/rasters/Arctic_arima_significant_clip.tif")
raster_as_ee(ALAN_cum)

# Define areas & products ----
area_burned <- ee$ImageCollection("MODIS/061/MCD64A1")

aoi <- ee$Geometry$Polygon(c(137.41580755062452,69.53198973831077,
                             150.95096380062452,69.53198973831077,
                             150.95096380062452,72.68361322664579,
                             137.41580755062452,72.68361322664579,
                             137.41580755062452,69.53198973831077))

# convert aoi to sf format for intersection with fire perimeters
aoi_sf <- ee_as_sf(aoi)

# import fire perimeteres from Talucci et al. (2022)
fire_perimeters <- read_sf('C:/data/3_fire_data/perimeters/siberia_talucci/data/SiberiaFires/SiberiaFires2001-2020_wgs.shp')

# filter out 2020 fires > 10 kha
fire_perimeters_20 <- fire_perimeters %>% 
  filter(FireYr == 2020 & SizeHa >= 10000) %>% 
  st_intersection(aoi_sf)

# Extract fires within AOI and store as ee.FeatureCollection
fire_perimeters_ee <- sf_as_ee(fire_perimeters_20[,c('SizeHa')])

# Filter area burned in AOI between dates selected
bands <- c('FirstDay')
ee_ba_20 <- area_burned$
  select(bands)$
  filterBounds(aoi)$
  filterDate("2020-05-01", "2020-10-31")

viz <- list(
  bands = bands,
  max = 366,
  min = 0,
  palette = scico(30, palette = 'lajolla')
)

# Show map https://www.rdocumentation.org/packages/rgee/versions/1.1.5/topics/Map
# Map$addLayers(eeObject = ee_ba_20,
#              vis = viz,
#              name = ee_ba_20$toBands()$bandNames()$getInfo()) +
Map$addLayer(eeObject = ee_ba_dsc,
             name = "Area burned (Copernicus)") +
  Map$addLayer(fire_perimeters_ee$style(color = "006600", width = 2,
                                        fillColor = "00000000"), {}, "Fire perimeters")+
  Map$addLegend(visParams = viz,name = 'Burn dates')
Map$centerObject(aoi)

