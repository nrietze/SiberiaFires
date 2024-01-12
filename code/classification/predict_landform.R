# Script to classify Arctic DEM in geomorphological units
# Nils Rietze nils.rietze@uzh.ch 16 November 2023

library(terra)
library(raster)
library(tidyverse)
library(tidyterra)
library(sf)
library(spatialEco)
library(caret)
library(randomForest)
library(cowplot)
library(scico)
library(parallel)
library(supercells)
library(pbapply)
library(rrrem)
library(future)

set.seed(10)

# Load data ----
PATH <- 'C:/data/4_geodata/arcticDEM/v3/UTM/'

dem_files <- list.files(PATH,pattern = 'aoi.*_dem_.*utm\\.tif$',full.names = T)
tpi_files <- list.files(PATH,pattern = 'aoi.*_tpi_.*\\.tif$',full.names = T)

pre_year <- 2019
planet_files_pre <- list.files(sprintf('C:/data/8_planet/%s/cropped/',pre_year),
                               pattern = '*composite.tif$',
                               full.names = T)

# load areas of interest
aois <- vect('./data/geodata/feature_layers/aoi_wv/aois.shp') %>% 
  project('EPSG:32655')

aoi_number <- 1
raster_index <- grep(paste0("aoi_", aoi_number, "_"), dem_files)

# load training polygons
training_polygons <- vect('data/geodata/feature_layers/landform/landform_training_polygons.gdb',
                          layer = 'aoi_1_training_polygons')
df_training_polygons <- as.data.frame(training_polygons)

# Load landforms by Huissteden 
gdb_path <- 'data/geodata/feature_layers/kytalyk_geomorphology.gdb'
landforms_kyt <- vect(gdb_path)
landforms_kyt_crop <- terra::intersect(aois[1],landforms_kyt)

# load DEM
dem <- rast(dem_files[raster_index])
names(dem) <- 'elevation'

# load focal DEM &TPI
dem_f99 <- rast(raster_files[2])
tpi_f199 <- rast(tpi_files[1])
names(tpi_f199) <- 'tpi_f199'

tpi_f99 <- rast(tpi_files[2])
names(tpi_f99) <- 'tpi_f99'

# compute slope
slope <- terrain(dem, "slope", unit = "degrees", neighbors = 8)

# A: Classify landform by thresholds ----
if(FALSE){
  # compute global statistics
  mean_dem <- global(dem,'mean', na.rm=TRUE)
  sd_dem <- global(dem,'sd', na.rm=TRUE)
  
  threshold_elevation <- 14
  threshold_std_deviations <- 2.5
  
  higher_than_13m <- dem > threshold_elevation
  higher_than_2sd <- dem > as.numeric(mean_dem + threshold_std_deviations * sd_dem)
  
  # Step 1a: If Yes, then "ridge"
  classified_dem <- rast(dem,vals = NA_integer_)
  classified_dem[higher_than_13m | higher_than_2sd] <- 2
  
  slope_threshold <- 5
  higher_than_slope_threshold <- slope > slope_threshold
  
  classified_dem[!is.na(classified_dem) & higher_than_slope_threshold] <- 1
  classified_dem <- focal(classified_dem,fun = 'modal',w = 5,na.policy = 'omit')
  
  plot(classified_dem)
  
  writeRaster(classified_dem,
              filename = paste0("data/geodata/raster/landform/",
                                "aoi_",aoi_number,"_landform.tif"),
              overwrite = T
  )

}



# B: Classify landform by histogram ----
if(FALSE){
  # Get DEM histogram 
  hist_values <- hist(dem, plot = FALSE)
  
  # number of classes
  num_classes <- 7
  
  ## Quantile breaks ----
  breaks <- quantile(hist_values$breaks, seq(0, 1, length.out = num_classes + 1))
  
  ## Equal breaks ----
  breaks <- seq(min(hist_values$breaks),max(hist_values$breaks), length.out = num_classes + 1)
  
  ## Geometric intervals breaks ----
  breaks <- exp(seq(log(max( 1, min(hist_values$breaks)) ), 
                    log(max(hist_values$breaks)), 
                    length.out = num_classes + 1))
  
  # Classify raster
  lower_bounds <- c(-Inf, head(breaks, -1))
  upper_bounds <- breaks
  labels <- 1:length(breaks)
  
  breaks_matrix <- cbind(lower_bounds, upper_bounds, labels)
  
  print(breaks_matrix)
  classified_dem <- classify(dem,breaks_matrix)
  
  # Plot classified raster
  ggplot() +
    geom_spatraster(data = classified_dem) +
    scale_fill_continuous('GnBu') +
    labs(title = "Classified DEM") +
    theme_cowplot()
  
  data_for_plot <- data.frame(values = spatSample(dem,size = 1e5))
  
  # Plot histogram with breaks as vertical lines
  ggplot(data_for_plot, aes(x = Band_1)) +
    geom_histogram(binwidth = .75, fill = "steelblue", color = "white", alpha = 0.7) +
    scale_x_continuous(breaks = breaks, labels = scales::number_format( accuracy = 0.01)) +
    geom_vline(xintercept = breaks, linetype = "dashed", color = "red") +
    labs(title = "DEM Histogram with Geometric Intervals", x = "Elevation (m)", y = "Frequency") +
    theme_cowplot()
}

# C: Classify landforms on segments ----

# Function to compute zonal statistics
get_combined_zonal_stats <- function(raster_list, polygons, summary_functions) {
  result_list <- lapply(raster_list, function(r) {
    sapply(summary_functions, function(fun) {
      zonal(r, polygons, fun = fun, na.rm = TRUE)
    })
  })
  
  combined_stats <- do.call(cbind, lapply(result_list, as.data.frame))
  return(combined_stats)
}

# define zonal statistics
summary_functions <- c('mean', 'sd', 'max', 'min')

# Set raster list
raster_list <- list(dem, slope,tpi_f99,tpi_f199)

## Generate supercells ----

# Save to shapefile
fname_slic <- sprintf('data/geodata/feature_layers/landform/aoi_%s_SLIC.shp',aoi_number)

if (!file.exists(fname_slic)){
  start.time <- Sys.time()
  print('No SLIC dataset for this AOI found. Computing new supercells...')
  
  plan(multisession, workers = 2)
  dem_slic <- supercells(dem, k = 1e3, compactness = .5,
                         chunks = TRUE, future = TRUE) 
  
  # convert to SpatVect
  v_slic <- vect(dem_slic)
  
  # compute slic feature area 
  dem_slic$area_km2 <- as.numeric(st_area(dem_slic)) / 1e6
  
  # get zonal stats in SLIC polygons
  combined_zonal_stats <- get_combined_zonal_stats(raster_list, v_slic, summary_functions)
  
  dem_slic <- cbind(dem_slic,combined_zonal_stats)
  
  # save to shapefile
  writeVector(vect(dem_slic),fname_slic,
              overwrite = T)
  
  print('done.')
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(time.taken)
  
} else {
  dem_slic <- st_read(fname_slic)
  v_slic <- vect(dem_slic)
}

# plot supercells
plot(v_slic,'supercells', lwd = 0.2)

# get zonal stats in training polygons
combined_zonal_stats <- get_combined_zonal_stats(raster_list, training_polygons, summary_functions)

df_training_polygons <- cbind(df_training_polygons,combined_zonal_stats)

## random forest for supercells ----

vars <- c("mean.elevation","sd.elevation","max.elevation","min.elevation",
          "mean.slope","sd.slope","max.slope","min.slope",
          # "mean.tpi_f99","sd.tpi_f99","max.tpi_f99","min.tpi_f99",
          "mean.tpi_f199","sd.tpi_f199","max.tpi_f199","min.tpi_f199")

# Prepare data
df_training_polygons$Classname <- factor(df_training_polygons$Classname)

train_indices <- createDataPartition(df_training_polygons$Classname, p = 0.8, list = FALSE)
training_data <- df_training_polygons[train_indices, ]
testing_data <- df_training_polygons[-train_indices, ]

rf_fit <- randomForest(formula(paste('Classname ~',paste(vars,collapse = "+"))),
                       data = select(df_training_polygons, c('Classname',vars))
                       )

top_5_gini <- varImp(rf_fit, scale = FALSE) %>% 
  top_n(5, Overall) %>%
  rownames()

print(top_5_gini)

fit_on_top5 <- T

if(fit_on_top5){
  rf_fit <- randomForest(formula(paste('Classname ~',paste(top_5_gini,collapse = "+"))),
                       data = select(df_training_polygons, c('Classname',top_5_gini))
)
}

# Predict on the testing set
test_preds <- predict(rf_fit, newdata = testing_data)

# Export confusion matrix
# fn_conf_matrix <- paste0("tables/rf_confusion_matrix_water_",site,"_",year,".txt")
# file_connection <- file(fn_conf_matrix)

confusionMatrix(data = test_preds, testing_data$Classname) %>%
  print() 

# Predict entire aoi
preds <- predict(rf_fit, newdata = st_drop_geometry(dem_slic),type = 'r')

v_slic$class <- preds
plot(v_slic,'class')

fname_pred <- sprintf('data/geodata/feature_layers/landform/aoi_%s_rf_classes.shp',aoi_number)
writeVector(v_slic,fname_pred,
            overwrite = T)


## k-means for supercells ----
df_sc <- st_drop_geometry(dem_slic['Band_1'])
km <- kmeans(df_sc,5)
dem_slic$cl <- km$cluster
plot(dem_slic['cl'])

