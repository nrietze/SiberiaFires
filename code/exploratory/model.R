library(terra)
library(tidyterra)
library(tidyverse)
library(scico)
library(cowplot)
library(GGally)
library(glcm)

# Load areas of interest
aois <- vect('./data/geodata/feature_layers/aoi_wv/aois.shp') %>%
  project('EPSG:32655') %>% 
  mutate(., id = 1:nrow(.))

aoi2 <- aois %>% filter(site == 'Kosukhino')

window_side_length <- 100 # in metres

# Load all data sets for model ----
# .................................

## DEM  ----
dem_path <- 'C:/data/4_geodata/arcticDEM/v3/UTM/'

dem_files <- list.files(dem_path,pattern = 'aoi.*_dem_.*utm\\.tif$',full.names = T)

aoi_number <- 2
raster_index <- grep(paste0("aoi_", aoi_number, "_"), dem_files)

dem <- rast(dem_files[raster_index]) %>% crop(aoi2)
names(dem) <- 'elevation'

ncells_resample <- window_side_length / 2 # e.g., if window is 100 m wide, a 2 m cell size gives 50 cells
dem_template <- aggregate(dem, fact = ncells_resample, fun = 'mean')
dem_100m <- resample(dem,dem_template, method = 'bilinear')

### Slope  ----
slope <- terrain(dem,'slope') %>% 
  resample(.,dem_template,method = 'bilinear') 

### Aspect ----
aspect <- terrain(dem,'aspect') %>% 
  resample(.,dem_template,method = 'bilinear') 

### Terrain roughness ----
tri <- terrain(dem,'TRI') %>% 
  resample(.,dem_template,method = 'bilinear')

### TPI ----
tpi_files <- list.files(dem_path,pattern = 'aoi.*_tpi_.*\\.tif$',full.names = T)
raster_index <- grep(paste0("aoi_", aoi_number, "_"), tpi_files)

tpi_500 <- rast(tpi_files[raster_index]) %>% 
  crop(aoi2) %>% 
  resample(.,dem_template,method = 'bilinear')
names(tpi_500) <- 'tpi_500'

## Burned area  ----
ba_path <- 'data/geodata/raster/burned_area/planet/Kosukhino_2020_burned_area.tif'

#### Percentage of burned or unburned in 100 m x 100 m windows ----
ba <- rast(ba_path) %>% crop(aoi2)

fraction_burned <- ifel(ba == 'burned', 1, 0) %>% 
  resample(.,dem_template,'sum') / (window_side_length / res(ba)[1])**2
names(fraction_burned) <- 'burned_fraction'

fraction_unburned <- ifel(ba == 'unburned', 1, 0) %>% 
  resample(.,dem_template,'sum') / (window_side_length / res(ba)[1])**2
names(fraction_burned) <- 'unburned_fraction'

plot(fraction_burned,
     range = c(0,1),
     col = scico(30,palette = 'bilbao',direction = -1),
     main = 'Burned fraction')

## LST ----
ls_path <- 'C:/data/9_landsat/2020/june/'

scale_factor <- 0.00341802 	
offset <- 149 - 273.15

lst <- rast(paste0(ls_path,'LC08_L2SP_115010_20200608_20200824_02_T1_ST_B10.TIF')) %>% 
  crop(aoi2) * scale_factor + offset
lst <- resample(lst,dem_template,'bilinear')
names(lst) <- 'LST'

plot(lst,col=scico(50,palette = 'lajolla'), main = "Surface temperature (°C)")

## NDVI ----
scale_factor <- 2.75e-05 	
offset <- -0.2

# Level 2
b4 <- rast(paste0(ls_path,'LC08_L2SP_115010_20200608_20200824_02_T1_SR_B4.TIF')) %>% 
  crop(aoi2) * scale_factor + offset
b5 <- rast(paste0(ls_path,'LC08_L2SP_115010_20200608_20200824_02_T1_SR_B5.TIF')) %>% 
  crop(aoi2) * scale_factor + offset

# Level 1
# b4 <- rast(paste0(ls_path,'LC08_L1TP_115010_20200608_20200824_02_T1_B4.TIF')) 
# b5 <- rast(paste0(ls_path,'LC08_L1TP_115010_20200608_20200824_02_T1_B5.TIF')) 

ndvi <- ((b5 - b4) / (b5 + b4)) %>% 
  resample(.,dem_template,'bilinear')
names(ndvi) <- 'NDVI'

plot(ndvi,
     range = c(-.5,.5),
     col = scico(30,palette = 'cork'),
     main = 'NDVI')

# GLCM of spectral bands
glcm_stats <- c("variance", "homogeneity","entropy")
glcm_nir <- glcm(as.array(b5[,,1]),statistics = glcm_stats)

glcm_nir_rast <- rast(glcm_nir, crs = crs(ndvi)) %>% 
  resample(.,dem_template,'bilinear')
names(glcm_nir_rast) <- glcm_stats

## water area ----
wa_path <- 'data/geodata/raster/water_area/planet/Kosukhino_2020_water_area.tif'
wa <- rast(wa_path) %>% resample(dem_template,'mode')

# Build model ----
# ................

# Data from random sample of raster values
extract_raster_values <- function(raster_list, random_points) {
  
  # Extract raster values at random points
  values_list <- lapply(raster_list, function(raster) {
    terra::extract(raster, random_points, ID = F)
  })
  
  # Combine the values into a data.frame
  df <- data.frame(do.call(cbind, values_list))
  
  return(df)
}

nsample <- 1e3

# Sample n burned pixels
random_points <- spatSample(dem_template,nsample, as.points = T,values = F)

stratified_random_points <- spatSample(resample(ba,dem_template,'mode'),
                                       nsample, method = 'stratified',
                                       as.points=TRUE)

raster_list <- list(fraction_burned,dem_100m,slope, aspect, tpi_500, lst,ndvi,glcm_nir_rast)
predictors <- c(fraction_burned,dem_100m,slope, aspect, tpi_500, lst,ndvi,glcm_nir_rast)

data <- extract_raster_values(raster_list, random_points) %>% 
  drop_na()

# All raster values
data <- data.frame(cbind(
  values(dem_100m),
  values(slope),
  values(aspect), 
  values(tpi_500),
  # values(tri),
  values(lst),
  values(ndvi),
  values(fraction_burned)
  )) %>% 
  drop_na()

gp <- ggpairs(data)
ggsave(gp,filename = 'figures/model/predictors_pairplot_aoi2.png',
       height = 8, width = 8)

brt_model <- gbm(unburned_fraction  ~ ., data = data,
                 distribution = 'gaussian',
                 n.trees = 1000, interaction.depth = 4)

# Summary and relative influence plots
summary(brt_model, plotit = TRUE)

# Partial dependence plots
plot(brt_model, i = 'LST', lwd = 2, main = "")
plot(brt_model, i = 'NDVI', lwd = 2, main = "")
plot(brt_model, i = 'tpi_500', lwd = 2, main = "")
plot(brt_model, i = 'slope', lwd = 2, main = "")
plot(brt_model, i = 'aspect', lwd = 2, main = "")

pred <- predict(predictors,brt_model)
values(pred) <- plogis(values(pred))
plot(pred,
     range = c(0,1),
     col = scico(30,palette = 'bilbao',direction = -1),
     main = 'Predicted burn fraction - BRT')

# Create a data frame with the extracted values
data_df <- data.frame(value1 = values(fraction_burned), 
                      value2 = values(pred)) %>% drop_na()

# Scatterplot
plot(data_df$class, data_df$lyr1, 
     pch = 16, col = "blue", main = "Scatterplot of Raster Values",
     xlab = "Raster 1 Values", ylab = "Raster 2 Values",ylim = c(0,1))

glm1 <- glm(unburned_fraction  ~ ., data = data)
summary(glm1)

pred <- predict(predictors,glm1)
plot(pred,
     range = c(0,1),
     col = scico(30,palette = 'bilbao',direction = -1),
     main = 'Predicted burn fraction - GLM')

rf_model <- randomForest::randomForest(unburned_fraction  ~ ., data = data)
rf_model

pred <- predict(predictors,rf_model)
plot(pred,
     range = c(0,1),
     col = scico(30,palette = 'bilbao',direction = -1),
     main = 'Predicted burn fraction - RF')
randomForest::varImpPlot(rf_model)

# Plots ----
# ..........

## Plot all predictors ----

# Define plotting parameters
parameters <- tibble(palettes = list(elevation = 'grayC',
                                     slope = "buda",
                                     aspect = "romaO",
                                     tpi_500 = "bam",
                                     LST = "lajolla",
                                     NDVI = "cork",
                                     unburned_fraction = "bilbao"),
                     labels = list(elevation = 'Elevation (m)',
                                   slope = 'Slope (°)',
                                   aspect = 'Aspect (°)',
                                   tpi_500 = expression(TPI[500~m]),
                                   LST = 'Land surface temperature (° C)',
                                   NDVI = 'NDVI (unitless)',
                                   unburned_fraction = NULL)
                     )

create_ggplot <- function(layer, parameters) {
  lyr_name <- names(layer)
  
  direction <- ifelse(lyr_name == 'unburned_fraction',-1,1)
  axis_labels <- ifelse(lyr_name == 'unburned_fraction','10 km',"")
  
  ggplot() +
    geom_spatraster(data = layer) +
    scale_fill_scico(palette = parameters$palettes[[lyr_name]], 
                     direction = direction) +
    theme_map(12) +
    theme(legend.position = "bottom",
          legend.box = "horizontal", 
          legend.justification = "center",
          legend.key.width = unit(0.05, "npc")) +
    guides(fill = guide_colourbar(title.position="bottom", title.hjust = 0.5)) +
    labs(title = names(layer),
         fill = parameters$labels[[lyr_name]])
}

ggplots <- lapply(names(predictors), function(layer) {
  create_ggplot(predictors[[layer]], parameters)
})

pred_plots <- plot_grid(plotlist = ggplots, ncol = 3, 
                        rel_heights = rep(1, length(ggplots)),
                        align = 'vh')

print(pred_plots)
ggsave2(pred_plots,filename = 'figures/model_predictors_aoi2.png',
        height = 16, width = 16)

## Scatterplot of observed vs. predicted burn fraction ----
data_df <- data.frame(obs = values(fraction_burned),
                      pred = values(pred))

colnames(data_df) <- c('obs','pred')

ggplot(data_df, aes(x = obs, y = pred)) +
  geom_point() +
  xlab('Burned fraction (Planet)')+
  xlab('Burned fraction (RF predicted)')

## Map of differences observed - predicted burn fractions ----
ggplot() +
  geom_spatraster(data=(fraction_burned - pred) ) +
  scale_fill_scico(palette = 'vik',
                   direction = -1,
                   guide = guide_colourbar(direction = "horizontal",
                                           title.position = 'top')) +
  labs(fill = 'Burned fraction residuals Observed - Predicted') +
  theme_minimal() +
  theme(legend.position = "bottom") 

## Plot burned area and burn fraction ----
red_colors <- brewer.pal(n = 3, name = "Reds")[-2]

p1 <- ggplot() +
  geom_spatraster(data = ba) +
  scale_fill_manual(values = red_colors) +
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.box = "horizontal", legend.justification = "center") +
  labs(fill = NULL,
       title = "Binary burn area (3 m x 3 m)")

p2 <- ggplot() +
  geom_spatraster(data = fraction_burned) +
  scale_fill_distiller(palette = "Reds", direction = 1,
                       guide = guide_colourbar(direction = "horizontal",
                                               title.position = 'bottom'),
                       limits = c(0,1)) +  
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.box = "horizontal", 
        legend.justification = "center",
        legend.key.width = unit(0.05, "npc")) +
  labs(fill = "burned fraction",
       title = "Burn fraction (100 m x 100 m)")

combined_plot <- plot_grid(p1, p2, ncol = 2,align = 'h')
print(combined_plot)
ggsave2(combined_plot,filename = 'figures/burn_aggregation_aoi2.png',
        bg = NULL,
        width = 12, height = 8)
