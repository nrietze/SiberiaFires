library(terra)
library(tidyterra)
library(tidyverse)
library(tidybayes)
library(scico)
library(RColorBrewer)
library(cowplot)
library(GGally)
library(glcm)
library(gbm)
library(caret)
library(car)
library(brms)
library(DHARMa)
library(sjPlot)

# Load areas of interest
aois <- vect('./data/geodata/feature_layers/aoi_wv/aois.shp') %>%
  project('EPSG:32655') %>% 
  mutate(., id = 1:nrow(.))

aoi2 <- aois %>% filter(site == 'Kosukhino')
aoi4 <- aois %>% filter(site == 'Berelech')

window_side_length <- 30 # in metres

# Create template grid from Landsat data
ls_path <- 'C:/data/9_landsat/2020/june/'
lst_og <- rast(paste0(ls_path,'LC08_L2SP_115010_20200608_20200824_02_T1_ST_B10.TIF')) %>% 
    crop(aoi2)

if (window_side_length == 30 ){ # Use landsat raster as template if grid cells are 30 x 30 m
  raster_grid_template <- rast(ext(lst_og), resolution=res(lst_og)) 
  crs(raster_grid_template) <- crs(lst_og)
} else {
  ncells_resample <- window_side_length / 2 # e.g., if window is 100 m wide, a 2 m cell size gives 50 cells
  raster_grid_template <- aggregate(lst_og, fact = ncells_resample, fun = 'mean')
}

# Load burn perimeter
bp <- vect('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_aoi2.shp') %>% 
  crop(aoi2)

# Load all data sets for model ----
# .................................

## DEM  ----
dem_path <- 'C:/data/4_geodata/arcticDEM/v3/UTM/'

dem_files <- list.files(dem_path,pattern = 'aoi.*_dem_.*utm\\.tif$',full.names = T)

aoi_number <- 2
raster_index <- grep(paste0("aoi_", aoi_number, "_"), dem_files)

dem_og <- rast(dem_files[raster_index]) %>% crop(aoi2)

dem <- resample(dem_og,raster_grid_template, method = 'cubicspline')
names(dem) <- 'elevation'

### Slope  ----
slope <- terrain(dem_og,'slope') %>% 
  resample(.,raster_grid_template,method = 'cubicspline') 

### Aspect ----
aspect <- terrain(dem_og,'aspect') %>% 
  resample(.,raster_grid_template,method = 'cubicspline') 

northness <- cos(aspect * pi/180)
names(northness) <- 'northness'

eastness <- sin(aspect * pi/180)
names(eastness) <- 'eastness'

### TPI ----
tpi_files <- list.files(dem_path,pattern = 'aoi.*_tpi_.*\\.tif$',full.names = T)
raster_index <- grep(paste0("aoi_", aoi_number, "_"), tpi_files)

tpi_500 <- rast(tpi_files[raster_index]) %>% 
  crop(aoi2) %>% 
  resample(.,raster_grid_template,method = 'cubicspline')
names(tpi_500) <- 'tpi_500'

## Burned area  ----
ba_path <- 'data/geodata/raster/burned_area/planet/Kosukhino_2020_burned_area.tif'

#### Percentage of burned or unburned in 100 m x 100 m windows ----
ba <- rast(ba_path) %>% crop(aoi2)

fraction_burned <- ifel(ba == 'burned', 1, 0) %>% 
  resample(.,raster_grid_template,'sum') / (window_side_length / res(ba)[1])**2
names(fraction_burned) <- 'burned_fraction'

fraction_unburned <- ifel(ba == 'unburned', 1, 0) %>% 
  resample(.,raster_grid_template,'sum') #/ (window_side_length / res(ba)[1])**2
names(fraction_unburned) <- 'unburned_fraction'

## LST ----
ls_path <- 'C:/data/9_landsat/2020/june/'
# ls_path <- 'C:/data/9_landsat/2018/july/LC08_L2SP_117010_20180719_20200831_02_T1_ST_B10.TIF'

scale_factor <- 0.00341802 	
offset <- 149 - 273.15

lst <- rast(paste0(ls_path,'LC08_L2SP_115010_20200608_20200824_02_T1_ST_B10.TIF')) %>% 
  crop(aoi2) * scale_factor + offset
names(lst) <- 'LST'

if (window_side_length != 30){
  lst <- resample(lst_og,raster_grid_template,'bilinear')
}

## NDVI ----

# Landsat
scale_factor <- 2.75e-05 	
offset <- -0.2

# Level 1
# b4 <- rast(paste0(ls_path,'LC08_L1TP_115010_20200608_20200824_02_T1_B4.TIF')) 
# b5 <- rast(paste0(ls_path,'LC08_L1TP_115010_20200608_20200824_02_T1_B5.TIF')) 

# Level 2
red <- rast(paste0(ls_path,'LC08_L2SP_115010_20200608_20200824_02_T1_SR_B4.TIF')) %>% 
  crop(aoi2) 
nir <- rast(paste0(ls_path,'LC08_L2SP_115010_20200608_20200824_02_T1_SR_B5.TIF')) %>% 
  crop(aoi2) 

ndvi <- ((nir - red) / (nir + red)) 
names(ndvi) <- 'NDVI'

if (window_side_length != res(ndvi)[1]){
  ndvi <-  resample(ndvi,raster_grid_template,'bilinear')
}

# Planet
planet_path <- 'C:/data/8_planet/2019/cropped/20190621_Kosukhino_PS2-SD_composite.tif'
# planet_path <- 'C:/data/8_planet/2020/cropped/20200615_Kosukhino_PS2-SD_composite.tif'

planet_msp <- rast(planet_path) %>% 
  crop(aoi2) 

get_ndvi_heterogeneity <- function(red, nir){
  ndvi <- ((nir - red) / (nir + red)) 
  
  # Compute patchiness metric before resampling
  ndvi_sd <- terra::aggregate(ndvi, fact = 10,fun = 'sd') %>% 
    resample(.,raster_grid_template,'bilinear')
  names(ndvi_sd) <- 'NDVI_sd'
  
  return(ndvi_sd)
}

ndvi_sd <- get_ndvi_heterogeneity(planet_msp$Red,planet_msp$NIR)

# GLCM of spectral bands
# glcm_stats <- c("variance", "homogeneity","entropy")
# 
# min_nir <- as.integer(global(nir,'min',na.rm = T))
# max_nir <- as.integer(global(nir,'max',na.rm = T))
# 
# glcm_nir <- glcm(as.array(nir)[,,1],statistics = glcm_stats,
#                  min_x = min_nir,max_x = max_nir,na_opt = 'ignore')
# 
# glcm_nir_rast <- rast(glcm_nir, crs = crs(ndvi), extent = ext(b5)) %>%
#   resample(.,raster_grid_template,'average')
# names(glcm_nir_rast) <- glcm_stats

## water area ----
wa_path <- 'data/geodata/raster/water_area/planet/Kosukhino_2020_water_area.tif'
wa <- rast(wa_path) %>% resample(raster_grid_template,'mode')

## gather all predictors and mask water areas ----
predictors_unscaled <- c(dem,slope, aspect,northness,eastness, tpi_500,
                         lst,ndvi,ndvi_sd) 
predictors <- c(scale(predictors),fraction_burned) %>% 
  mask(wa,maskvalue = 2, updatevalue = NA)

# Build model ----
# ................

# Get data from random sample of raster values
nsample <- 2e3

# Sample n pixels and mask with burn perimeter
random_points <- spatSample(raster_grid_template,nsample, as.points = T,values = F) %>% 
  mask(bp)

# Extract data at random points 
data <- terra::extract(predictors, random_points,ID = F) %>% 
  drop_na()

# Define predictor variables for the model
mod_vars <- c('elevation','slope','northness','eastness','tpi_500', 
              'LST','NDVI','NDVI_sd')

mod_formula <- formula(paste('burned_fraction ~',paste(mod_vars, collapse = '+')))
# mod_formula_zi <- formula(paste('zoi ~',paste(mod_vars, collapse = '+')))
# mod_formula_coi <- formula(paste('coi ~',paste(mod_vars, collapse = '+')))

# Check for multicollinearity using a linear model
lm_mod <- lm(mod_formula, data = data)
performance::multicollinearity(lm_mod)

# Pair plot of predictors
(gp <- ggpairs(data[c('burned_fraction',mod_vars)]))
ggsave(gp,filename = 'figures/model/scaled_predictors_pairplot_aoi2_v3.png',
       height = 8, width = 8)

## Zero-One inflated beta regression ----
zoi_beta_model <- brm(
  bf(mod_formula),
  data = data,
  family = zero_one_inflated_beta(),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234,
  file = paste0("zoi_beta_model_",today())
)

summary(zoi_beta_model)

# Compare models with and without NDVI:
zoi_beta_model_nondvi <- readRDS('zoi_beta_model_excl_ndvi_30m.rds')
zoi_beta_model_withndvi <- readRDS('zoi_beta_model_incl_ndvi_30m.rds')

# Use WAIC criterion:
waic(zoi_beta_model_nondvi,zoi_beta_model_withndvi,compare = T)

# elpd = expected log pointwise predictive density

# Diagnostic plots:
plot(zoi_beta_model)
pp_check(zoi_beta_model)

residuals <- createDHARMa(simulatedResponse = t(posterior_predict(zoi_beta_model, ndraws = 100)),
                          observedResponse = data$burned_fraction,    # response variable
                          fittedPredictedResponse = apply(t(posterior_epred(zoi_beta_model, ndraws = 100)), 1, mean),
                          integerResponse = TRUE)

plot(residuals) 

# Check individual predictors: 
plotResiduals(residuals, form = data$LST)

# extract conditional effects
ce <- conditional_effects(zoi_beta_model)
plot(ce,points = T)

# Extract posterior estimates
dat_for_plot <- zoi_beta_model %>%
  # The `b_.*` tells it to take any estimated parameter starting with 'b_'
  gather_draws(`b_.*`, regex = TRUE) %>%
  # You might want to keep the intercept, depends on the analysis: 
  filter(!`.variable` %in% c("b_Intercept")) %>%
  
  # Reorder the levels:
  mutate(`.variable` = factor(`.variable`,
                              levels = get_variables(zoi_beta_model)[c(2:(length(mod_vars)+1))])) %>% 
  mutate(label = tools::toTitleCase(sub("^b_", "", .variable)) )

# Plot the posterior estimates, more info here: https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html 
plot <- ggplot(data = dat_for_plot, aes(y = label, x = .value, fill = after_stat(x < 0))) +
  stat_halfeye(size = 0.5) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_y_discrete(labels = mod_labs) +
  scale_fill_manual(values = coef_cols) +
  theme_minimal_hgrid() + 
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(x = "Posterior estimates", y = ""); plot

ggsave(plot,filename = sprintf('figures/model/zoib_model_aoi2_%s.png',today()),bg = 'white')

# Predict burn fraction in entire map
preds_zoib <- predict(zoi_beta_model,predictors)

preds_zoib_rast_df <- as.data.frame(predictors$elevation,xy = T) 
preds_zoib_rast_df$burned_fraction_pred <- preds_zoib[,1]

preds_zoib_rast <- rast(preds_zoib_rast_df,type = 'xyz',
           crs = "EPSG:32655",extent = ext(predictors)) %>% 
  subset('burned_fraction_pred')

plot(preds_zoib_rast,col = scico(30, palette = 'bilbao',direction = -1))
plot(mask(preds_zoib_rast,bp),col = scico(30, palette = 'bilbao',direction = -1))

# scatterplot obs. vs. pred
data_df <- data.frame(value1 = values(mask(fraction_burned,bp)), 
                      value2 = values(mask(preds_zoib_rast,bp))) %>% 
  drop_na()

# Scatterplot
ggplot() + 
  geom_point(data = data_df, aes (x = data_df[,1], y = data_df[,2])) +
  labs(x = "observed burn fraction",y = "predicted burn fraction") +
  lims(x = c(0,1),y = c(0,1)) +
  theme_minimal()

## Boosted regression tree ----
brt_model <- gbm(mod_formula, data = data,
                 n.trees = 1000, interaction.depth = 4)

# Summary and relative influence plots
summary(brt_model, plotit = TRUE)

# Partial dependence plots
plot(brt_model, i = 'LST', lwd = 2, main = "")
plot(brt_model, i = 'tpi_500', lwd = 2, main = "")
plot(brt_model, i = 'slope', lwd = 2, main = "")

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


## GLM ----

glm1 <- glm(burned_fraction ~ .,
            # family = quasibinomial('logit'),
            data = data)
summary(glm1)

library(lmtest)
library(sandwich)

# robust std. errors
coeftest(glm1, vcov.=vcovHC(glm1, type="HC0"))

plot(glm1)

pred <- predict(predictors,glm1)
plot(pred,
     range = c(0,1),
     col = scico(30,palette = 'bilbao',direction = -1),
     main = 'Predicted burn fraction - GLM')

## Random forest ----
rf_model <- randomForest::randomForest(burned_fraction  ~ ., data = data)
rf_model

pred <- predict(predictors,rf_model)
plot(pred,
     range = c(0,1),
     col = scico(30,palette = 'bilbao',direction = -1),
     main = 'Predicted burn fraction - RF')
randomForest::varImpPlot(rf_model)

## GAM ----
library(mgcv)
gam1 <- gam(burned_fraction ~  s(LST) + s(NDVI) + 
              s(elevation) +s(tpi_500) + s(aspect) + s(slope),
            data = data)
summary(gam1)

pred <- predict(predictors,gam1)
plot(pred,
     range = c(0,1),
     col = scico(30,palette = 'bilbao',direction = -1),
     main = 'Predicted burn fraction - RF')


# Plots ----
# ..........

## Plot all predictors ----

# Define plotting parameters
parameters <- tibble(palettes = list(elevation = 'grayC',
                                     slope = "buda",
                                     aspect = "romaO",
                                     northness = 'vik',
                                     eastness = 'vik',
                                     tpi_500 = "bam",
                                     LST = "lajolla",
                                     NDVI = "cork",
                                     burned_fraction = "bilbao"),
                     labels = list(elevation = 'Elevation (m)',
                                   slope = 'Slope (°)',
                                   aspect = 'Aspect (°)',
                                   northness = 'Northness (unitless)',
                                   eastness = 'Eastness (unitless)',
                                   tpi_500 = expression(TPI[500~m]),
                                   LST = 'Land surface temperature (° C)',
                                   NDVI = 'NDVI (unitless)',
                                   burned_fraction = NULL)
                     )

create_ggplot <- function(layer, parameters) {
  lyr_name <- names(layer)
  
  direction <- ifelse(lyr_name == 'burned_fraction',-1,1)
  axis_labels <- ifelse(lyr_name == 'burned_fraction','10 km',"")
  
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

ggplots <- lapply(mod_vars, function(layer) {
  create_ggplot(predictors_unscaled[[layer]], parameters)
})

pred_plots <- cowplot::plot_grid(plotlist = ggplots, ncol = 3, 
                                 rel_heights = rep(1, length(ggplots)),
                                 align = 'vh'); print(pred_plots)

ggsave2(pred_plots,filename = 'figures/model/model_predictors_aoi2_v3.png',
        height = 16, width = 16)

## Scatterplot of observed vs. predicted burn fraction ----
data_df <- data.frame(obs = values(fraction_burned,bp),
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

## Compare histograms of LST products ----
clrs <- c("#5A4A6F",  "#EBB261")

pl <- ggplot() +
  geom_density(data = lst_og, aes(fill = clrs[1], x = LST_og),
               color = clrs[1],alpha = .4) +
  geom_density(data = lst, aes(fill = clrs[2],x = LST),
               color = clrs[2],alpha = .4) +
  scale_fill_manual(values = clrs,name = NULL, 
                    labels = c("Level 2 product (30 m x 30 m)", "Resampled to 100 m x 100 m")) +
  labs(title = "Pre-fire Landsat LST",
       x = "Surface temperature (°C)",
       y = "Density") +
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.box = "horizontal", legend.justification = "center") 

ggsave(pl,filename = 'figures/aoi2_LST_densities.png',bg = 'white')
  

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

## Plot water area products ----
blue_colors <- brewer.pal(n = 3, name = "Blues")[-2]

# Load water mask in native resolution
wa_og <- rast(wa_path) 

# Set categories of water mask
cls <- data.frame(id=1:2, cover=c('no_water','water'))
levels(wa) <- cls

p1 <- ggplot() +
  geom_spatraster(data = wa_og) +
  scale_fill_manual(values = blue_colors) +
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.box = "horizontal", legend.justification = "center") +
  labs(fill = NULL,
       title = "Water area (3 m x 3 m)")

p2 <- ggplot() +
  geom_spatraster(data = wa) +
  scale_fill_manual(values = blue_colors) +  
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.box = "horizontal", legend.justification = "center")  +
  labs(fill = NULL,
       title = "Water area (100 m x 100 m)")

combined_plot <- cowplot::plot_grid(p1, p2, ncol = 2,align = 'h');print(combined_plot)

ggsave2(combined_plot,filename = 'figures/water_mask_aggregation_aoi2.png',
        bg = 'white',
        width = 12, height = 8)
