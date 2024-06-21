library(terra)
library(tidyverse)
library(tidyterra)
library(tidybayes)
library(brms)
library(broom.mixed)

# 1. Configuration ----
window_side_length <- 30 # in metres

model_all_aois <- TRUE  # model all AOIs or just one?
use_planet_wa <- FALSE    # use Planet water mask (=TRUE) or Landsat's QA mask (=FALSE)?
predict_bf <- FALSE     # run prediction analysis over LarcgeScarControl?

# Load areas of interest
aois <- vect('./data/geodata/feature_layers/aoi_wv/aois_analysis.geojson') %>%
  project('EPSG:32655') %>% 
  mutate(., id = 1:nrow(.))

# Provide a part of the AOI name for processing
if (model_all_aois){
  aoi_names <- aois$site
} else {
  aoi_names <- c('LargeScarCenter')
}

# 2. Prepare data ----
prepare_data <- function(aoi_name, window_side_length){
  
  ## Get current AOI outline ----
  aoi <- aois[aois$site == aoi_name]
  
  cat(sprintf('Preparing model data for %s ... \n',aoi$site) )
    
  # Get processing mask in current AOI if it's AOI 4 or 5
  if (aoi$site %in% c("Berelech","LargeScarCenter")){
    poly_mask <- vect('data/geodata/feature_layers/planet_masks.shp') %>% 
      terra::intersect(aoi)
    
    has_mask <- TRUE

  } else{has_mask <- FALSE}
  
  # Load burn perimeter
  bp <- vect(sprintf('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_%s.shp',aoi$site) ) %>%
    crop(aoi)

  ## Load Landsat data ----
  cat('Loading Landsat data... \n')
  ls_path <- 'data/geodata/raster/landsat/'
  
  clear_value <- 21824
  
  ### LST ----
  lst_scale_factor <- 0.00341802 	
  lst_offset <- 149 - 273.15
  
  if (aoi$site %in% c("LargeScarControl","LargeScarCenter","Libenchik")){
    qa <- rast(paste0(ls_path,'LC08_L2SP_116010_20200615_20200824_02_T1_QA_PIXEL.TIF')) %>% 
      crop(aoi)
    
    lst <- rast(paste0(ls_path,'LC08_L2SP_116010_20200615_20200824_02_T1_ST_B10.TIF')) %>% 
      crop(aoi) * lst_scale_factor + lst_offset
    
    # Mask only clear observations
    lst <- mask(lst,qa,maskvalues = clear_value,inverse = TRUE)
    
  } else {
    qa <- rast(paste0(ls_path,'LC08_L2SP_115010_20200608_20200824_02_T1_QA_PIXEL.TIF')) %>% 
      crop(aoi)
    
    lst <- rast(paste0(ls_path,'LC08_L2SP_115010_20200608_20200824_02_T1_ST_B10.TIF')) %>% 
      crop(aoi) * lst_scale_factor + lst_offset
    
    # Mask only clear observations
    lst <- mask(lst,qa,maskvalues = clear_value,inverse = TRUE)
  }
  names(lst) <- 'LST'
  
  # Create raster grid template for PLanetScope & ArcticDEM data
  if (window_side_length == 30 ){ # Use landsat raster as template if grid cells are 30 x 30 m
    raster_grid_template <- rast(ext(lst), resolution=res(lst)) 
    crs(raster_grid_template) <- crs(lst)
  } else {
    ncells_resample <- window_side_length / 2 # e.g., if window is 100 m wide, a 2 m cell size gives 50 cells
    raster_grid_template <- aggregate(lst, fact = ncells_resample, fun = 'mean')
  }
  
  # Resample LST if resolution is not 30 m
  if (window_side_length != 30){
    lst <- resample(lst,raster_grid_template,'bilinear')
  }
  
  ### NDVI ----
  msp_scale_factor <- 2.75e-05 	
  msp_offset <- -0.2
  
  # Level 2
  if (aoi$site %in% c("LargeScarControl","LargeScarCenter")){
    qa <- rast(paste0(ls_path,'LC08_L2SP_116010_20200615_20200824_02_T1_QA_PIXEL.TIF')) %>% 
      crop(aoi)
    red <- rast(paste0(ls_path,'LC08_L2SP_116010_20200615_20200824_02_T1_SR_B4.TIF')) %>% 
      crop(aoi) 
    nir <- rast(paste0(ls_path,'LC08_L2SP_116010_20200615_20200824_02_T1_SR_B5.TIF')) %>% 
      crop(aoi)
  } else {
    qa <- rast(paste0(ls_path,'LC08_L2SP_115010_20200608_20200824_02_T1_QA_PIXEL.TIF')) %>% 
      crop(aoi)
    red <- rast(paste0(ls_path,'LC08_L2SP_115010_20200608_20200824_02_T1_SR_B4.TIF')) %>% 
      crop(aoi) 
    nir <- rast(paste0(ls_path,'LC08_L2SP_115010_20200608_20200824_02_T1_SR_B5.TIF')) %>% 
      crop(aoi)
  }

  ndvi <- ((nir - red) / (nir + red)) 
  names(ndvi) <- 'NDVI'
  
  # Mask only clear observations
  ndvi <- mask(ndvi,qa,maskvalues = clear_value,inverse = TRUE)
  
  # Resample NDVI if resolution is not 30 m
  if (window_side_length != res(ndvi)[1]){
    ndvi <-  resample(ndvi,raster_grid_template,'bilinear')
  }
  
  ## Load DEM data ----
  cat('Loading ArcticDEM data... \n')
  dem_path <- 'data/geodata/raster/arcticDEM/'
  
  dem_files <- list.files(dem_path,pattern = 'aoi.*_dem_.*utm\\.tif$',full.names = T)
  
  raster_index <- grep(paste0("aoi_", aoi$site, "_"), dem_files)
  
  dem_og <- rast(dem_files[raster_index]) %>% crop(aoi)
  
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
  raster_index <- grep(paste0("aoi_", aoi$site, "_"), tpi_files)
  
  tpi_500 <- rast(tpi_files[raster_index]) %>% 
    crop(aoi) %>% 
    resample(.,raster_grid_template,method = 'cubicspline')
  names(tpi_500) <- 'tpi_500'
  
  ## PlanetScope data ----
  cat('Loading PlanetScope data... \n')
  
  ### NDVI heterogeneity ----
  planet_path <- list.files('C:/data/8_planet/2019/cropped/',
                            pattern = paste0(aoi$site,'.*composite\\.tif$'),
                            full.names = TRUE)
  
  planet_msp <- rast(planet_path) %>% 
    crop(aoi) 
  
  get_ndvi_heterogeneity <- function(red, nir){
    ndvi <- ((nir - red) / (nir + red)) 
    
    # Compute patchiness metric before resampling
    ndvi_sd <- terra::aggregate(ndvi, fact = 10,fun = 'sd') %>% 
      resample(.,raster_grid_template,'bilinear')
    names(ndvi_sd) <- 'NDVI_sd'
    
    return(ndvi_sd)
  }
  
  ndvi_sd <- get_ndvi_heterogeneity(planet_msp$Red,planet_msp$NIR)
  
  ### Burned area  ----
  ba_path <- list.files('data/geodata/raster/burned_area/planet/',
                        pattern = paste0(aoi$site,'.*burned_area_top5TD\\.tif$'),
                        full.names = TRUE)
  # calculate burned fraction
  ba <- rast(ba_path) %>% 
    crop(aoi)
  
  fraction_burned <- ifel(ba == 'burned', 1, 0) %>% 
    resample(.,raster_grid_template,'sum') / (window_side_length / res(ba)[1])**2
  names(fraction_burned) <- 'burned_fraction'
  
  fraction_unburned <- ifel(ba == 'unburned', 1, 0) %>% 
    resample(.,raster_grid_template,'sum') #/ (window_side_length / res(ba)[1])**2
  names(fraction_unburned) <- 'unburned_fraction'
  
  ### water mask ----
  if(use_planet_wa){
    wa_path <- list.files('data/geodata/raster/water_area/planet/',
                        pattern = paste0(aoi$site,'.*water_area_top5TD\\.tif$'),
                        full.names = TRUE)
    wa <- rast(wa_path) %>% 
      crop(aoi) %>% 
      resample(raster_grid_template,'mode')
    
    writeRaster(wa,
                filename = paste0("data/geodata/raster/water_area/",aoi$site,
                                  "_resampledPS_mask.tif"),
                overwrite = T
    )
  } else{
    wa <- ifel(qa == 21824,1,2)       # if no clear observation, pixel is masked out
    plot(wa)
    
    writeRaster(wa,
                filename = paste0("data/geodata/raster/water_area/",aoi$site,
                                  "_Landsat_mask.tif"),
                overwrite = T
    )
  }
  
  ## gather all predictors and mask water areas ----
  predictors <- c(fraction_burned,
                  dem,slope, aspect,northness,eastness, tpi_500,
                  lst,ndvi,ndvi_sd) %>% 
    mask(.,wa,maskvalue = 2, updatevalue = NA)
  
  cat("Writing raster of all predictors ...\n")
  writeRaster(predictors,
              filename = paste0("data/geodata/raster/predictors/",aoi$site,
                                "_predictors_",window_side_length,"m.tif"),
              overwrite = T
  )
  
  # Get data from random sample of raster values
  nsample <- 3e3
  
  # Sample n pixels and mask with burn perimeter
  random_points <- spatSample(raster_grid_template,nsample, 
                              as.points = TRUE,values = FALSE) %>% 
    mask(bp) 
  
  if (has_mask){
    random_points <- random_points %>% 
      mask(poly_mask,inverse = T)
  }
  
  # Extract data at random points 
  data <- terra::extract(predictors, random_points,ID = FALSE, xy = TRUE) %>% 
    drop_na() 
  data$site <- aoi$site
  
  cat('Exporting data to csv... \n')
  
  # Export dataframe to csv
  write.csv(data, file = paste0("tables/predictors/",aoi$site,
                                "_predictors_",window_side_length,"m.csv"))
  
  cat('done. \n')
  
  return(NULL)
}

# call function
load_data <- function(aoi_name, path){
  csv_path <- paste0(path,aoi_name,"_predictors_",window_side_length,"m.csv")
  
  # Check if file exists
  if (file.exists(csv_path)){
    cat(sprintf('Data found for %s, loading data ... \n',aoi_name))
    
    data <- read.csv(csv_path)
    
    return(data)
  } else{
    cat('No data frame could be found for this site. Preparing data first... \n')
    prepare_data(aoi_name,window_side_length)
    
    data <- read.csv(csv_path)
    return(data)
  }
 
}

# 3. Apply model ----
use_all_points <- TRUE

if (use_all_points){
  # Merge all dataframes together
  model_data <- do.call(rbind,
                        lapply(aoi_names, FUN = load_data, path = "tables/predictors/")) %>% 
    dplyr::select(-X)
  suffix <- ''
} else {
  # load model data with 400 m distance contstraint
  model_data <- read.csv("tables/predictors/model_data_400m.csv")
  suffix <- 'thinned_data_'
}

## a) Scale data ----
model_data_scaled <- model_data %>% 
  dplyr::select(-c(burned_fraction,x,y)) %>% 
  mutate(across(where(is.numeric), \(x) scale(x, center = TRUE, scale = TRUE)) ) %>% 
  mutate(burned_fraction = model_data$burned_fraction)

print(summary(model_data_scaled))

model_data_scaled$site <- as.factor(model_data_scaled$site)

# Define predictor variables for the model
mod_vars <- c('elevation','slope','northness','eastness','tpi_500', 
              'LST','NDVI_sd','NDVI')

mod_formula <- formula(paste('burned_fraction ~ (1| site) +',paste(mod_vars, collapse = '+')))
phi_formula <- formula(paste('phi ~ (1| site) +',paste(mod_vars, collapse = '+')))
zoi_formula <- formula(paste('zoi ~ (1| site) +',paste(mod_vars, collapse = '+')))
coi_formula <- formula(paste('coi ~ (1| site) +',paste(mod_vars, collapse = '+')))

# Check for multicollinearity using a linear model
lm_mod <- lme4::lmer(mod_formula, data = model_data_scaled)
performance::multicollinearity(lm_mod)

## b) Run zero-one inflated beta regression ----
n_iter <- 10000
n_thin <- ifelse(n_iter > 5000, 10, 1)

zoi_beta_model <- brm(
  bf(mod_formula,
     phi_formula,
     zoi_formula,
     coi_formula),
  data = model_data_scaled,
  family = zero_one_inflated_beta(),
  chains = 4, 
  iter = n_iter, 
  # warmup = 1000,        # use default warmup (iter/2)
  thin = n_thin,
  cores = 4, seed = 1234,
  file = paste0("zoi_beta_model_",suffix,today())
)

# fit Null model
null_model <- brm(
  bf(burned_fraction ~ 1), 
  family = zero_one_inflated_beta(), 
  data = model_data_scaled,chains = 4, 
  iter = n_iter, 
  thin = n_thin,
  cores = 4, seed = 1234,
  file = paste0("null_model_",today())
)

## c) Diagnose model ----
summary(zoi_beta_model)
summary(null_model)

# 4. Predict burned fraction in Large Scar Control AOI ----
if (predict_bf){
  aoi_names_subset <- c('LargeScarCenter','Kosukhino','Berelech','DrainedThawlake')
  
  # load model data
  model_data_subset <- do.call(rbind,
                               lapply(aoi_names_subset, FUN = load_data, path = "tables/predictors/"))
  
  model_data_subset_scaled <- model_data_subset %>% 
    select(-X) %>% 
    dplyr::select(-burned_fraction) %>% 
    mutate(across(where(is.numeric), scale, center = TRUE, scale = TRUE)) %>% 
    mutate(burned_fraction = model_data_subset$burned_fraction) %>% 
    mutate(site = as.factor(site)) 
  
  
  print(summary(model_data_subset_scaled))
  
  # use RandomForest for prediciton?
  use_RF <- TRUE
  
  if(use_RF){
    library(randomForest)
    pred_model <- randomForest(burned_fraction ~ elevation + slope + LST + NDVI_sd + northness,
                               data = model_data_subset_scaled)
    print(pred_model)
    
  } else {
    pred_model <- brm(
      bf(burned_fraction ~ elevation + slope + LST + NDVI_sd),
      data = model_data_subset,
      family = zero_one_inflated_beta(),
      chains = 4, 
      iter = n_iter, 
      # use default warmup (iter/2)
      # warmup = 1000, 
      thin = n_thin,
      cores = 4, seed = 1234,
      file = paste0("mod_no_LSControl_",today())
    )
    
    # assess model
    summary(pred_model)
    
    plot(pred_model)
  }
  
  # load predictors
  LSControl_predictors <- rast(paste0("data/geodata/raster/predictors/LargeScarControl_predictors_",
                                      window_side_length,"m.tif") )
  
  predictors_scaled <- LSControl_predictors %>% 
    tidyterra::select(-burned_fraction) %>% 
    mutate(site = 'LargeScarControl') %>% 
    scale()
  
  # predict burned fraction in entire map
  preds <- terra::predict(pred_model,predictors_scaled,cores = 4)
  
  preds_rast_df <- as.data.frame(predictors_scaled$elevation,xy = T) 
  
  if(use_RF){
    preds_rast_df$burned_fraction_pred <- preds
  } else {
    # preds_rast_df$burned_fraction_pred <- plogis(preds[,1])
    preds_rast_df$burned_fraction_pred <- preds[,1]
  }
  
  preds_rast <- rast(preds_rast_df,type = 'xyz',
                     crs = "EPSG:32655",extent = ext(predictors_scaled)) %>% 
    subset('burned_fraction_pred') %>% 
    mask(is.na(predictors_scaled$NDVI_sd),maskvalues = TRUE)
  
  # Assign observed burned fraction
  preds_rast$burned_fraction_observed <- LSControl_predictors$burned_fraction
  
  # export predicted burned fraction
  # writeRaster(preds_rast,
  #             filename = "data/geodata/raster/burned_area/planet/predicted_BF_LSControl.tif",
  #             overwrite = T
  # )
  
  # Load burn perimeter
  bp <- vect('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_LargeScarControl.shp') %>% 
    crop(aois[aois$site == 'LargeScarControl' ])
  
  # Plot results of prediction
  ggplot() + 
    geom_spatraster(data = preds_rast) +
    facet_wrap(~lyr) +
    geom_spatvector(data = bp,color = 'black', fill = NA, size = 2) +
    scale_fill_scico(palette = 'bilbao',direction = -1,na.value="white") +
    theme_cowplot()
  
  ggsave('figures/model/LargeScarControl_pred_obs_bf_RF.png',
         width = 16, height = 8,bg = 'white')
  
  # scatterplot obs. vs. pred
  data_df <- data.frame(value1 = values(mask(LSControl_predictors$burned_fraction,bp)), 
                        value2 = values(mask(preds_rast,bp))) %>% 
    drop_na()
  
  # Scatterplot
  ggplot() + 
    geom_point(data = data_df, aes (x = data_df[,1], y = data_df[,2])) +
    labs(x = "observed burn fraction",y = "predicted burn fraction") +
    lims(x = c(0,1),y = c(0,1)) +
    theme_cowplot()
  
  ggsave('figures/model/LargeScarControl_pred_obs_bf_scatter_RF.png',
         width = 8, height = 8,bg = 'white')
}

# Logitsic regression and random forest model on 0s and 1s ----
library(randomForest)
model_data_scaled_zo <- model_data_scaled %>% 
  filter(burned_fraction == 1 | burned_fraction == 0) %>% 
  mutate(burned_fraction = as.factor(burned_fraction))

rf_zo <- randomForest(burned_fraction ~ ., data = model_data_scaled_zo, importance=TRUE)
rf_zo

n_iter <- 10000
n_thin <- ifelse(n_iter > 5000, 10, 1)

lr_mod_zo <- brm(
  bf(mod_formula),
  data = model_data_scaled_zo,
  family = bernoulli(),
  chains = 4, 
  iter = n_iter, 
  thin = n_thin,
  cores = 4, seed = 1234,
  file = paste0("log_reg_model_zo_",suffix,today())
)
summary(lr_mod_zo)
performance::r2_bayes(lr_mod_zo)

varImpPlot(rf_zo,main = 'Variable importances of RF on 0/1 burned fractions')

df_rf_zo <- data.frame(predicted = predict(rf_zo,model_data_scaled_zo),
                          observed = model_data_scaled_zo$burned_fraction)
rmse <- sqrt(mean((df_rf_zo$observed - df_rf_zo$predicted) ^ 2))
print(rmse)
r2 <- cor(df_rf_zo$predicted,df_rf_zo$observed) ^2
print(r2)

ggplot(df_rf_zo) +
  geom_abline(intercept = 0,slope = 1,color = "gray",linewidth = 1) +
  geom_point(aes(x = observed, y = predicted)) +
  labs(x = 'Observed fractions either 0 and 1',y = 'Predicted fractions either 0 and 1') +
  coord_equal() +
  theme_cowplot()

# Beta and random forest regression on values in between 0 & 1 ----
model_data_scaled_no_zo <- model_data_scaled %>% 
  filter(burned_fraction < 1 & burned_fraction > 0) 

beta_mod <- brm(
  bf(mod_formula,
     phi ~ (1 | site) + elevation + slope + northness + 
       eastness + tpi_500 + LST + NDVI_sd + NDVI),
  data = model_data_scaled_no_zo,
  family = Beta(),
  chains = 4, 
  iter = n_iter, 
  thin = n_thin,
  cores = 4, seed = 1234,
  file = paste0("beta_reg_model_",suffix,today())
)

summary(beta_mod)
performance::r2_bayes(beta_mod)


rf_no_zo <- randomForest(burned_fraction ~ ., data = model_data_scaled_no_zo,
                         importance=TRUE)
rf_no_zo

varImpPlot(rf_no_zo,main = 'Variable importances of RF on continuous burned fractions')

df_rf_no_zo <- data.frame(predicted = predict(rf_no_zo,model_data_scaled_no_zo),
                          observed = model_data_scaled_no_zo$burned_fraction)

rmse <- sqrt(mean((df_rf_no_zo$observed - df_rf_no_zo$predicted) ^ 2))
print(rmse)
r2 <- cor(df_rf_no_zo$predicted,df_rf_no_zo$observed) ^2
print(r2)

ggplot(df_rf_no_zo) +
  geom_abline(intercept = 0,slope = 1,color = "gray",linewidth = 1) +
  geom_point(aes(x = observed, y = predicted)) +
  labs(x = 'Observed fractions between 0 and 1',y = 'Predicted fractions between 0 and 1') +
  coord_equal() +
  theme_cowplot()


# RF classifier of 0 or 1 or in-between fraction
model_data_scaled_cl <- model_data_scaled %>% 
  mutate(class_bf = ifelse(burned_fraction < 1 & burned_fraction > 0,2,burned_fraction) ) %>% 
  mutate(class_bf = as.factor(class_bf))

model_data_scaled_cl$id <- 1:nrow(model_data_scaled_cl)
training <- model_data_scaled_cl %>%
  slice_sample(prop = 0.8)
validation <- filter(model_data_scaled_cl, !(id %in% training$id)) 

rf_cl <- randomForest(class_bf ~ ., 
                      data = training %>% 
                        select(-c(burned_fraction,id)),
                      importance=TRUE)
rf_cl
varImpPlot(rf_cl,main = 'Variable importances of RF 3 classes of burned fractions')

caret::confusionMatrix(predict(rf_cl,validation),validation$class_bf)

model_data_scaled_cl$class_pred <- predict(rf_cl,model_data_scaled_cl)
training2 <- filter(model_data_scaled_cl, (id %in% training$id)) 
validation2 <- filter(model_data_scaled_cl, !(id %in% training$id)) 

rf_cl_inbetween <- randomForest(burned_fraction ~ ., 
                      data = model_data_scaled_cl %>% 
                        filter(class_pred == 2) %>% 
                        select(-c(class_bf,class_pred,id)),
                      importance=TRUE)
rf_cl_inbetween
varImpPlot(rf_cl_inbetween,main = 'Variable importances of RF for in-between burned fractions')

df_rf_cl_inbetween <- data.frame(predicted = predict(rf_cl_inbetween,
                                                     filter(model_data_scaled_cl,class_pred == 2)),
                                 observed = model_data_scaled_cl %>% 
                                   filter(class_pred == 2) %>% 
                                   pull(burned_fraction) )

rmse <- sqrt(mean((df_rf_cl_inbetween$observed - df_rf_cl_inbetween$predicted) ^ 2))
print(rmse)
r2 <- cor(df_rf_cl_inbetween$predicted,df_rf_cl_inbetween$observed) ^2
print(r2)

ggplot(df_rf_cl_inbetween) +
  geom_abline(intercept = 0,slope = 1,color = "gray",linewidth = 1) +
  geom_point(aes(x = observed, y = predicted)) +
  labs(x = 'Observed fractions between 0 and 1',y = 'Predicted fractions between 0 and 1') +
  coord_equal() +
  theme_cowplot()
