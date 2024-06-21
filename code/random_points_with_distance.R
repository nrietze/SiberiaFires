# Load required library
library(terra)
library(sp)
library(OSMscale)
library(tictoc)
# Assuming your dataframe is named model_data and has columns named 'x' for longitude and 'y' for latitude

# Convert dataframe to a simple feature (sf) object
model_sf <- vect(model_data, geom=c("x", "y"),crs = 'epsg:32655')

set.seed(42)  # For reproducibility

# using a distance matrix
is_far_enough <- function(point, selected_points,mindist) {
  distances <- terra::distance(selected_points, point)
  all(distances > mindist)
}

get_random_points <- function(aoi_name,data,sample_size,mindist){
  tic()
  
  max_iterations <- 1000
  iterations <- 0
  
  print(aoi_name)
  
  initial_point <- data %>% 
    filter(site == aoi_name) %>% 
    as.data.frame() %>% 
    sample_n(1) %>% 
    vect(., geom=c("x", "y"),crs = 'epsg:32655')
  
  selected_points <- initial_point
  
  while (nrow(selected_points) < sample_size && iterations < max_iterations) {
    random_point <- data %>% 
      filter(site == aoi_name) %>% 
      as.data.frame() %>% 
      sample_n(1) %>% 
      vect(., geom=c("x", "y"),crs = 'epsg:32655')
    
    if (is_far_enough(random_point, selected_points,mindist)) {
      selected_points <- rbind(selected_points, random_point)
    }
    
    iterations <- iterations + 1
  }
  
  if (iterations >= max_iterations) {
    print("Maximum iterations reached. Could not find enough points.")
  }
  
  selected_points <- as.data.frame(selected_points) %>% 
    mutate(x = data.frame(geom(selected_points))$x,
           y = data.frame(geom(selected_points))$y)
  toc()
  return(selected_points)
}


mindist <- 400
sample_size <- 50

model_data_400m <- do.call(rbind,
                          lapply(aoi_names,get_random_points,
                                 data = model_data,
                                 sample_size = sample_size,
                                 mindist=mindist))

write.csv(model_data_400m,file = "tables/predictors/model_data_400m.csv")

residuals_400m <- do.call(rbind,
                          lapply(aoi_names,get_random_points,
                                 data = df_res,
                                 sample_size = sample_size,
                                 mindist=mindist))

# Plot points in space
model_data_400m %>% 
  ggplot() + 
  geom_point(aes(x = x,y = y,color = site)) + 
  facet_wrap(~site,scales = 'free') +
  theme_cowplot()

bp <- vect('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_Kosukhino.shp')

model_data_400m %>% 
  filter(site =='Kosukhino') %>% 
  ggplot() + 
  geom_point(aes(x = x,y = y,color = site)) + 
  geom_spatvector(data = bp,fill = NA) +
  theme_cowplot()

# Plot semivariograms
residuals_400m %>% 
  group_by(site) %>% 
  group_map(~geostats::semivariogram(x = .x$x,y =.x$y, 
                                     z = .x$residuals,
                                     main = .y,
                                     bw = 10,nb = 100) ) 

# convert to spatial dataframe
for (name in aoi_names){
  print(name)
  spatial_residuals <-  residuals_400m %>% 
    filter(site == name) %>% 
    sf::st_as_sf(coords = c('x','y'),crs = st_crs(32655))

  lstw  <- nb2listw(knn2nb(knearneigh(spatial_residuals, k = 8)), # create neighbors
                    style="W")                                     # row standardized weights
  print(moran.test(spatial_residuals$residuals, lstw))
}
