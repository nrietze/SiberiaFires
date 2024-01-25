library(terra)
library(tidyverse)
library(tidyterra)
library(sf)
library(spatialEco)
library(caret)
library(randomForest)
library(cowplot)
library(scico)
library(glcm)
library(raster)
set.seed(10)

# Load ... ----

# ...areas of interest
aois <- read_sf('./data/geodata/feature_layers/aoi_wv/aois.shp') %>%
  mutate(., id = 1:nrow(.))

# ... list of raster files
pre_year <- 2019
post_year <- 2020

raster_files_pre <- list.files(sprintf('C:/data/8_planet/%s/cropped/',pre_year),
                               pattern = '*composite.tif$',
                               full.names = T
)
raster_files_post <- list.files(sprintf('C:/data/8_planet/%s/cropped/',post_year),
                                pattern = '*composite.tif$',
                                full.names = T
)

regex_name <- 'Kosukhino'
fn_raster_pre <- raster_files_pre[grepl(regex_name,raster_files_pre)]

# ... raster
scaling_factor <- 1e4

r <- rast(fn_raster_pre) / scaling_factor

v <- vect('data/geodata/feature_layers/landform/histogram_aois.shp') %>% 
  mutate(., ID = 1:nrow(.))
r_c <- crop(r, v[1])

# Plot results ----
ggplot() +
  geom_spatraster(data=r_c) +
  facet_wrap(~lyr)+
  scale_fill_scico(palette = 'grayC') +
  theme_cowplot()

f <- hist(r_c, breaks=30,plot = F)
dat <- do.call(rbind, lapply(f, function(hist_data) {
  data.frame(counts = hist_data$counts, breaks = hist_data$mids, layer = hist_data$xname)
}))

ggplot(dat, aes(x = breaks, y = counts)) + 
  facet_wrap(~layer)+
  geom_bar(stat = "identity",fill='blue',alpha = 0.8)+
  xlab("Reflectance") + ylab("Frequency")+
  scale_x_continuous(breaks = seq(0,1,.1), 
                     labels = seq(0,1,.1)) +
  theme_cowplot()

ndvi <-  (r['NIR'] - r['Red'])/(r['NIR'] + r['Red'])
ndvi_c <- crop(ndvi, v[1])

ggplot() +
  geom_spatraster(data=ndvi) +
  scale_fill_scico(palette = 'bam') +
  theme_cowplot()

ndvi_coma <- glcm(raster(ndvi_c), window = c(33, 33), shift = c(1, 0), na_opt = "any")

ndvi_coma_rast <- rast(ndvi_coma)



# ***
raster_data <- terra::extract(r['NIR'], v,df = T) %>% 
  pivot_longer(cols = -ID, names_to = "Bandname", values_to = "Value") %>% 
  left_join(.,as.data.frame(v), by = 'ID') %>% 
  mutate(.,ID=as.factor(ID))

ggplot(raster_data, aes(x = Value, color = landform,fill = landform)) +
  # geom_histogram(binwidth = 0.01, position = "identity", alpha = 0.7) +
  geom_density(size = 1,alpha = .1) +
  facet_wrap(~ Bandname) +
  theme_cowplot()
