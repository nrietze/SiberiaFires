library(terra)
library(tidyterra)
library(cowplot)
library(tidyverse)
library(landscapemetrics)
library(motif)
library(scico)

## load burned area ----
ba <- rast('data/geodata/raster/burned_area/planet/Berelech_2020_burned_area.tif')
ba <- rast('data/geodata/raster/burned_area/planet/Kosukhino_2020_burned_area.tif')
ba <- rast('data/geodata/raster/burned_area/planet/LargeScarCenter_2020_burned_area.tif')

## Fire perimeteres from Talucci et al. (2022) ----
fire_perimeter <- vect('C:/data/3_fire_data/burned_area/siberia_talucci/data/SiberiaFires/SiberiaFires2001-2020_wgs.shp') %>% 
  filter(FireId == 'RU_2020_NSCT_20200615') %>% 
  project('EPSG:32655')
fire_perimeter <- fire_perimeter[1]

## load water area ----
wa_path <- 'data/geodata/raster/water_area/planet/Kosukhino_2020_water_area.tif'
wa <- rast(wa_path) 

ba <- mask(ba, wa, maskvalues = 2,updatevalue = NA)

check_landscape(ba)

# Plot burned area
plot(ba)

# Get individual patches and plot it
patched_raster <-  get_patches(ba,return_raster = T)
plot(patched_raster$layer_1$class_1,main = 'unburned patches')
plot(patched_raster$layer_1$class_2,main = 'burned patches')

# Compute patch areas
patch_area <- lsm_p_area(ba)
patch_area$area_m2 <- patch_area$value * 1e4

# Compute patch density in same areas


# Plot inverese cumulative distribution 
(patch_plot <- patch_area %>% 
  ggplot(aes(x = area_m2,color = factor(class))) +
  stat_ecdf(geom = "line", pad = FALSE, 
            aes(y = after_stat(1-y)), size = 1) +
  scale_color_manual(values = c("green", "black"), labels = c("unburned", "burned")) +
  scale_x_continuous(expand = c(0, 0),limits = c(-10,1e3),trans='log10') +
  scale_y_continuous(expand = c(0, 0.05), limits = c(0,1)) +
  labs(x = "Patch area (m²)", y = "Inverse Cumulative Proportion",color = NULL) +
  theme_cowplot())

ggsave(patch_plot,filename = 'figures/aoi2_patch_area.png')

# Patch distance from unburned to burned
distance_to_seed <- distance(ba,target = 1,exclude = NA) %>% 
  mask(bp)

dist_rast <- ggplot() + 
    geom_spatraster(data = distance_to_seed) +
    scale_fill_scico(palette = 'navia') + 
    labs(fill = 'distance (m)')+
    theme_cowplot()

dist_plot <- distance_to_seed %>% 
    as.data.frame() %>% 
    ggplot(aes(x = class)) +
    stat_ecdf(geom = "line", pad = FALSE, 
            aes(y = after_stat(1-y)), size = 1) +
    lims(x = c(0,400),y = c(0,1)) +
    labs(x = "Distance to unburned pixel (m)", y = "Inverse Cumulative Proportion",
         color = NULL) +
    theme_cowplot()

(dist_plots <- plot_grid(plotlist = list(dist_rast,dist_plot),
                         rel_widths = c(1,1),
                         ncol = 2,align = 'v',axis = 'b'))
