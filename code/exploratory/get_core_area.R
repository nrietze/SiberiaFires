library(terra)
library(tidyterra)
library(cowplot)
library(scales)
library(tidyverse)
library(landscapemetrics)
library(motif)
library(scico)

# 1. Configuration & loading data ----
plot_all_aois <- TRUE

# Provide a part of the AOI name for processing
if (plot_all_aois){
  aoi_names <- c('LargeScarCenter','LargeScarControl','Kosukhino','Berelech','DrainedThawlake')
} else {
  aoi_names <- c('LargeScarCenter')
}

# Load areas of interest
aois <- vect('./data/geodata/feature_layers/aoi_wv/aois_analysis.geojson') %>%
  project('EPSG:32655') %>% 
  mutate(., id = 1:nrow(.))

get_patch_area <- function(aois,aoi_name){
  
  # Get current AOI
  aoi <- aois[aois$site == aoi_name]
  
  cat(sprintf('Gathering patch areas for %s ... \n',aoi$site) )
  
  # Load burned area
  ba_path <- sprintf('data/geodata/raster/burned_area/planet/%s_burned_area_top5TD.tif',aoi_name)
  ba <- rast(ba_path)
  
  # Load water mask
  wa_path <- sprintf('data/geodata/raster/water_area/planet/%s_water_area_top5TD.tif',aoi_name)
  wa <- rast(wa_path) 
  
  # Load burn perimeter
  bp <- vect(sprintf('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_aoi%s.shp',aoi$site) ) %>%
    crop(aoi)
  
  # Mask out water pixels
  ba <- ba %>% 
    mask(wa, maskvalues = 2,updatevalue = NA) %>% 
    mask(bp,updatevalue = NA)
  
  # Load invalid data masks
  if (aoi$id %in% c(4,5)){
    poly_mask <- vect('data/geodata/feature_layers/planet_masks.shp') %>% 
      terra::intersect(aoi)
    
    ba <- mask(ba,poly_mask,inverse = TRUE)
    } 
  
  print(check_landscape(ba))
  
  # Compute patch areas
  patch_area <- lsm_p_area(ba)
  
  # Convert to square metres
  patch_area$area_m2 <- patch_area$value * 1e4
  
  patch_area$site <- aoi_name
  
  return(patch_area)
}

# Merge all dataframes together
all_patch_areas <- do.call(rbind,
                           lapply(aoi_names, FUN = get_patch_area, aois = aois))
all_patch_areas <- mutate(all_patch_areas, site = factor(site))

# 2. Plot patch area proportions ----
class_to_plot <- 2 # ( 1 = unburned, 2 = burned)

site_colors <- c("#a6cee3", "#1f78b4", "#33a02c", "#b2df8a")

(patch_plot <- all_patch_areas %>% 
    filter(class == class_to_plot) %>% 
    ggplot(aes(x = area_m2,colour = site)) +
    geom_vline(xintercept = 9, linewidth = .5,linetype = "dashed") +
    annotate("text", x = 10, y = 30, label = "1 planet pixel", vjust = 0, hjust = 0, size = 3) +
    geom_vline(xintercept = 81, linewidth = .5,linetype = "dashed") +
    annotate("text", x = 85, y = 40, label = "3 x 3 planet pixels", vjust = 0, hjust = 0, size = 3) +
    geom_vline(xintercept = 900, linewidth = .5,linetype = "dashed") +
    annotate("text", x = 920, y = 60, label = "1 landsat pixel", vjust = 0, hjust = 0, size = 3) +
    stat_ecdf(aes(y = after_stat(y*100)),
              geom = "line",
              pad = FALSE, size = 1) +
    scale_x_continuous(expand = c(1e-2, 0),trans='log10',
                       limits = c(8,1e5),
                       labels = label_number()) +
    scale_color_manual(values = site_colors) +
    labs(x = "Patch area (m²)", y = "Cumulative Percentage (%)",color = 'Site') +
    theme_cowplot() +
    theme(
      legend.position = "bottom",  
      legend.title = element_blank(),
      legend.direction = "horizontal",
    ) 
  )

ggsave(patch_plot,filename ='figures/ECDF_all_patch_areas.png',
       bg = 'white',width = 10, height = 10)


landsat_area_log10 <- log10(30**2)

ecdf_unburned <- all_patch_areas %>%
  filter(class == 1 & area_m2 > 81) %>%
  pull(area_m2) %>%
  log10() %>%
  ecdf()


print(ecdf_unburned(landsat_area_log10)*100)

# Compute patch density in same areas
ggplot(patch_area, aes(x = log10(area_m2), fill = factor(class))) +
  geom_histogram(binwidth = 0.25, alpha = 0.5, position = "identity") +  # Add histogram with transparency
  scale_fill_manual(values = c("red", "green")) +  # Define colors for classes
  labs(x = "Log10 Area (m^2)", y = "Frequency") +  # Label axes
  theme_minimal() +
  xlim(c(0.5, max(log10(patch_area$area_m2))))  # Set x-axis limits to include 0

# Patch distance from unburned to burned
bp <- vect(sprintf('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_aoi%s.shp',5) ) %>%
  crop(aois[5])
distance_to_seed <- distance(ba,target = 1,exclude = NA) %>% 
  mask(bp)

(dist_rast <- ggplot() + 
    geom_spatraster(data = distance_to_seed) +
    scale_fill_scico(palette = 'navia') + 
    labs(fill = 'distance (m)')+
    theme_cowplot())

(dist_plot <- distance_to_seed %>% 
    as.data.frame() %>% 
    ggplot(aes(x = class)) +
    stat_ecdf(geom = "line", pad = FALSE, 
            aes(y = after_stat(1-y)), linewidth = 1) +
    lims(x = c(0,400),y = c(0,1)) +
    labs(x = "Distance to unburned pixel (m)", y = "Inverse Cumulative Proportion",
         color = NULL) +
    theme_cowplot())

(dist_plots <- plot_grid(plotlist = list(dist_rast,dist_plot),
                         rel_widths = c(1,1),
                         ncol = 2,align = 'v',axis = 'b'))


aoi_5_rasters <- rast(paste0("data/geodata/raster/predictors/",
                             aoi_names[1],"_predictors_",
                             window_side_length,"m.tif") )
fraction_burned <- aoi_5_rasters$burned_fraction
extremes <- ifel((fraction_burned<1 & fraction_burned>0), 
                 NA, fraction_burned)
distance_to_unburned_seed <- distance(extremes,target = 1,exclude = NA,unit = 'm')
distance_to_lesserburned_seed <- distance(extremes,target = 1,unit = 'm')


plot(distance_to_seed)
hist(distance_to_unburned_seed[distance_to_unburned_seed>0])
hist(distance_to_lesserburned_seed[distance_to_lesserburned_seed>0],add = T)


writeRaster(distance_to_seed,
            filename = "data/geodata/raster/burned_area/planet/distance_between_complete_and_unburned_aoi2.tif",
            overwrite = T
)
writeRaster(fraction_burned,
            filename = "data/geodata/raster/burned_area/planet/aoi2_burned_fraction.tif",
            overwrite = T
)
