# Script to plot unburned patch metrics
# Nils Rietze: nils.rietze@uzh.ch 
# 26 March 2024

library(terra)
library(tidyterra)
library(cowplot)
library(scales)
library(tidyverse)
library(landscapemetrics)
library(motif)
library(scico)
library(colorspace)

# 1. Configuration & loading data ----
plot_all_aois <- TRUE

# Load areas of interest
aois <- vect('./data/geodata/feature_layers/aoi_wv/aois_analysis.geojson') %>%
  project('EPSG:32655') %>% 
  mutate(., id = 1:nrow(.))

# Provide a part of the AOI name for processing
if (plot_all_aois){
  aoi_names <- aois$site
} else {
  aoi_names <- c('LargeScarCenter')
}

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
  bp <- vect(sprintf('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_%s.shp',aoi$site) ) %>%
    crop(aoi)
  
  # Mask out water pixels
  ba <- ba %>% 
    mask(wa, maskvalues = 2,updatevalue = NA) %>% 
    mask(bp,updatevalue = NA)
  
  # Load invalid data masks
  if (aoi$site %in% c('Berelech','LargeScarCenter')){
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
  
  return(list(patch_area = patch_area, ba = ba))
}

# Initialize lists to store results
all_patch_areas <- list()
all_ba_objects <- list()

# Loop over AOIs
for (aoi_name in aoi_names) {
  result <- get_patch_area(aois = aois, aoi_name = aoi_name)
  all_patch_areas[[aoi_name]] <- result$patch_area
  all_ba_objects[[aoi_name]] <- result$ba
}

# Combine patch areas into one dataframe if needed
all_patch_areas <- do.call(rbind, all_patch_areas) %>% 
  mutate(site = factor(site))

# 2. Plot patch area proportions ----
class_to_plot <- 1 # ( 1 = unburned, 2 = burned)

## a) Plot Histograms of unburned islands ----
n_steps <- 15

# For AOI Kosukhino
patch_area_rast_kosuk <- spatialize_lsm(ifel(all_ba_objects$Kosukhino == class_to_plot,
                                             all_ba_objects$Kosukhino,NA),
                                        what = "lsm_p_area",
                                        directions = 8)
raster_kosukh <- log10(patch_area_rast_kosuk$layer_1$lsm_p_area * 1e4)

raster_df <- raster_kosukh %>% 
  rename(x = value) %>% 
  as.data.frame()

bp_kosukh <- vect('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_Kosukhino.shp')
wa <- rast('data/geodata/raster/water_area/planet/Kosukhino_water_area_top5TD.tif')
wa <- ifel(wa == 'water',wa,NA)

(p1 <- ggplot() +
    geom_spatraster(data = raster_kosukh,aes(fill = ..value..)) +
    scale_fill_viridis_c(option = 'inferno', na.value = "white") +
    theme(legend.position = 'none') +
    ggnewscale::new_scale_fill() +
    geom_spatraster(data = wa,inherit.aes = FALSE,
                    aes(fill = ..value..)) +
    scale_fill_gradientn(colours = c("transparent", "steelblue4"),
                         na.value = "transparent", values = c(0, 1)) +
    geom_spatvector(data = bp_kosukh,
                    fill = NA, linewidth = 1,
                    color = 'grey80') +
    theme(legend.position = 'none') +
    theme_cowplot())

# Plot histogram
res <- all_patch_areas %>% 
  filter(site == 'Kosukhino') %>% 
  mutate(area_m2 = log10(area_m2)) %>%  
  pull(area_m2) %>% 
  hist(plot = FALSE)

res_df <- data.frame(counts = res$counts,
                     mids = res$mids)
ggplot(res_df) +
  geom_bar(aes(x = mids,y = counts,fill = mids),
           stat = 'identity') +
  geom_rect(aes(xmin = min(res$breaks), xmax = max(res$breaks),
                ymin = -max(res$counts)*.2, ymax = 0,
                fill = mids),
            size = 1) +
  annotate("rect", 
           xmin = min(res$breaks),
           xmax = max(res$breaks),
           ymin = -max(res$counts)*.2,
           ymax = 0,
           fill = NA, 
           colour = "black",
           size = 0.75) +
  scale_fill_viridis_c(option = 'inferno') +
  scale_color_viridis_c(option = 'inferno') +
  scale_x_continuous(
    expand = c(-1e-2,0),
    breaks = c(1:5),  # Set the breaks to the desired values
    labels = c("10", "100", "1000", "10000", "100000") ) +
  labs(x = expression(Patch~Area~(m^2)), y = "Frequency") +
  theme_cowplot(font_size) +
  theme(legend.position = "bottom",
        plot.margin = margin(1, 1, 1, 0.25, "inch"),
        axis.title.y = element_text(hjust = 0.7))

# For all AOIs
mean_area_m2_per_site <- all_patch_areas %>%
  filter(class == class_to_plot & area_m2 > 81) %>%
  group_by(site) %>%
  summarise(mean_area_m2 = median(area_m2, na.rm = TRUE))

all_patch_areas %>% 
  filter(class == class_to_plot) %>% 
  ggplot(aes(x = area_m2, fill = site)) +
  geom_histogram(binwidth = 0.25, alpha = 0.5, 
                 position = "identity") + 
  geom_vline(data = mean_area_m2_per_site, aes(xintercept = mean_area_m2),
             linetype = "dashed", color = "gray") +
  scale_fill_viridis_d(option = 'inferno') +
  scale_x_continuous(
    expand = c(-1e-2,0),
    trans='log10',
    limits = c(.9,1e5),
    labels = label_number()) +
  labs(x = expression(Patch~Area~(m^2)), y = "Frequency") + 
  facet_wrap(~site,ncol = 3) +
  theme_cowplot() + 
  theme(legend.position = 'none')

## b) Plot ECDF ----
(patch_plot <- all_patch_areas %>% 
    filter(class == class_to_plot) %>% 
    # group_by(site) %>%
    # arrange(desc(area_m2)) %>%
    # slice(-(1:2)) %>%
    # filter(area_m2 != max(area_m2)) %>% 
    # ungroup() %>% 
    ggplot(aes(x = area_m2,colour = site)) +
    geom_vline(xintercept = 9, linewidth = .5,linetype = "dashed") +
    annotate("text", x = 10, y = 30, label = "1 planet pixel", vjust = 0, hjust = 0, size = 3) +
    geom_vline(xintercept = 81, linewidth = .5,linetype = "dashed") +
    annotate("text", x = 85, y = 40, label = "3 x 3 planet pixels", vjust = 0, hjust = 0, size = 3) +
    geom_vline(xintercept = 900, linewidth = .5,linetype = "dashed") +
    annotate("text", x = 920, y = 60, label = "1 landsat pixel", vjust = 0, hjust = 0, size = 3) +
    stat_ecdf(aes(y = after_stat(y*100)),
              geom = "line",
              pad = FALSE, linewidth = 1) +
    scale_x_continuous(expand = c(1e-2, 0),trans='log10',
                       limits = c(8,1e5),
                       labels = label_number()) +
    # scale_color_manual(values = site_colors) +
    scale_color_viridis_d(option = "inferno") +
    labs(x = expression(Patch~area~(m^2)),
         y = "Cumulative Percentage (%)",color = 'Site') +
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

# Patch distance from unburned to burned
aoi_name <- 'DrainedThawlake'
bp <- vect(sprintf('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_%s.shp',aoi_name) ) %>%
  crop(aois[aois$site == aoi_name])
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
