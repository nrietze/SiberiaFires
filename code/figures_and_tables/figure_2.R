# Script to plot unburned patch metrics
# Nils Rietze: nils.rietze@uzh.ch
# 26 March 2024

library(terra)
library(tidyterra)
library(patchwork)
library(cowplot)
library(scales)
library(tidyverse)
library(landscapemetrics)
library(motif)
library(gt)
library(scico)
library(colorspace)
library(extrafont)
loadfonts(device = "win")

# 0. User functions ----
## a) Function to determine patch area ----
get_patch_area <- function(aois, aoi_name){
  
  # Get current AOI
  aoi_names_new <- c("DrainedThawlake" = "Ebelyakh",
                     "Kosukhino" = "Kosukhino",
                     "Berelech" = "Berelech",
                     "LargeScarCenter" = "Lapcha",
                     "LargeScarControl" = "Keremesit",
                     "Libenchik" = "Sala")
  
  aoi_name_new <- aoi_names_new[[aoi_name]]
  aoi_name_old <- aoi_name
  aoi <- aois[aois$site == aoi_name_new]
  
  cat(sprintf('Gathering patch areas for %s ... \n',aoi_name_new) )
  
  # Load burned area
  ba_path <- sprintf('data/geodata/raster/burned_area/planet/%s_burned_area_top5TD.tif',aoi_name_old)
  ba <- rast(ba_path)
  
  if (use_planet_wa){
    wa <- rast(
      sprintf('data/geodata/raster/water_area/planet/%s_water_area_top5TD.tif',aoi_name_old)
    )
    wa <- ifel(wa == 'water',wa,NA)
  } else {
    wa <- rast(
      sprintf('data/geodata/raster/water_area/%s_Landsat_mask.tif',aoi_name_old)
    )
    
    wa <- resample(wa,ba)
  }
  
  # Load burn perimeter
  bp <- vect(sprintf('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_%s.shp',aoi_name_old) ) %>%
    crop(aoi)
  
  # Mask out water pixels
  ba <- ba %>% 
    mask(wa, maskvalues = 2,updatevalue = NA) %>% 
    mask(bp,updatevalue = NA)
  
  # Load invalid data masks
  if (aoi_name_old %in% c('Berelech','LargeScarCenter')){
    poly_mask <- vect('data/geodata/feature_layers/planet_masks.shp') %>% 
      terra::intersect(aoi)
    
    ba <- mask(ba,poly_mask,inverse = TRUE)
    } 
  
  print(check_landscape(ba))
  
  # Compute patch areas
  patch_area <- lsm_p_area(ba)
  
  # Convert to square metres
  patch_area$area_m2 <- patch_area$value * 1e4
  
  patch_area$site <- aoi_name_new
  
  return(list(patch_area = patch_area, ba = ba))
}

## b) Function to apply on both burned and unburned patches ----
plot_patch_size_histogram <- function(df_patch_sizes, 
                                      CLASS_TO_PLOT,
                                      N_BREAKS = 30,
                                      FONT_SIZE = 18){
  
  # get histogram bins and values
  res <- df_patch_sizes %>%
    filter(class == CLASS_TO_PLOT) %>% 
    mutate(area_m2 = log(area_m2)) %>% 
    pull(area_m2) %>%
    hist(plot = FALSE,
         breaks = N_BREAKS)
  
  
  res_df <- data.frame(counts = res$counts,
                       mids = res$mids)
  
  x_res <- (res_df$mids[2] - res_df$mids[1])
  y_max <- max(res_df$counts, na.rm = T) * 1.1
  
  # prepare plotting parameters
  res_df <- mutate(res_df,
                   inv_perc = c(1, 1 - head(cumsum(counts)/sum(counts), -1)),
                   x_min = mids - (x_res / 2),
                   x_max = mids + (x_res / 2),
                   bar_lower = y_max * -0.2,
                   col = as.factor(1:nrow(res_df)))
  
  colour_ramp <- sequential_hcl(nrow(res_df), palette = 'inferno')
  
  # get % of patches that are large than 1 Landsat pixel
  landsat_area_log10 <- log10(30**2)
  
  ecdf_unburned <- df_patch_sizes %>%
    filter(class == CLASS_TO_PLOT) %>% 
    pull(area_m2) %>% 
    log10() %>%
    ecdf()
  
  # get % if we exclude unburned patches < 81 m2
  ecdf_unburned_81 <- df_patch_sizes %>%
    filter(class == CLASS_TO_PLOT) %>% 
    filter(area_m2 > 81) %>%
    pull(area_m2) %>%
    log10() %>%
    ecdf()
  cat(sprintf("%.2f %% of unburned patches are larger than Landsat
            (excluding single and 3x3-pixel patches)",
              100 - ecdf_unburned_81(landsat_area_log10) *100) )
  
  pct_landsat <- 100 - ecdf_unburned(landsat_area_log10)*100
  s <- bquote(.(round(pct_landsat,1))~"% of patches are larger than 900"~m^2)
  s <- sprintf("%.1f%% of patches \n are larger than 900 m²",round(pct_landsat,1))
  
  # plot with x-axis log10, y-axis counts 
  p_hist <- ggplot(res_df, aes(fill = col, colour = col)) +
    # histogram bars
    geom_rect(aes(xmin = x_min, xmax = x_max,
                  ymin = 0, ymax = counts),
              colour = "#505050", size = 0.005) +
    geom_rect(aes(xmin = x_min, xmax = x_max,
                  ymin = 0,
                  ymax = -.06 * max(counts)),
              size = 0.005) +
    # vertical line with text at 900m2
    geom_vline(xintercept = log(900),color = "grey40",
               linetype = "dashed", size = 2) +
    annotate("text", x = log(1e5), y = .4* max(res_df$counts),
             label = s,
             size = 7, color = "grey40") +
    geom_curve(aes(x = log(1e4), y = .3 * max(counts), 
                   xend = log(900), yend = .03 * max(counts)),
               arrow = arrow(length = unit(0.08, "inch")), linewidth = 1,
               color = "grey40", curvature = -0.3) +
    # black outline around colorramp
    annotate("rect",
             xmin = min(res$breaks),
             xmax = max(res_df$x_max),
             ymin = 0,
             ymax = -.06 * max(res_df$counts),
             fill = NA,
             colour = "black",
             size = 0.75) +
    # plot styling
    scale_fill_manual(values = colour_ramp) +
    scale_colour_manual(values = colour_ramp) +
    scale_x_continuous(
      expand = c(0,0),
      limits = c(min(res_df$x_min),max(res_df$x_max)),
      breaks = log(as.integer(10^(1:7))),
      labels = c(10,100,1000,"1e4","1e5","1e6","1e7")
    ) +
    labs(x = expression(Unburned~patch~size~(m^2)), 
         y = "Counts",
         # subtitle = "b)"
    ) +
    theme_cowplot(FONT_SIZE) +
    theme(legend.position = "none",
          plot.margin = ggplot2::margin(1, 1, 1, 0.25, "cm"),
          plot.subtitle=element_text(face='bold',size = FONT_SIZE))
  
  return(p_hist)
}

## c) Function to plot patch size histograms for all fire scars----
plot_all_histograms <- function(all_patch_areas,CLASS_TO_PLOT, FONT_SIZE = 18, color_ramp = "inferno"){
  if (CLASS_TO_PLOT == 1){
    label = expression(Unburned~patch~size~(m^2))
  } else {
    label = expression(Burned~patch~size~(m^2))
  }
  
  median_area_m2_per_site <- all_patch_areas %>%
    filter(class == CLASS_TO_PLOT) %>% 
    group_by(site) %>%
    summarise(md = median(area_m2, na.rm = TRUE))
  
  
  fig <- all_patch_areas %>% 
    filter(class == CLASS_TO_PLOT) %>% 
    ggplot(aes(x = area_m2, fill = site)) +
    geom_histogram(binwidth = 0.25, position = "identity") + 
    geom_vline(data = median_area_m2_per_site, aes(xintercept = md),
               linetype = "dashed", color = "gray") +
    scale_fill_viridis_d(option = color_ramp, end = .8) +
    scale_x_continuous(
      expand = c(-1e-2,0),trans='log10',limits = c(.9,1e5),
      labels = label_number()) +
    labs(x = label, y = "Counts") + 
    facet_wrap(~site,ncol = 3) +
    theme_cowplot(FONT_SIZE) + 
    theme(legend.position = 'none') 
  
  return(fig)
} 

## d) Function to plot patch size ECDFs for all fire scars ----
plot_all_ecdf <- function(all_patch_areas,CLASS_TO_PLOT, FONT_SIZE = 18, color_ramp = "inferno"){
  if (CLASS_TO_PLOT == 1){
    label = expression(Unburned~patch~size~(m^2))
  } else {
    label = expression(Burned~patch~size~(m^2))
  }
  
  fig <- all_patch_areas %>% 
    filter(class == CLASS_TO_PLOT) %>% 
    ggplot(aes(x = area_m2,colour = site)) +
    geom_vline(xintercept = 9,
               color = "grey40",linewidth = 1,linetype = "dashed") +
    annotate("text", x = 10, y = 30, 
             label = "1 PlanetScope pixel",
             color = "grey40",
             vjust = 0, hjust = 0, size = 5) +
    geom_vline(xintercept = 81, 
               color = "grey40",linewidth = 1,linetype = "dashed") +
    annotate("text", x = 85, y = 40, 
             label = "3 x 3 PlanetScope pixels", 
             color = "grey40",
             vjust = 0, hjust = 0, size = 5) +
    geom_vline(xintercept = 900,
               color = "grey40",linewidth = 1,linetype = "dashed") +
    annotate("text", x = 920, y = 60, 
             label = "1 Landsat pixel", 
             color = "grey40",
             vjust = 0, hjust = 0, size = 5) +
    stat_ecdf(aes(y = after_stat(y*100)),
              geom = "line",
              pad = FALSE, linewidth = 1) +
    scale_x_continuous(
      expand = c(0, 0),
      trans='log',
      limits = c(8,exp(15)),
      breaks = as.integer(10^(1:7)),
      labels = c(10,100,1000,"1e4","1e5","1e6","1e7")
      ) +
    scale_color_viridis_d(option = color_ramp, end = .8) +
    labs(x = label,
         y = "Cumulative Percentage (%)",color = 'Site') +
    theme_cowplot(FONT_SIZE) +
    theme(
      legend.position = "bottom",  
      legend.title = element_blank(),
      legend.direction = "horizontal",
    ) 
  
  return(fig)
}

## e) Function to load and mask burned area data ----
GetBurnedArea <- function(aoi_name, 
                          product = "descals", 
                          return_all = TRUE,
                          mask_lsat = FALSE){
  cat(sprintf("processing %s ... \n", aoi_name))
  
  # load burn perimeter
  bp <- vect(
    sprintf('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_%s.shp',aoi_name)
  )
  
  # Get current AOI
  aoi_names_new <- c("DrainedThawlake" = "Ebelyakh",
                     "Kosukhino" = "Kosukhino",
                     "Berelech" = "Berelech",
                     "LargeScarCenter" = "Lapcha",
                     "LargeScarControl" = "Keremesit",
                     "Libenchik" = "Sala")
  
  aoi_name_new <- aoi_names_new[[aoi_name]]
  aoi_name_old <- aoi_name
  aoi <- aois[aois$site == aoi_name_new]
  
  if (aoi_name %in% c("Berelech","LargeScarCenter")){
    has_mask <- TRUE
  } else{has_mask <- FALSE}
  
  # load water mask
  wa <- rast(
    sprintf('data/geodata/raster/water_area/%s_Landsat_mask.tif',aoi_name)
  )
  
  # load PlanetScope burned area
  ba <- rast(
    sprintf('data/geodata/raster/burned_area/planet/%s_burned_area_top5TD.tif',aoi_name)
  ) 
  wa_3m <- resample(wa,ba)
  
  # mask water areas in burned area map
  ba <- mask(ba,wa_3m,maskvalues = 2,updatevalue = NA) %>% 
    mask(bp)
  
  # load comparison Landsat burned areas
  if (product == "descals"){
    fn_comp <- 'data/geodata/raster/burned_area/ba_descals_landsat_2020_utm_shifted.tif'
    burn_val <- 30
    
  } else if (product == "GABAM"){
    fn_comp <- 'data/geodata/raster/burned_area/N75E145_burn_class_UTM_55N.tif'
    burn_val <- 255
  }
  
  ba_comp <- rast(fn_comp) %>% 
    crop(bp) %>%
    rename(class = 1)
  
  ba_comp_bin <- ifel(ba_comp == burn_val, 1, 0) %>%  #convert values to 0 & 1
    mask(ba_comp,maskvalues = NA, updatevalue = 0) %>% # convert NA values to 0 (unburned), important for GABAM
    mask(crop(wa,bp),maskvalues = 2, updatevalue = NA) %>%  # mask out water areas
    mask(bp) %>% 
    as.factor()
  levels(ba_comp_bin) <- data.frame(id=0:1, class= c('unburned','burned'))
  
  # load burned fractions
  r_preds <- rast(
    sprintf('data/geodata/raster/predictors/%s_predictors_30m.tif',aoi_name)
  ) %>% 
    # select(burned_fraction) %>% 
    crop(bp) %>% 
    mask(crop(wa,bp),maskvalues = 2, updatevalue = NA) %>% # mask out water areas
    mask(bp)
  
  if (has_mask){
    poly_mask <- vect('data/geodata/feature_layers/planet_masks.shp') %>% 
      terra::intersect(aoi) 
    
    ba_comp_bin <- ba_comp_bin %>% 
      mask(poly_mask,inverse = T)
    
    ba <- ba %>% 
      mask(poly_mask,inverse = T)
    
    r_preds <- r_preds %>% 
      mask(poly_mask,inverse = T)
  }
  
  if (mask_lsat){
    # resample Landsat burned area to 3 m (replicate pixels)
    ba_comp_bin_3m <- resample(ba_comp_bin,ba) %>% 
      as.int()
    
    # Mask Planetscope burned area with Landsat burned area to see what's new
    ba <- terra::mask(ba, ba_comp_bin_3m,
                      maskvalues = 1, updatevalue = NA)
  }
  
  # 2. Extract and format data
  df_bf_in_comp <- rbind(
    # extract burned fractions from Landsat 'burned' (Descals et al. 2022)
    comp_burned <- ifel(ba_comp_bin == 'burned',r_preds$burned_fraction,NA) %>% 
      as.data.frame() %>% 
      rename(burned_fraction = 1) %>%
      mutate(group = "burned"),
    # Extract burned fractions from Landsat 'unburned' (Descals et al. 2022)
    comp_unburned <- ifel(ba_comp_bin == 'unburned',r_preds$burned_fraction,NA) %>% 
      as.data.frame() %>% 
      rename(burned_fraction = 1) %>%
      mutate(group = "unburned")
  ) %>% 
    mutate(group = as.factor(group))
  
  if (return_all){
    return(list(crop(ba,bp), 
                r_preds,
                ba_comp_bin,
                df_bf_in_comp,
                bp))
  } else{
    return(df_bf_in_comp)
  }
  
}

# 1. Load data ----
plot_all_aois <- TRUE
use_planet_wa <- FALSE

# Load areas of interest
aois <- vect('./data/geodata/feature_layers/aoi_wv/aois_analysis.geojson') %>%
  project('EPSG:32655') %>% 
  mutate(., id = 1:nrow(.))

# set which fire scars should be analyzed (all or just one?)
if (plot_all_aois){
  aoi_names <- c("Berelech",
                 "LargeScarCenter",
                 "LargeScarControl",
                 "Libenchik",
                 "Kosukhino",
                 "DrainedThawlake"
  )
} else {
  aoi_names <- c('LargeScarCenter')
}

CLASS_TO_PLOT <- 1 # ( 1 = unburned, 2 = burned)
FONT_SIZE  <- 18

# Initialize lists to store results
all_patch_areas <- list()
all_ba_objects <- list()

# Loop over AOIs
aoi_names_new <- c("DrainedThawlake" = "Ebelyakh",
                   "Kosukhino" = "Kosukhino",
                   "Berelech" = "Berelech",
                   "LargeScarCenter" = "Lapcha",
                   "LargeScarControl" = "Keremesit",
                   "Libenchik" = "Sala")

# Calculate patch areas of burned and unburned patches in all fire scars
for (aoi_name in aoi_names) {
  result <- get_patch_area(aois = aois, aoi_name = aoi_name)
  
  aoi_name_new <- aoi_names_new[[aoi_name]]
  all_patch_areas[[aoi_name_new]] <- result$patch_area
  all_ba_objects[[aoi_name_new]] <- result$ba
}

# Combine patch areas into one dataframe if needed
all_patch_areas <- do.call(rbind, all_patch_areas) %>% 
  mutate(site = factor(site,levels = aoi_names_new))

# 2. Create Figure 2 ----

# Extract unburned patch sizes in the Kosukhino scar
patch_area_rast_kosuk <- spatialize_lsm(ifel(all_ba_objects$Kosukhino == CLASS_TO_PLOT,
                                             all_ba_objects$Kosukhino,NA),
                                        what = "lsm_p_area",
                                        directions = 8)
# convert to log m2 (from ha)
raster_kosukh <- log(patch_area_rast_kosuk$layer_1$lsm_p_area * 1e4)

# convert to dataframe
raster_df <- raster_kosukh %>% 
  rename(x = value) %>% 
  as.data.frame()

# Load geodata for Map of unburned patches
bp_kosukh <- vect('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_Kosukhino.shp')
wa <- rast('data/geodata/raster/water_area/Kosukhino_Landsat_mask.tif')

# resample Landsat water mask to PlanetScope resolution
wa <- resample(wa,raster_kosukh)

# Define position of text and drone symbol
xtext <- 544500
ytext <- 7871400
textsize <- 7
  
drone_lat <- 7867419
drone_lon <- 545492
drone_angle <- 114

# Plot map of unburned islands
fig_2a <- ggplot() +
    # plot unburned areas
    geom_spatraster(data = raster_kosukh,aes(fill = ..value..),
                    show.legend = FALSE) +
    scale_fill_viridis_c(option = 'inferno', na.value = "white",
                         guide = guide_colorbar(
                           title = expression(Patch~area~(m^2)),
                           ),
                         breaks = log(as.integer(10^(1:6))),
                         labels = as.integer(10^(1:6))
                         ) +
    # plot water areas
    ggnewscale::new_scale_fill() +
    geom_spatraster(data = wa,inherit.aes = FALSE,
                    aes(fill = ..value..),alpha = 0.4,
                    show.legend = FALSE) +
    scale_fill_gradientn(colours = c("transparent", "steelblue4"),
                         na.value = "transparent", values = c(0, 1)) +
    # plot burn perimeter
    geom_spatvector(data = bp_kosukh,
                    fill = NA, linewidth = 1,
                    color = 'grey40',
                    show.legend = FALSE) +
    # scale bar
    geom_rect(aes(xmin = xtext - 2500, xmax = xtext - 500, ymin = ytext, ymax = ytext + 200),
              fill = 'gray20') + 
    geom_text(aes(x = xtext - 1500, y = ytext + 600,
                  label = '2 km',fontface = 'bold'),
              size = textsize,
              colour = 'gray20') +
    # north arrow
    geom_text(aes(x = xtext + 40, y = ytext + 100,
                  label = 'N',fontface = 'bold'),
              size = textsize, colour = 'gray20') +
    geom_text(aes(x = xtext + 65, y = ytext + 560,
                  label = '$', angle = 90,
                  family = 'ESRI arrowhead',
                  fontface = 'bold'),
              size = textsize, colour = 'gray20') +
    # drone marker
    geom_text(aes(x = drone_lon, y = drone_lat,
                label = 'D', angle = drone_angle,
                family = 'ESRI arrowhead'),
            size = textsize, colour = 'gray20') +
    # water text
    geom_text(aes(x = xtext - 2000, y = ytext - 900,
                  label = 'w a t e r',fontface = 'bold'),
              size = textsize - .2, angle = 20,
              colour = 'steelblue4') +
    # plot styling
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_map(FONT_SIZE) +
    theme(plot.subtitle=element_text(face='bold',size = FONT_SIZE),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) 

# Plot Histograms of unburned patches
df_patch_sizes <- all_patch_areas %>% 
  filter(site == 'Kosukhino') 

fig_2b <- plot_patch_size_histogram(df_patch_sizes = df_patch_sizes,
                                CLASS_TO_PLOT = 1)

# Plot drone image
fig_2c <- ggdraw() +
  # draw_image("figures/Figure_2c.png", width = 1) +
  draw_image("figures/Figure_2c_square.PNG", width = 1) + # square version of drone image
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

# Plot ECDFs of unburned patch sizes for all fire scars
fig_2d <- plot_all_ecdf(all_patch_areas, 1, color_ramp = "mako")

# create plot_grid
pp <- fig_2a + fig_2b + fig_2c + fig_2d + 
  plot_layout(ncol = 2,widths = c(1, 1)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") &
  theme(plot.tag = element_text(size = FONT_SIZE))

# add annotation for drone location
final_plot <- ggdraw(pp) + 
  draw_plot_label(
    label = "D", 
    x = 0.06, y = 0.55, angle = drone_angle,
    hjust = 0.5, vjust = 0.5,
    family = "ESRI arrowhead", size = FONT_SIZE, color = "gray20"
  ) +
  draw_plot_label(
    label = "Drone viewshed in (c) (southeast)", 
    x = 0.07, y = 0.55,
    hjust = 0, vjust = 0.5, 
    fontface = "bold", size = FONT_SIZE, color = "gray20"
  )

# export figure
ggsave2(final_plot, filename ='figures/Figure_2.png',
        device = png, type = "cairo",
        bg = 'white',width = 16, height = 16)

# 3. Create Figure s5: Patch size histograms for all fire scars ----
fig_s5a <- plot_all_histograms(all_patch_areas, 1, color_ramp = "mako") # unburned patch histograms
fig_s5b <- plot_all_histograms(all_patch_areas, 2, color_ramp = "mako") # burned patch histograms

# create plot grid and export
fig_s5 <- fig_s5a + fig_s5b + 
  plot_layout(nrow = 2,heights = c(1,1)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") &
  theme(plot.tag = element_text(size = FONT_SIZE))

ggsave(fig_s5, filename ='figures/Figure_S5.png',
       bg = 'white',width = 12, height = 8)

# 4. Create Figure S6: Patch size ECDFs for all fire scars ----
fig_s6 <- plot_all_ecdf(all_patch_areas, 2, color_ramp = "mako") # burned patch ECDFs

# export 
ggsave(fig_s6, filename ='figures/Figure_S6.png',
       bg = 'white',width = 10, height = 10)

# 5. Create Table s11 & s12: Table of patch size statistics ----

# Calculate largest patch index
df_lpi <- data.frame()
for (aoi_name_new in aoi_names_new){
  lpi_result <- lsm_c_lpi(all_ba_objects[[aoi_name_new]])
  
  lpi_result <- lpi_result %>%
    mutate(site = aoi_name_new) %>% 
    select(-c(layer, level, id, metric)) %>% 
    rename(c("LPI" = "value"))
  
  df_lpi <- bind_rows(df_lpi, lpi_result) # add to dataframe
}

# Calculate statistics per burn class and site
mean_area_m2_per_site <- all_patch_areas %>%
  group_by(class,site) %>%
  summarise(md = median(area_m2, na.rm = TRUE),
            mdc = median(area_m2[area_m2 > 81], na.rm = TRUE),
            mu = mean(area_m2, na.rm = TRUE),
            area = sum(area_m2) / 1e6) %>% 
  pivot_wider(names_from = class, values_from = c(md,mdc,mu, area)) %>% 
  pivot_longer(cols = c(-site),
               names_to = c(".value", "class"),
               names_sep = "_") 

mean_area_m2_per_site <- df_lpi %>% 
  mutate(class = as.character(class),
         site = as.factor(site)) %>% 
  left_join(mean_area_m2_per_site, by = c("class","site")) %>% 
  select(class, site, md, mdc, mu, LPI, area)

# Reformat table and format to nice looking table
t_s11 <- mean_area_m2_per_site %>%
  filter(class != 2) %>% select(-class) %>% 
  mutate(site = str_replace_all(site, "(?<=[a-z])(?=[A-Z])", " ")) %>%
  gt() %>% 
  cols_label(
    site = "Site",
    md = html("Median patch size <br>m<sup>2</sup>"),
    mdc = html("Conservative median patch size <br>m<sup>2</sup>"),
    mu = html("Mean patch size <br>m<sup>2</sup>"),
    LPI = html("LPI <br> %"),
    area = html("Total unburned area <br>km<sup>2</sup>"),
  ) %>% 
  fmt_number(decimals = 1) %>%
  tab_style(
    style = cell_text(v_align = "top", weight = 'bold'),
    locations = cells_column_labels())

# export
gtsave(t_s11,filename = "tables/TableS11.html")

# Reformat table and format to nice looking table
t_s12 <- mean_area_m2_per_site %>%
  filter(class != 1) %>% select(-class) %>% 
  mutate(site = str_replace_all(site, "(?<=[a-z])(?=[A-Z])", " ")) %>%
  gt() %>% 
  cols_label(
    site = "Site",
    md = html("Median patch size <br>m<sup>2</sup>"),
    mdc = html("Conservative median patch size <br>m<sup>2</sup>"),
    mu = html("Mean patch size <br>m<sup>2</sup>"),
    LPI = html("LPI <br> %"),
    area = html("Total burned area <br>km<sup>2</sup>"),
  ) %>% 
  fmt_number(decimals = 1) %>%
  tab_style(
    style = cell_text(v_align = "top", weight = 'bold'),
    locations = cells_column_labels())

# export
gtsave(t_s12,filename = "tables/TableS12.html")

# 6. Calculate statistics of newly detected burned patches ----
all_ba_objects <- lapply(aoi_names, function(aoi_name) {
  L <- GetBurnedArea(aoi_name, mask_lsat = TRUE)
  L[[1]]  # Extract the first item from the list
})

names(all_ba_objects) <- aoi_names

# Calculate LPI for all patches in Landsat unburned
df_lpi <- data.frame()
for (aoi_name_new in aoi_names_new){
  aoi_name_old <- names(aoi_names_new[aoi_names_new == aoi_name_new])
  lpi_result <- lsm_c_lpi(all_ba_objects[[aoi_name_old]])
  
  lpi_result <- lpi_result %>%
    mutate(site = aoi_name_new) %>% 
    select(-c(layer, level, id, metric)) %>% 
    rename(c("LPI" = "value"))
  
  df_lpi <- bind_rows(df_lpi, lpi_result) # add to dataframe
}

plot(all_ba_objects[["Kosukhino"]])

# Calculate patch size statistics per burn class and site (Table S13)
df_pa <- data.frame()

for (aoi_name_new in aoi_names_new) {
  aoi_name_old <- names(aoi_names_new[aoi_names_new == aoi_name_new])
  patch_area <- lsm_p_area(all_ba_objects[[aoi_name_old]]) %>% 
    mutate(area_m2 = value * 1e4) %>% 
    group_by(class) %>%
    summarise(md = median(area_m2, na.rm = TRUE),
              mdc = median(area_m2[area_m2 > 81], na.rm = TRUE),
              mu = mean(area_m2, na.rm = TRUE),
              area = sum(area_m2) / 1e6) %>% 
    pivot_wider(names_from = class, values_from = c(md,mdc,mu, area)) %>% 
    pivot_longer(cols = everything(),
                 names_to = c(".value", "class"),
                 names_sep = "_")
  
  patch_area$site <- aoi_name_new
  
  df_pa <- bind_rows(df_pa, patch_area) # add to dataframe
}

(t_s99 <- df_lpi %>% 
  mutate(class = as.character(class),
         site = as.factor(site)) %>% 
  left_join(df_pa, by = c("class","site")) %>% 
  mutate(class = if_else(class == 1, "unburned","burned")) %>% 
  select(class, site, md, mdc, mu, LPI, area) %>% 
  # filter(class != 1) %>% select(-class) %>% 
  mutate(site = str_replace_all(site, "(?<=[a-z])(?=[A-Z])", " ")) %>%
  gt() %>% 
  cols_label(
    class = "Burn class",
    site = "Fire scar",
    md = html("Median patch size <br>m<sup>2</sup>"),
    mdc = html("Conservative median patch size <br>m<sup>2</sup>"),
    mu = html("Mean patch size <br>m<sup>2</sup>"),
    LPI = html("LPI <br> %"),
    area = html("Total area <br>km<sup>2</sup>"),
  ) %>% 
  fmt_number(decimals = 1) %>%
  tab_style(
    style = cell_text(v_align = "top", weight = 'bold'),
    locations = cells_column_labels()))

gtsave(t_s99,filename = "tables/TableS13.html")
