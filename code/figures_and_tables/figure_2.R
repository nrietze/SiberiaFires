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
library(extrafont)
loadfonts(device = "win")

# 1. Configuration & loading data ----
plot_all_aois <- TRUE
use_planet_wa <- FALSE

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
  
  if (use_planet_wa){
    wa <- rast(
      sprintf('data/geodata/raster/water_area/planet/%s_water_area_top5TD.tif',aoi_name)
    )
    wa <- ifel(wa == 'water',wa,NA)
  } else {
    wa_ps <- rast(
      sprintf('data/geodata/raster/water_area/planet/%s_water_area_top5TD.tif',aoi_name)
    )
    
    wa <- rast(
      sprintf('data/geodata/raster/water_area/%s_Landsat_mask.tif',aoi_name)
    )
    
    wa <- resample(wa,wa_ps)
    rm(wa_ps)
  }
  
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
  mutate(site = factor(site,levels = aoi_names))

# 2. Plot patch area proportions ----
class_to_plot <- 1 # ( 1 = unburned, 2 = burned)
font_size  <- 18

## a) Plot map of unburned islands ----
n_steps <- 15

# For AOI Kosukhino
patch_area_rast_kosuk <- spatialize_lsm(ifel(all_ba_objects$Kosukhino == class_to_plot,
                                             all_ba_objects$Kosukhino,NA),
                                        what = "lsm_p_area",
                                        directions = 8)
raster_kosukh <- log(patch_area_rast_kosuk$layer_1$lsm_p_area * 1e4)

raster_df <- raster_kosukh %>% 
  rename(x = value) %>% 
  as.data.frame()

bp_kosukh <- vect('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_Kosukhino.shp')
wa_ps <- rast('data/geodata/raster/water_area/planet/Kosukhino_water_area_top5TD.tif')
wa <- rast('data/geodata/raster/water_area/Kosukhino_Landsat_mask.tif')

wa <- resample(wa,wa_ps)
rm(wa_ps)

xtext <- 544500
ytext <- 7871400
textsize <- 7
  
drone_lat <- 7867419
drone_lon <- 545492
drone_angle <- 114
  
p1 <- ggplot() +
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
    theme_map(font_size) +
    theme(plot.subtitle=element_text(face='bold',size = font_size),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) 

## b) Plot Histograms of unburned islands ----
n_breaks <- 30

# get histogram bins and values
res <- all_patch_areas %>% 
  filter(site == 'Kosukhino') %>% 
  mutate(area_m2 = log(area_m2)) %>%  
  pull(area_m2) %>% 
  hist(plot = FALSE,
       breaks = n_breaks)

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

ecdf_unburned <- all_patch_areas %>%
  filter(class == 1 & site  == 'Kosukhino') %>%
  pull(area_m2) %>%
  log10() %>%
  ecdf()

# get % if we exclude unburned patches <= 81 m2
ecdf_unburned_81 <- all_patch_areas %>%
  filter(class == 1 & site  == 'Kosukhino' & area_m2 > 81) %>%
  pull(area_m2) %>%
  log10() %>%
  ecdf()
cat("% of unburned patches larger than Landsat 
      excluding single and 3x3-pixel patches:",
      100 - ecdf_unburned_81(landsat_area_log10) *100)

all_patch_areas %>% 
  filter(site == 'Kosukhino') %>%
  arrange(desc(area_m2)) %>%
  mutate(
    top_10 = row_number() <= 10,
    sum_top_10 = sum(area_m2[top_10]),
    sum_remaining = sum(area_m2[!top_10])
  ) %>%
  summarise(
    sum_top_10 = first(sum_top_10)/1e6,
    sum_remaining = first(sum_remaining)/1e6,
    percent = sum_remaining / (sum_top_10 + sum_remaining) *100
  )

pct_landsat <- 100 - ecdf_unburned(landsat_area_log10)*100
s <- bquote(.(round(pct_landsat,1))~"% of patches are larger than 900"~m^2)
s <- sprintf("%.1f%% of patches \n are larger than 900 m²",round(pct_landsat,1))

# plot x-axis log10, y-axis inverse cumulative percentage
(p2 <- ggplot(res_df, aes(fill = col, colour = col)) +
    # histogram bars
    geom_rect(aes(xmin = x_min, xmax = x_max,
                  ymin = 0, ymax = inv_perc),
              colour = "#505050", size = 0.005) +
    geom_rect(aes(xmin = x_min, xmax = x_max,
                  ymin = 0,
                  ymax = -.06),
              size = 0.005) +
    # vertical line with text at 900m2
    geom_vline(xintercept = log(900),color = "grey40",
               linetype = "dashed", size = 2) +
    annotate("text", x = log(1e5), y = .4,
             label = s,
             size = 7, color = "grey40") +
    geom_curve(aes(x = log(1e4), y = .35, xend = log(900), yend = .03),
      arrow = arrow(length = unit(0.08, "inch")), linewidth = 1,
      color = "grey40", curvature = -0.3) +
    # black outline around colorramp
    annotate("rect",
             xmin = min(res$breaks),
             xmax = max(res_df$x_max),
             ymin = 0,
             ymax = -.06,
             fill = NA,
             colour = "black",
             size = 0.75) +
    # plot styling
    scale_fill_manual(values = colour_ramp) +
    scale_colour_manual(values = colour_ramp) +
    scale_y_continuous(labels = label_percent()) +
    scale_x_continuous(
    expand = c(0,0),
    limits = c(min(res_df$x_min),max(res_df$x_max)),
    breaks = log(as.integer(10^(1:7))),
    labels = c(10,100,1000,"1e4","1e5","1e6","1e7")
  ) +
  labs(x = expression(Unburned~patch~size~(m^2)), 
       y = "Inverse cumulative percentage",
       # subtitle = "b)"
       ) +
  theme_cowplot(font_size) +
  theme(legend.position = "none",
        plot.margin = ggplot2::margin(1, 1, 1, 0.25, "cm"),
        plot.subtitle=element_text(face='bold',size = font_size),
        axis.title.y = element_text(hjust = 0.7))
)

### plot x-axis log10, y-axis counts -----
(p2 <- ggplot(res_df, aes(fill = col, colour = col)) +
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
    theme_cowplot(font_size) +
    theme(legend.position = "none",
          plot.margin = ggplot2::margin(1, 1, 1, 0.25, "cm"),
          plot.subtitle=element_text(face='bold',size = font_size))
)

## c) create plot_grid and export ----
top_row <- cowplot::plot_grid(p1,p2,
                         ncol = 2,
                         labels = c("a)","b)"),label_size = font_size,
                         align = 'h', axis = 'tb',
                         rel_widths = c(1.5,1.8))

p_pic <- ggdraw() +
  draw_image("figures/Figure_2c.png", width = 1) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

pg <- cowplot::plot_grid(top_row,p_pic,
                         labels = c('','c)'),label_size = font_size,
                         ncol = 1)  

pg1 <- pg +
  geom_text(aes(x = 0.03, y = 0.6,
                label = 'D', angle = drone_angle,
                family = 'ESRI arrowhead'),
            size = textsize, colour = 'gray20') +
  geom_text(aes(x = 0.05, y = 0.6,hjust = "left",
                label = 'Drone position (looking southeast)',fontface = 'bold'),
            size = textsize, colour = 'gray20')


ggsave2(pg1, filename ='figures/Figure_2_nn.png',
        device = png, type = "cairo",
        bg = 'white',width = 12, height = 16)

# export individual figures
ggsave(p1,filename ='figures/map_kosukh_patch_areas.png',
       bg = 'white',width = 10, height = 10)
ggsave(p2,filename ='figures/inv_cum_pct_kosukh_patch_areas.png',
       bg = 'white',width = 10, height = 6)

# 3. Area statistics for all sites ----

## a) create table of all statistics ----

# total burned areas per site
bd <- c(22.6,11.3,0.9,23.8,3.5,8.5)

# Get mean and median patch areas per site
mean_area_m2_per_site <- all_patch_areas %>%
  group_by(class,site) %>%
  summarise(md = median(area_m2, na.rm = TRUE),
            mdc = median(area_m2[area_m2 > 81], na.rm = TRUE),
            mu = mean(area_m2, na.rm = TRUE),
            area = sum(area_m2) / 1e6) %>% 
  pivot_wider(names_from = class, values_from = c(md,mdc,mu, area)) %>% 
  mutate(ratio_ubd_bd = area_1  / area_2) %>% 
  pivot_longer(cols = c(-site, -ratio_ubd_bd),
               names_to = c(".value", "class"),
               names_sep = "_") %>% 
  relocate(ratio_ubd_bd, .after = last_col()) %>% 
  filter(class != 2) %>% select(-class)
  
# Reformat table and format to nice looking table
mean_area_m2_per_site %>%
  mutate(site = str_replace_all(site, "(?<=[a-z])(?=[A-Z])", " ")) %>%
  gt() %>% 
  cols_label(
    site = "Site",
    md = html("Median area <br>m<sup>2</sup>"),
    mdc = html("Conservative median area <br>m<sup>2</sup>"),
    mu = html("Mean area <br>m<sup>2</sup>"),
    area = html("Total unburned area <br>km<sup>2</sup>"),
    ratio_ubd_bd = html("Ratio <br> Unburned:Burned")
  ) %>% 
  fmt_number(decimals = 1) %>%
  tab_style(
    style = cell_text(v_align = "top", weight = 'bold'),
    locations = cells_column_labels()) %>% 
  gtsave(filename = "tables/Table4.html")

## b) Plot all site's patch area distributions ----
### i. separate histogram ----
(ps1 <- all_patch_areas %>% 
  filter(class == class_to_plot) %>% 
  ggplot(aes(x = area_m2, fill = site)) +
  geom_histogram(binwidth = 0.25, 
                 position = "identity") + 
  geom_vline(data = mean_area_m2_per_site, aes(xintercept = md),
             linetype = "dashed", color = "gray") +
  scale_fill_viridis_d(option = 'inferno', end = .8) +
  scale_x_continuous(
    expand = c(-1e-2,0),
    trans='log10',
    limits = c(.9,1e5),
    labels = label_number()) +
  labs(x = expression(Patch~Area~(m^2)), y = "Frequency") + 
  facet_wrap(~site,ncol = 3) +
  theme_cowplot() + 
  theme(legend.position = 'none') )

### ii. combined ECDFs ----
(ps2 <- all_patch_areas %>% 
    filter(class == class_to_plot) %>% 
    # group_by(site) %>%
    # arrange(desc(area_m2)) %>%
    # slice(-(1:2)) %>%
    # filter(area_m2 != max(area_m2)) %>% 
    # ungroup() %>% 
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
    scale_x_continuous(expand = c(1e-2, 0),trans='log10',
                       limits = c(8,1e5),
                       labels = label_number()) +
    # scale_color_manual(values = site_colors) +
    scale_color_viridis_d(option = "inferno", end = .8) +
    labs(x = expression(Patch~area~(m^2)),
         y = "Cumulative Percentage (%)",color = 'Site') +
    theme_cowplot() +
    theme(
      legend.position = "bottom",  
      legend.title = element_blank(),
      legend.direction = "horizontal",
    ) 
  )

### iii. create plot grid and export ----
(pgs <- cowplot::plot_grid(ps1,ps2,
                          labels = c("a)","b)"),
                          nrow = 2,
                          align = 'h', axis = 'tb',
                          rel_heights = c(1.5,1.8)) )
ggsave2(pgs, filename = "figures/Figure_S1.png",
        bg = "white",width = 8,height = 14)

# export individual plots
ggsave(ps1, filename ='figures/Figure_S1.png',
       bg = 'white',width = 14, height = 10)
ggsave(ps2,filename ='figures/ECDF_all_patch_areas.png',
       bg = 'white',width = 10, height = 10)