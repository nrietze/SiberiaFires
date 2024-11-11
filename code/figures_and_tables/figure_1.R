# Script to plot binary burned area maps and burned fractions
# Nils Rietze: nils.rietze@uzh.ch 
# 26 March 2024

library(terra)
library(tidyterra)
library(cowplot)
library(patchwork)
library(scales)
library(tidyverse)
library(scico)
library(colorspace)
library(gt)
library(extrafont)
loadfonts(device = "win")

# load geom_flat_violin
source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")

# 1. Configure and load stuff ----
set.seed(1234)

use_planet_wa <- FALSE    # use Planet water mask (=TRUE) or Landsat's QA mask (=FALSE)?

# Load study regions
aois <- vect('./data/geodata/feature_layers/aoi_wv/aois_analysis.geojson') %>%
  project('EPSG:32655') %>% 
  mutate(., id = 1:nrow(.))

poly_mask <- vect('data/geodata/feature_layers/planet_masks.shp') 

GetBurnedArea <- function(aoi_name, product = "descals", return_all = TRUE){
  cat(sprintf("processing %s ... \n", aoi_name))
  
  # load burn perimeter
  bp <- vect(
    sprintf('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_%s.shp',aoi_name)
  )
  
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
  fraction_burned <- rast(
    sprintf('data/geodata/raster/predictors/%s_predictors_30m.tif',aoi_name)
  ) %>% 
    select(burned_fraction) %>% 
    crop(bp) %>% 
    mask(crop(wa,bp),maskvalues = 2, updatevalue = NA) %>% # mask out water areas
    mask(bp)
  
  if (has_mask){
    ba_comp_bin <- ba_comp_bin %>% 
      mask(poly_mask,inverse = T)
    
    ba <- ba %>% 
      mask(poly_mask,inverse = T)
    
    fraction_burned <- fraction_burned %>% 
      mask(poly_mask,inverse = T)
  }
  
  # 2. Extract and format data ----
  df_bf_in_comp <- rbind(
    # extract burned fractions from Landsat 'burned' (Descals et al. 2022)
    comp_burned <- ifel(ba_comp_bin == 'burned',fraction_burned,NA) %>% 
      as.data.frame() %>% 
      rename(burned_fraction = 1) %>%
      mutate(group = "burned"),
    # Extract burned fractions from Landsat 'unburned' (Descals et al. 2022)
    comp_unburned <- ifel(ba_comp_bin == 'unburned',fraction_burned,NA) %>% 
      as.data.frame() %>% 
      rename(burned_fraction = 1) %>%
      mutate(group = "unburned")
  ) %>% 
    mutate(group = as.factor(group))
  
  if (return_all){
    return(list(crop(ba,bp), 
                fraction_burned,
                ba_comp_bin,
                df_bf_in_comp,
                bp))
  } else{
    return(df_bf_in_comp)
  }
  
}

L <- GetBurnedArea("Kosukhino")

ba <- L[[1]]
fraction_burned <- L[[2]]
ba_comp_bin <- L[[3]]
df_bf_in_comp <- L[[4]]
bp <- L[[5]]

# 3. Plots ----

# figure options
font_size <- 22

# binary_colors <- c("burned" = "#EBB261","unburned" = "#5A4A6F")
magma_2 <- viridis::magma(2,begin = .85,end = .15)
binary_colors <- c("unburned" = magma_2[1],"burned" = magma_2[2])

## Subfigure 1 - burned area map ----
xtext <- 544500
ytext <- 7871400

p1 <- ggplot() +
  geom_spatraster(data = ba, show.legend = FALSE) +
    scale_fill_manual(values = binary_colors,
                      na.value = "white",
                      labels = c('Unburned','Burned')) +
     geom_spatvector(data = bp,
                     fill = NA, linewidth = 1,
                     color = 'black') +
     # scale bar
     geom_rect(aes(xmin = xtext - 2500, xmax = xtext - 500, ymin = ytext, ymax = ytext + 200),
               fill = 'gray20') + 
     geom_text(aes(x = xtext - 1500, y = ytext + 800,
                   label = '2 km',fontface = 'bold'),
               size = 7,
               colour = 'gray20') +
     # north arrow
     geom_text(aes(x = xtext, y = ytext,
                   label = 'N',fontface = 'bold'),
               size = 7, colour = 'gray20') +
     geom_text(aes(x = xtext + 60, y = ytext + 750,
                   label = '$', angle = 90,
                   family = 'ESRI arrowhead',
                   fontface = 'bold'),
               size = 7, colour = 'gray20') +
     # plot styling
     scale_x_continuous(expand=c(0,0)) +
     scale_y_continuous(expand=c(0,0)) +
     theme_map(font_size) +
     theme(plot.subtitle = element_text(hjust = 0.5),
           plot.margin = unit(c(0,0,0,0), "cm"),
           legend.position = "right",
           legend.box = "vertical", 
           legend.justification = "center",
           panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
     labs(fill = NULL,
          subtitle = "Burned area PlanetScope \n (3 m x 3 m)")

## Subfigure 2 - burned fraction map ----
p2 <- ggplot() +
  geom_spatraster(data = fraction_burned, show.legend = FALSE) +
  scale_fill_viridis_c(option = 'rocket',
                       direction = -1,
                       # rescale color ramp to highlight intermediate fractions
                       rescaler = function(x, to = c(0, .7), from = NULL) {
                         ifelse(x<1, 
                                scales::rescale(x,
                                                to = to,
                                                from = c(min(x, na.rm = TRUE), 0.99)),
                                1)},
                         # begin = .15, end = .85,
                         na.value = "white",
                         guide = guide_colorbar(
                           # direction = "horizontal",
                           title.position = 'bottom',
                           title = NULL),
                         limits = c(0,1),
                         labels =label_percent()
    ) +
    geom_spatvector(data = bp,
                    fill = NA, linewidth = 1,
                    color = 'black') +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0)) +
    theme_map(font_size) +
    theme(plot.subtitle = element_text(hjust = 0.5),
          plot.margin = unit(c(0,0,0,0), "cm"),
          legend.position = "right",
          legend.direction = "vertical",
          legend.justification = "center",
          legend.key.width = unit(0.05, "npc"),
          legend.key.height = unit(.152, "npc"),
          panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
    labs(fill = NULL,subtitle = "Aggregated burned fraction \n (30 m x 30 m)")

## Subfigure 3 - Landsat burned area map ----
p3 <- ggplot() +
   geom_spatraster(data = ba_comp_bin, show.legend = FALSE) +
   scale_fill_manual(values = binary_colors,
                     na.value = "white") +
   geom_spatvector(data = bp,
                   fill = NA, linewidth = 1,
                   color = 'black') +
  # annotate perimeter
  annotate("segment", 
           x = xtext - 500, xend = 545500, y = ytext - 900, yend = 7869200, 
           colour = "black",linewidth = 1) +
  # annotate("text", x = 5, y = 2.15, label = "Some text")
  geom_text(aes(x = xtext - 800, y = ytext,
                label = 'fire \nperimeter',fontface = 'bold'),
            size = 7,
            colour = 'black') +
   scale_x_continuous(expand=c(0,0)) +
   scale_y_continuous(expand=c(0,0)) +
   theme_map(font_size) +
   theme(plot.subtitle = element_text(hjust = 0.5),
         plot.margin = unit(c(0,1,0,0), "cm"),
         panel.border = element_rect(colour = "black", fill=NA, linewidth=1)) +
   labs(fill = NULL,
        subtitle = "Burned area Landsat \n (30 m x 30 m)")

## Subfigure 4 - plot distributions of burned fraction per Landsat class ----
df_bf_in_comp <- GetBurnedArea("Kosukhino",
                               # product = "GABAM",
                               return_all = FALSE) %>% 
  group_by(group) %>% 
  sample_frac(size = .3)

p4 <- ggplot(data = df_bf_in_comp,
             aes(x = group,y = burned_fraction,fill = group)) +
  geom_flat_violin(position = position_nudge(x = 0.2, y = 0), 
                   scale = "width",
                   size = 0.2,alpha = 0.8) +
  geom_point(aes(y = burned_fraction, color = group),
             position = position_jitter(width = 0.15), 
             size = 1, alpha = 0.1) +
  # geom_boxplot(lwd = 0.3, width = .2,outlier.shape = NA, alpha = 0.6) +
  labs(y = NULL, x = NULL) +
  scale_y_continuous(labels = label_percent(),
                     expand = c(0,0)) +
  scale_x_discrete(labels = c('Landsat \n Burned',
                              'Landsat \n Unburned') ) +
  scale_fill_manual(values = binary_colors) +
  scale_color_manual(values = binary_colors) +
  theme_cowplot(font_size) + 
  theme(legend.position = 'none',
        plot.margin = unit(c(0,1,0,0), "cm"),
        # plot.background = element_rect(fill = 'white'),
        aspect.ratio=1) 

## make legend for binary maps ----
# create dummy plot
gplot <- ggplot(data.frame(X = c(1,1),
                           Y = c("burned","unburned"),
                           Class = c("burned","unburned")),
                aes(X, Y, color = Class)) +    
  geom_point(size = 7,shape = 15) +
  scale_color_manual(values = binary_colors,
                     labels = c("Burned","Unburned")) +
  theme_cowplot(font_size) + 
  theme(legend.title = element_blank(),
        plot.margin = unit(c(0,1,0,0), "cm"))

# Grab legend from gplot and draw
leg <- ggdraw(get_legend(gplot))    

## make color bar ----
cb <- ggplot() +
  # 0 - 99 % tile
  geom_tile(aes(x = 1, y = seq(0,1,.01),
                fill = seq(0,1,.01)),
            show.legend = FALSE ) +
  # 100 % tile
  geom_tile(aes(x = 1, y = seq(0.995,1.1,.01)),
            fill = "black",
            show.legend = FALSE ) +
  # adjust y axis styling
  scale_y_continuous(
    # limits = c(0,.99),
    labels=label_percent(),
    breaks = c(0,.25,.5,.75,.99),
    position = "left",
    expand = c(0,0)) +
  # apply colorscheme from map
  scale_fill_viridis_c(option = 'rocket',
                       direction = -1,
                       # rescale color ramp to highlight intermediate fractions
                       rescaler = function(x, to = c(0, .7), from = NULL) {
                         ifelse(x<1, 
                                scales::rescale(x,
                                                to = to,
                                                from = c(min(x, na.rm = TRUE), 0.99)),
                                1)},
                       # begin = .15, end = .85,
                       na.value = "white",
                       guide = guide_colorbar(
                         title.position = 'bottom',
                         title = NULL),
                       limits = c(0,1),
                       labels = label_percent() ) +
  # annotate within bar
  annotate("text", x = 1, y = .5,
           label = "Burned fraction", color = "white",
           size = 8,angle = 90, hjust = 0.5, vjust = 0.5) +
  annotate("text", x = 1, y = 1.05,
           label = "100%", color = "white",
           size = 5, hjust = 0.5, vjust = 0.5) +
  theme_cowplot(font_size) +
  theme(plot.margin = unit(c(0,1,0,0), "cm"),
        # adjust x axis
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.line.x=element_blank(),
        # adjust y axis
        axis.title.y=element_blank(),
        axis.ticks.length = unit(0.3, "cm"))

## make plot grid ----
pmap <- ggdraw() + 
   draw_image('figures/Fig_1_v3.png', scale = 0.8) +
  theme(plot.margin = unit(c(0,0,0,0), "cm"))

bottom_four <- cowplot::plot_grid(
                  p1,leg,NULL,NULL,p3, 
                  p2,NULL,cb,NULL,p4,
                  labels = c('b)','','','','c)',
                             'd)','','','','e)'),
                  label_size = font_size, 
                  hjust = 0, label_x = .1,
                  greedy = FALSE,
                  ncol = 5,
                  rel_widths = c(1,-0.1,.4,-.05,1),
                  align = 'hv', axis = 'trbl')

pg <- plot_grid(pmap,NULL, bottom_four, 
                rel_heights = c(1,-.08,1),
                labels = c('a)', ''),
                label_x = 0.05,label_y = .9,
                label_size = font_size,
                ncol = 1)

ggsave2(pg,filename = 'figures/Figure_1.png',
        device = png, type = "cairo",
        bg = 'white',width = 10, height = 18)


# export figures to png
ggsave(p1,filename = 'figures/Fig_1a.png',
       bg = NULL,width = 12, height = 8)
ggsave(p2,filename = 'figures/Fig_1b.png',
       bg = NULL,width = 12, height = 8)
ggsave(p3,filename = 'figures/Fig_1c.png',
       bg = NULL,width = 12, height = 8)
ggsave(p4,filename = 'figures/Fig_1d_GABAM_Kosukhino.png',
       bg = NULL,width = 12, height = 8)

# 4. Statistical tests ----
tabl <- tibble()
# product <- "GABAM"
product <- "descals"

for (aoi_name in aois$site){
  L <- GetBurnedArea(aoi_name,product = product) 
  
  ba <- L[[1]]
  fraction_burned <- L[[2]]
  ba_comp_bin <- L[[3]]
  df_bf <- L[[4]]
  bp <- L[[5]]
  
  # (p4 <- df_bf %>% 
  #     group_by(group) %>% 
  #     sample_frac(size = .3) %>% 
  #     ggplot(aes(x = group,y = burned_fraction,fill = group)) +
  #     geom_flat_violin(position = position_nudge(x = 0.2, y = 0),
  #                      scale = "width",
  #                      size = 0.2,alpha = 0.8) +
  #     geom_point(aes(y = burned_fraction, color = group),
  #                position = position_jitter(width = 0.15), size = 1, alpha = 0.1) +
  #     geom_boxplot(lwd = 0.3, width = .2,outlier.shape = NA, alpha = 0.6) +
  #     labs(y = "Burned fraction \n", x = NULL) +
  #     scale_y_continuous(labels = label_percent(),
  #                        expand = c(0,0)) +
  #     scale_x_discrete(labels = c('Burned', 'Unburned') ) +
  #     scale_fill_manual(values = binary_colors) +
  #     scale_color_manual(values = binary_colors) +
  #     theme_cowplot(font_size) +
  #     theme(legend.position = 'none'))

  # export to png
  # ggsave(p4,filename = sprintf('figures/Fig_1d_%s.png',aoi_name),
  #        bg = NULL,width = 12, height = 8)

  stats <- df_bf %>% 
    group_by(group) %>% 
    dplyr::summarise(md = median(burned_fraction, na.omit = T),
                     ps_burned_area = sum(burned_fraction)* 900/1e6,
                     ps_unburned_area = sum(1 - burned_fraction)* 900/1e6,
                     ls_area = n() * 900 / 1e6)
    
  # Wilcoxon signed-rank test: are burned fractions in groups are different from one?
  wcx_greater <- df_bf %>% 
    group_by(group) %>% 
    sample_frac(size = .3) %>% 
    group_map(~ 
                broom::tidy(wilcox.test(.x$burned_fraction, 
                                 mu = 0.95, alternative = "greater")) %>% 
                mutate(group = .y$group)
    ) %>% 
    bind_rows()
  
  # Wilcoxon signed-rank test: are burned fractions in groups are different from zero?
  wcx_less <- df_bf %>% 
    group_by(group) %>% 
    sample_frac(size = .3) %>% 
    group_map(~ 
                broom::tidy(wilcox.test(.x$burned_fraction, 
                                 mu = 0.05, alternative = "less")) %>% 
                mutate(group = .y$group)
    ) %>% 
    bind_rows()
  
  tabl <- rbind(tabl,
               tibble(site = aoi_name,
                      group = wcx_greater$group,
                      md = stats$md * 100,
                      ls_area = stats$ls_area,
                      ps_unburned_area  = stats$ps_unburned_area,
                      ps_burned_area  = stats$ps_burned_area,
                      wcx_gt = wcx_greater$p.value,
                      wcx_lt = wcx_less$p.value) )
}


## a) make Table S6 & 7 ----
tabl %>%
  select(-c(md,ls_area,ps_unburned_area,ps_burned_area)) %>% 
  mutate(wcx_gt = if_else(wcx_gt < 0.01,"Yes","No"),
         wcx_lt = if_else(wcx_lt < 0.01,"Yes","No")) %>%
  pivot_wider(names_from = group, values_from = c(wcx_gt, wcx_lt)) %>%
  mutate(site = str_replace_all(site, "(?<=[a-z])(?=[A-Z])", " ")) %>%
  gt() %>% 
  tab_spanner(
    label = "Is burned fraction >= 95%?",
    columns  = c(wcx_gt_burned,wcx_gt_unburned),
    id = "gt"
  ) %>% 
  tab_spanner(
    label = "Is burned fraction <= 5%?",
    columns  = c(wcx_lt_burned,wcx_lt_unburned),
    id = "lt"
  ) %>% 
  cols_label(
    site = "Site",
    wcx_gt_burned = "burned",
    wcx_gt_unburned = "unburned",
    wcx_lt_burned = "burned",
    wcx_lt_unburned = "unburned"
  ) %>% 
  tab_style(
    style = cell_text(v_align = "top", weight = 'bold'),
    locations = cells_column_labels()) %>% 
  tab_style(
    style = cell_text(v_align = "top", weight = 'bold'),
    locations = cells_column_spanners()) %>% 
  tab_style(
    style = cell_text(style = "italic"),
    locations = cells_body(rows = wcx_gt_burned =="Yes",
                           columns = wcx_gt_burned)
  ) %>% 
  tab_footnote(
    footnote = "Wilcoxon signed rank test (one-sided, H1: median is less than 0.95).",
    locations = cells_column_spanners(spanners = "gt")
  ) %>% 
  tab_footnote(
    footnote = "Wilcoxon signed rank test (one-sided, H1: median is greater than 0.05).",
    locations = cells_column_spanners(spanners = "lt")
  ) %>% 
  gtsave(filename =  sprintf("tables/Table2_%s.html",product))

## b) make Table S8 & 9 ----
tabl %>%
  select(-c(wcx_gt,wcx_lt)) %>% 
  pivot_wider(names_from = group, values_from = c(md, ls_area,ps_unburned_area,ps_burned_area)) %>%
  mutate(site = str_replace_all(site, "(?<=[a-z])(?=[A-Z])", " ")) %>%
  gt() %>% 
  tab_spanner(
    label = "Landsat burned",
    columns  = c(md_burned,ls_area_burned, ps_unburned_area_burned,ps_burned_area_burned ),
    id = "bd"
  ) %>% 
  tab_spanner(
    label = "Landsat unburned",
    columns  = c(md_unburned,ls_area_unburned, ps_unburned_area_unburned,ps_burned_area_unburned ),
    id = "ubd"
  ) %>% 
  cols_label(
    site = "Site",
    md_burned = "Median burned fraction (%)",
    md_unburned = "Median",
    ls_area_burned = "LS8 burned area (km2)",
    ps_unburned_area_burned = "PS unburned area (km2)",
    ps_burned_area_burned = "PS burned area (km2)",
    ls_area_unburned = "LS8 unburned area (km2)",
    ps_unburned_area_unburned = "PS unburned area (km2)",
    ps_burned_area_unburned = "PS burned area (km2)"
  ) %>% 
  fmt_number(decimals = 1) %>% 
  tab_style(
    style = cell_text(v_align = "top", weight = 'bold'),
    locations = cells_column_labels()) %>% 
  tab_style(
    style = cell_text(v_align = "top", weight = 'bold'),
    locations = cells_column_spanners()) %>% 
  gtsave(filename = sprintf("tables/Table3_%s.html",product))
