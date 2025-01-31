library(landscapemetrics)
library(terra)
library(tidyterra)
library(tidyverse)
library(ggpubr)
library(patchwork)
library(cowplot)
library(colorspace)
library(gt)
set.seed(1234)

# 0. Load data ----
aoi_names_new <- c("Berelech" = "Berelech",
                   "LargeScarCenter" = "Lapcha",
                   "LargeScarControl" = "Keremesit",
                   "Libenchik" = "Sala",
                   "Kosukhino" = "Kosukhino",
                   "DrainedThawlake" = "Ebelyakh"
)
df_all <- data.frame()

for (aoi_name in names(aoi_names_new)){
  aoi_name_new <- aoi_names_new[[aoi_name]]
  
  cat(sprintf("Loading data for %s fire scar \n",aoi_name_new))
  
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
  
  product <- "descals"
  
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
  
  # Load predictors
  r_preds <- rast(paste0("data/geodata/raster/predictors/",aoi_name_new,"_predictors_30m.tif")) %>% 
    crop(bp) %>% 
    mask(bp)
  
  # Calculate area discrepancy (pixel-wise difference between burned area Landsat - Planetscope)
  r_area_landsat <- ifel(ba_comp_bin == "burned", 900, 0) # convert binary to area burned
  r_area_planetscope <- r_preds$burned_fraction * 900 # covnert fraction to area burned
  
  r_diff <- r_area_landsat - r_area_planetscope
  names(r_diff) <- "AreaDifference"
  
  # convert raster stack to dataframe
  df_aoi <- c(ba_comp_bin,r_preds,r_diff) %>% 
    as.data.frame() %>% 
    drop_na() %>% 
    group_by(class) %>%
    sample_frac(0.1) %>% 
    ungroup()
  
  # Calculate clusters on k-means
  df_aoi <- df_aoi %>% 
    mutate(AreaDifferenceClass = cut(
      AreaDifference,
      breaks = c(-Inf, -300, 300, Inf), 
      labels = c("Underestimated", "Moderate", "Overestimated")),
      site = aoi_name_new
    )

  df_all <- bind_rows(df_all,df_aoi)
}

# 1. Create Figure s4: Relationship between burned area underestimation vs. NDVI ----
FONT_SIZE <- 22
annotation_size <- 7

# Define the breakpoints and colors
breaks <- c(-Inf, -300, 300, Inf)
colors <- c("lightblue", "lightgreen", "lightpink")
darker_colors <- c("#2A4365", "#2E8B57", "#8B475D")
labels <- c("Underestimated", "Moderate", "Overestimated")

# only use Kosukhino for this analysis
df_single_scar <- df_all %>% 
  filter(site == "Kosukhino")

# Histogram of burned area discrepancy
(fig_s4a <- ggplot(data = df_single_scar) +
    # Add background colors
    geom_rect(aes(
      xmin = breaks[1], xmax = breaks[2], ymin = -Inf, ymax = Inf),
      fill = colors[1], alpha = 0.2) +
    geom_rect(aes(
      xmin = breaks[2], xmax = breaks[3], ymin = -Inf, ymax = Inf),
      fill = colors[2], alpha = 0.2) +
    geom_rect(aes(
      xmin = breaks[3], xmax = breaks[4], ymin = -Inf, ymax = Inf),
      fill = colors[3], alpha = 0.2) +
    geom_histogram(aes(x = AreaDifference)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0,3100)) +
    scale_x_continuous(breaks = seq(-900,900,300),
                       labels = seq(-900,900,300), 
                       limits = c(-900,900)) +
    labs(x = "Burned area discrepancy \n Landsat - Planetscope (m²)",
         y = "Counts") + 
    # Add labels for background colors
    annotate("text", x = -600, y = 3000, label = "Under-\nestimated", 
             vjust = 1, color = darker_colors[1],
             size = annotation_size, fontface =2) +
    annotate("text", x = 0, y = 3000, label = labels[2], 
             vjust = 1, color = darker_colors[2],
             size = annotation_size, fontface =2) +
    annotate("text", x = 600, y = 3000, label = "Over-\nestimated", 
             vjust = 1, color = darker_colors[3],
             size = annotation_size, fontface =2) +
    coord_cartesian(clip = 'off') + 
    theme(plot.margin = unit(c(3,1,1,1), "lines")) +
    theme_cowplot(FONT_SIZE))

# Plot scatter NDVI vs. sd NDVI
(fig_s4b <- ggplot(data = df_single_scar, aes(x = NDVI, y = NDVI_sd,
                                     color = AreaDifference)) +
    scale_color_continuous_diverging(
      name = "Burned area discrepancy\nLandsat - Planetscope (m²)",
      breaks = seq(-900,900,300),
      labels = seq(-900,900,300)) +
    geom_point(alpha = 0.1) + 
    theme_cowplot(FONT_SIZE) +
    labs(x = expression(Greenness~(NDVI[Landsat])), 
         y = expression(atop(Greenness~heterogeneity, 
                             paste("(", sigma~NDVI[Planet], ")")))) +
    theme(legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title.align = 0.5, 
          legend.justification = c(0.5, 0),
          legend.box.spacing = unit(0.5, "cm"),
          legend.key.width = unit(2, "cm"),
          legend.title = element_text(hjust = 0.5)) +
    guides(color = guide_colorbar(title.position = "top")))

# Boxplot clusters vs. NDVI
kw_ndvi <- kruskal.test(NDVI ~ AreaDifferenceClass, data = df_single_scar)
(fig_s4c <- ggboxplot(df_single_scar, x = "AreaDifferenceClass", y = "NDVI",
                     color = "AreaDifferenceClass", palette = colors,
                     label.select = ) +
    stat_compare_means() +
    rremove("legend") +
    labs(y = expression(Greenness~(NDVI[Landsat])),
         x = "") +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs_pubr(FONT_SIZE))

# Boxplot clusters vs. NDVI_sd
kw_ndvi_sd <- kruskal.test(NDVI_sd ~ AreaDifferenceClass, data = df_single_scar)
(fig_s4d <- ggboxplot(df_single_scar, x = "AreaDifferenceClass", y = "NDVI_sd",
                     color = "AreaDifferenceClass", palette = colors,
                     show.legend = FALSE) +
    stat_compare_means() +
    rremove("legend") +
    labs(y = expression(atop(Greenness~heterogeneity, 
                             paste("(", sigma~NDVI[Planet], ")"))),
         x = "") +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs_pubr(FONT_SIZE))

# build grid and export
fig_s4 <- fig_s4a + fig_s4b + fig_s4c + fig_s4d +
  plot_layout(nrow = 2,heights = c(1,1)) +
  plot_annotation(tag_levels = 'a', tag_suffix = ")") &
  theme(plot.tag = element_text(size = FONT_SIZE, face = "bold"))

ggsave(fig_s4, filename ='figures/Figure_S4.png',
       bg = 'white',width = 14, height = 14)

  # 2. Fit regression line NDVI and Area Discrepancy ----
df_lm <- df_all %>% 
  filter(site == "Kosukhino") %>% 
  filter(AreaDifference <= 0)
m1 <- lm(AreaDifference ~ NDVI, df_lm) %>% summary()
m2 <- lm(AreaDifference ~ NDVI_sd, df_lm) %>% summary()
m3 <- lm(AreaDifference ~ NDVI + NDVI_sd, df_lm) %>% summary()

extract_model_info <- function(summary_model, terms) {
  # get coefficients
  coefs <- as.data.frame(summary_model$coefficients)
  coefs <- coefs[terms, , drop = FALSE] %>%
    select(Estimate = Estimate)
  
  # Add normal and adjusted R-squared
  r2 <- data.frame(
    Term = "R-squared",
    Estimate = summary_model$r.squared * 100
  )
  adj_r2 <- data.frame(
    Term = "Adjusted R-squared",
    Estimate = summary_model$adj.r.squared * 100
  )
  
  # Combine coefficients and adjusted R²
  coefs <- coefs %>%
    rownames_to_column("Term") %>%
    bind_rows(r2) %>% 
    bind_rows(adj_r2)
  
  coefs
}

# Specify terms of interest
terms_of_interest <- c("(Intercept)", "NDVI", "NDVI_sd")

# Extract info for all models
model1_info <- extract_model_info(m1, terms_of_interest) %>%
  rename(Model1 = Estimate)
model2_info <- extract_model_info(m2, terms_of_interest) %>%
  slice(-3) %>% 
  rename(Model2 = Estimate)
model3_info <- extract_model_info(m3, terms_of_interest) %>%
  rename(Model3 = Estimate)

# Combine results into one table
results_table <- model1_info %>%
  full_join(model2_info, by = "Term") %>%
  full_join(model3_info, by = "Term") %>% 
  slice(-3) %>% 
  mutate(Term = factor(Term, 
                       levels = c("(Intercept)", "NDVI", "NDVI_sd",
                                  "R-squared", "Adjusted R-squared"))) %>%
  arrange(Term)

# Format the table with gt
t_s8 <- results_table %>%
  gt() %>%
  sub_missing(
    columns = everything(),
    rows = everything(),
    missing_text = "---"  
  ) %>% 
  fmt_number(
    columns = c(Model1, Model2, Model3),
    decimals = 1
  ) %>%
  cols_label(
    Term = "",
    Model1 = md("D ~ *NDVI*"),
    Model2 = md("D ~ &sigma; *NDVI*"),
    Model3 = md("D ~ &sigma; *NDVI* + *NDVI*")
  ) %>%
  tab_options(
    table.font.size = "small",
    heading.align = "center"
  )
t_s8
# export
# gtsave(t_s8,filename = "tables/TableS8.html")

