# Script to plot model results
# Nils Rietze: nils.rietze@uzh.ch 
# 26 March 2024

library(terra)
library(tidyterra)
library(tidyverse)

library(scico)
library(RColorBrewer)
library(cowplot)
library(scales)

library(brms)
library(tidybayes)
library(DHARMa)
library(sjPlot)
library(marginaleffects)
library(broom.mixed)

# 1. Configure and load stuff ----
aois <- vect('./data/geodata/feature_layers/aoi_wv/aois_analysis.geojson') %>%
  project('EPSG:32655') %>% 
  mutate(., id = 1:nrow(.))

# zoib_model <- readRDS("zoi_beta_model_2024-03-20.rds") # zoi, phi & coi are identities
mod <- readRDS("zoi_beta_model_2024-03-28.rds") # zoi, phi & coi are fitted with predictors
summary(mod)
model_data_scaled <- mod$data

# 2. Model quality reporting ----
# report ESS
bayestestR::effective_sample(mod)

# Generate pp_check densities
sim <- pp_check(mod)
print(sim) # plot pp_check

# Proportions of 0/1 in data
prop_ones <- length(which(model_data_scaled$burned_fraction == 1))/nrow(model_data_scaled)
prop_zeros <- length(which(model_data_scaled$burned_fraction == 0))/nrow(model_data_scaled)
prop_zeros_or_ones <- prop_ones + prop_zeros

# Proportions of 0/1 in simulated data
sim_data <- sim$data
ones_sim <- length(which(sim_data$value == 1))
zeros_sim <- length(which(sim_data$value == 0))
zeros_or_ones_sim <- ones_sim + zeros_sim

# report proportion of obs. & fitted extremes
data.frame(obs = c(prop_zeros,prop_ones),
           fit = c(zeros_sim/nrow(sim_data),ones_sim/nrow(sim_data)),
           row.names = c("0s","1s")) %>% 
  print()

# Binomial Tests
# Zeros & Ones
b_ZO <- binom.test(zeros_or_ones_sim, nrow(sim_data), p = prop_zeros_or_ones,
                   alternative = c("two.sided"),
                   conf.level = 0.95)
# If the p-value is high (typically greater than 0.05), it suggests that the observed 
# probability of zeros and ones is not significantly different from the expected probability, 
# that your model's prediction is accurate
b_ZO$p.value

# Zeros
b_Z <- binom.test(zeros_sim, nrow(sim_data), p = prop_zeros,
                  alternative = c("two.sided"),
                  conf.level = 0.95)
# If the p-value is high (typically greater than 0.05), it suggests that the observed 
# probability of zeros is not significantly different from the expected probability, 
# that your model's prediction is accurate
b_Z$p.value

# Ones
b_O <- binom.test(ones_sim, nrow(sim_data), p = prop_ones,
                  alternative = c("two.sided"),
                  conf.level = 0.95)
# If the p-value is high (typically greater than 0.05), it suggests that the observed 
# probability of ones is not significantly different from the expected probability, 
# that your model's prediction is accurate
b_O$p.value

# 3. Make plots ----

## a) Plot histogram icons ----

# create new column indicating the 0s or 1s
histogram_data <- model_data_scaled %>% 
  filter(site == "Kosukhino") %>% 
  mutate(is_zero_or_one = (burned_fraction == 0) | (burned_fraction == 1)) %>% 
  mutate(is_zero = burned_fraction == 0,
         burned_fraction = ifelse(is_zero, burned_fraction - 0.01, burned_fraction)) %>% 
  mutate(is_one = burned_fraction == 1,
         burned_fraction = ifelse(is_one & !is_zero, burned_fraction + 0.01, burned_fraction)) 

# histogram for ZOI component
hist_zoi <- ggplot(histogram_data, aes(x = burned_fraction )) +
  geom_histogram(aes(fill = is_zero_or_one), binwidth = 0.1, 
                 boundary = 0,show.legend = FALSE) +
  # geom_vline(xintercept = 0,linewidth = 1) +
  # geom_vline(xintercept = 1,linewidth = 1) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_fill_manual(values = c("gray","gray12")) +
  theme_map()

# histogram for COI component
hist_coi <- ggplot(histogram_data, aes(x = burned_fraction )) +
  geom_histogram(aes(fill = is_one), binwidth = 0.1, 
                 boundary = 0,show.legend = FALSE) +
  # geom_vline(xintercept = 0,linewidth = 1) +
  # geom_vline(xintercept = 1,linewidth = 1) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_fill_manual(values = c("gray","gray12")) +
  theme_map()

# histogram for intermediate fractions
hist_mu <- ggplot(histogram_data, aes(x = burned_fraction )) +
  geom_histogram(aes(fill = !is_zero_or_one), binwidth = 0.1, 
                 boundary = 0,show.legend = FALSE) +
  # geom_vline(xintercept = 0,linewidth = 1) +
  # geom_vline(xintercept = 1,linewidth = 1) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_fill_manual(values = c("gray","gray12")) +
  theme_map()

## b) Plot estimated effects ----
# define labels for plot
mod_labs <-c(b_elevation = 'Elevation',
             b_slope = 'Slope',
             b_northness = 'Northness',
             b_eastness = 'Eastness',
             b_tpi_500 = expression(Landform~(TPI[500~m])),
             b_LST = expression(Land~surface~cooling~(LST[Landsat])),
             b_NDVI = expression(Greenness~(NDVI[Landsat])), 
             b_NDVI_sd = expression(Greenness~heterogeneity~(sigma~NDVI[Planet]))
)

# Plot the posterior estimates, more info here: 
# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html
font_size <- 22

### i. Plot posterior estimates for zero-one inflation model ----
p1 <- mod %>%
   gather_draws(`b_.*`, regex = TRUE) %>%          # filter estimated effects
   filter(!grepl("Intercept", `.variable`)) %>%    # remove Intercept from list
   filter(grepl("zoi", `.variable`)) %>%           # remove precision parameters
   mutate(`.variable` = gsub("zoi_", "", `.variable`) ) %>%
   mutate(`.variable` = factor(`.variable`,levels = names(mod_labs))) %>% 
   ggplot(aes(y = .variable, x = .value)) +
    stat_halfeye(.width = c(0.05,0.95),size = 0.5,fill = "#E8CEB6") +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_y_discrete(labels = mod_labs) +
    labs(x = "Coefficient estimate \n(scaled)",
         y = "",
         subtitle = "<span style='color: #E8CEB6'>**a) Comp. 1: extreme or <br>
         intermediate fractions?**</span>") +
   theme_minimal_hgrid(font_size) +
    theme(
      plot.subtitle = ggtext::element_markdown(size = font_size),
      axis.title.x = element_text(hjust = 0.5),
      legend.position = "none") 

### ii. Plot posterior estimates for conditional one inflation model ----
p2 <- mod %>%
   gather_draws(`b_.*`, regex = TRUE) %>%          # filter estimated effects
   filter(!grepl("Intercept", `.variable`)) %>%    # remove Intercept from list
   filter(grepl("coi", `.variable`)) %>%           # remove precision parameters
   mutate(`.variable` = gsub("coi_", "", `.variable`) ) %>%
   mutate(`.variable` = factor(`.variable`,levels = names(mod_labs))) %>% 
   ggplot(aes(y = .variable, x = .value)) +
   stat_halfeye(.width = c(0.05,0.95),size = 0.5,fill = "#BF96AB") +
   geom_vline(xintercept = 0, linewidth = 0.3) +
   scale_y_discrete(labels = mod_labs) +
   labs(x = "Coefficient estimate \n(scaled)",
        y = "",
        subtitle = "<span style='color: #BF96AB'>**b) Comp. 2: effect on <br>
        full burned fraction**</span>") +
   theme_minimal_hgrid(font_size) +
   theme(
     axis.title.x = element_text(hjust = 0.5),
     plot.subtitle = ggtext::element_markdown(size = font_size),
     legend.position = "none") 
   
### iii. Plot posterior estimates for mean model ----
p3 <- mod %>%
   gather_draws(`b_.*`, regex = TRUE) %>%          # filter estimated effects
   filter(!grepl("Intercept", `.variable`)) %>%    # remove Intercept from list
   filter(!grepl("phi", `.variable`)) %>%          # remove precision parameters
   filter(!grepl("zoi", `.variable`)) %>%          # remove zero-one-inflation
   filter(!grepl("coi", `.variable`)) %>%          # remove conditional one inflation
   mutate(`.variable` = factor(`.variable`,levels = names(mod_labs))) %>% 
   ggplot(aes(y = .variable, x = .value)) +
   stat_halfeye(.width = c(0.05,0.95),size = 0.5,fill = "#b5ccb9") +
   geom_vline(xintercept = 0, linewidth = 0.3) +
   scale_y_discrete(labels = mod_labs) +
   labs(x = "Coefficient estimate \n(scaled)", 
        y = "",
        subtitle = "<span style='color: #b5ccb9'>**c) Comp. 3: effect on <br>
        intermediate fractions**</span>") +
   theme_minimal_hgrid(font_size) +
   theme(
     axis.line.y = element_line(colour = 'gray'),
     plot.subtitle = ggtext::element_markdown(size = font_size),
     axis.title.x = element_text(hjust = 0.5),
     legend.position = "none") 

### iv. Plot grid and export ----
pg <- cowplot::plot_grid(ggdraw(p1) +
                     draw_plot(hist_zoi, .9, .85, .12, .1) , 
                   ggdraw(p2 + 
                            theme(axis.text.y = element_blank(),
                                  # axis.line.y = element_blank(),
                                  axis.line.y = element_line(colour = 'gray'),
                                  axis.title.y= element_blank(),
                                  axis.ticks.y= element_blank())) +
                     draw_plot(hist_coi, .7, .85, .25, .1) ,
                   ggdraw(p3 + 
                            theme(axis.text.y = element_blank(),
                                  # axis.line.y = element_blank(),
                                  axis.line.y = element_line(colour = 'gray'),
                                  axis.title.y= element_blank(),
                                  axis.ticks.y= element_blank())) +
                     draw_plot(hist_mu, .7, .85, .25, .1),
                   nrow = 1,
                   rel_widths = c(1.6, .8, .8),
                   align = 'h', axis = 'tb')

ggsave2(pg, filename = 'figures/Figure_3.png',
        bg = 'white',width = 18, height = 8)

ggsave(p1,filename = sprintf('figures/model/zoib_mean_model_%s.png',today()),
       bg = 'white',width = 12, height = 8)
ggsave(p3,filename = sprintf('figures/model/zoib_coi_model_%s.png',today()),
       bg = 'white',width = 12, height = 8)
ggsave(p2,filename = sprintf('figures/model/zoib_zoi_model_%s.png',today()),
       bg = 'white',width = 12, height = 8)

### v. Plot random effects for mean model ----
(p4 <- mod %>%
  gather_draws(`r_site\\[.*`, regex = TRUE) %>%          # filter estimated effects
  mutate(`.variable` = factor(`.variable`)) %>% 
  ggplot(aes(y = .variable, x = .value)) +
  stat_halfeye(.width = c(0.05,0.95),size = 0.5,fill = '#b5ccb9') +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_y_discrete(labels = mod_labs) +
  theme_minimal_hgrid() + 
  theme(panel.grid = element_blank(),
        legend.position = "none") +
  labs(x = "Posterior estimates (standardized)", y = ""))

ggsave(p4,filename = sprintf('figures/model/RE_zoib_model_%s.png',today()),bg = 'white')


## d) Visualize conditional effects----
ce <- conditional_effects(mod,dpar = 'mu')
# plot(ce,points = F,rug = T)

labels = c(elevation = 'Elevation',
           slope = 'Slope',
           aspect = 'Aspect',
           northness = 'Northness',
           eastness = 'Eastness',
           tpi_500 = expression(Landform~(TPI[500~m])),
           LST = expression(Land~surface~cooling~(LST[Landsat])),
           NDVI = expression(Greenness~(NDVI[Landsat])), 
           NDVI_sd = expression(Greenness~heterogeneity~(sigma~NDVI[Planet])),
           burned_fraction = NULL)

# Function to plot conditional effects
plot_conditional_effects <- function(df, model_data_scaled, labels, save_plots = FALSE){
  
  # get variable name
  var_name <- colnames(df)[1]
  xlabel <- labels[var_name]
  
  # extract scaling factors
  m <-attr(model_data_scaled[[var_name]],'scaled:center')
  ss <-attr(model_data_scaled[[var_name]],'scaled:scale')
  
  p <- ggplot(df) + 
    geom_ribbon(aes(x = effect1__* ss + m,ymin = lower__,ymax = upper__),
                fill = 'grey70',alpha = .7) +
    geom_line(aes(x = effect1__* ss + m ,y = estimate__)) +
    geom_rug(data = attr(df,which = 'points'),
             aes(x = !!sym(var_name)* ss + m),
             sides = 'b', inherit.aes = F) +
    scale_y_continuous(labels = label_percent()) +
    scale_x_continuous(labels = label_number_auto()) +
    labs(x = xlabel, y = 'Burned fraction') +
    theme_cowplot()
  
  print(p)
  
  if(save_plots){
    ggsave(sprintf('figures/model/ce_%s.png',var_name),
           width = 10, height = 10, bg = 'white')
  }
}

lapply(ce, plot_conditional_effects,labels = labels, 
       model_data_scaled = model_data_scaled,save_plots = FALSE)


## f) plot marginal effects ----

# create new data for predictions
new_tibble <- tibble(row_id = 1:100)

for (col_name in names(model_data_scaled)[1:9]) {
  # get range of values in the original dataframe
  range_values <- range(model_data_scaled[[col_name]], na.rm = TRUE)
  # make sequence of 100 points within the range
  seq_values <- seq(from = range_values[1], to = range_values[2], length.out = 100)
  # add  sequence of values to the new tibble
  new_tibble[[col_name]] <- seq_values
}

# expand tibbles with all 6 AOI names
new_tibble <- expand_grid(new_tibble,site = aoi_names)

### i. plot average of the posterior predictive distribution of the model ----
# Mean process
mod %>% 
  epred_draws(newdata = new_tibble,dpar = 'mu') %>% 
  ggplot(aes(y = .epred, x = NDVI_sd)) +
  stat_lineribbon() +
  scale_y_continuous(labels = label_percent()) +
  scale_fill_viridis_d(option = 'magma',begin = .6) +
  labs(y = "Predicted burned fraction", x = 'NDVI_sd',
       fill = "Credible interval") +
  theme_cowplot() +
  theme(legend.position = "bottom")

# ZOI process
mod %>% 
  epred_draws(newdata = new_tibble, dpar = 'zoi') %>% 
  ggplot(aes(y = zoi, x = NDVI_sd)) +
  stat_lineribbon() +
  scale_fill_viridis_d(option = 'magma',begin = .6) +
  scale_y_continuous(labels = label_percent()) +
  labs(y = "Predicted burned fraction", x = "NDVI_sd",
       fill = "Credible interval",
       title = "ZOI process") +
  theme_cowplot() +
  theme(legend.position = "bottom")

# COI process
mod %>% 
  epred_draws(newdata = new_tibble,dpar = 'coi') %>% 
  ggplot(aes(y = coi, x = NDVI_sd)) +
  stat_lineribbon() +
  scale_fill_viridis_d(option = 'magma',begin = .6) +
  scale_y_continuous(labels = label_percent()) +
  labs(y = "Predicted burned fraction", x = "NDVI_sd",
       fill = "Credible interval",
       title = "COI process") +
  theme_cowplot() +
  theme(legend.position = "bottom")

### ii. predict in successions ----
p_zoi <- mod %>% 
  posterior_predict(dpar = 'zoi') 

p_zoi_mean <- colMeans(p_zoi)


### iii. plot conditional effects ----

# create table where only NDVI_sd has values
range_values <- range(model_data_scaled$NDVI_sd, na.rm = TRUE)
# make sequence of 100 points within the range
seq_values <- seq(from = range_values[1], to = range_values[2], length.out = 100)

fe_only_ndvisd <- tibble(NDVI_sd = modelr::seq_range(model_data_scaled$NDVI_sd, n = 100),
                       elevation = 0,
                       slope = 0,
                       aspect = 0,
                       eastness = 0,
                       northness = 0,
                       tpi_500 = 0,
                       LST = 0, 
                       NDVI = 0) %>% 
  expand_grid(site = aoi_names) %>% 
  add_epred_draws(mod,
                  dpar = 'zoi',
                  ndraws = 100)

ggplot(data = fe_only_ndvisd,aes(x = NDVI_sd, y = .epred,
                                 fill = site,col = site)) +
  stat_lineribbon(alpha = .5) +
  facet_grid(~ site) +
  scale_color_viridis_d(option = 'magma') +
  scale_fill_viridis_d(option = 'magma') +
  labs(x = "NDVI_sd", y = "burned_fraction") +
  theme_cowplot() 
  

### iv. plot predicted burned fractions per site ----
zoib_pred <- mod %>% 
  predicted_draws(newdata = new_tibble) %>% 
  mutate(is_zero_or_one = (.prediction == 0) | (.prediction == 1)) %>% 
  mutate(is_zero = .prediction == 0,
         .prediction = ifelse(is_zero, .prediction - 0.01, .prediction)) %>% 
  mutate(is_one = .prediction == 1,
         .prediction = ifelse((is_one) & (!is_zero), .prediction + 0.01, .prediction)) 

ggplot(zoib_pred, aes(x = .prediction)) +
  geom_histogram(aes(fill = is_zero_or_one), binwidth = 0.025, 
                 boundary = 0,
                 color = "white") +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 1) +
  scale_x_continuous(labels = label_percent()) +
  scale_fill_viridis_d(option = "plasma", end = 0.5,
                       guide = guide_legend(reverse = TRUE)) +
  labs(x = "Predicted burned fractions", 
       y = "Count", fill = "Is zero or one?") +
  facet_wrap(~site, ncol = 2) +
  theme_cowplot() +
  theme(legend.position = "bottom")

# Plot the prediction intervals
mod %>% 
  predicted_draws(newdata = new_tibble) %>% 
  ggplot(aes(y = .prediction, x = LST)) +
    stat_lineribbon() +
    scale_y_continuous(labels = label_percent()) +
    scale_fill_brewer(palette = "Purples") +
    labs(y = "Predicted burned fraction", x = 'LST',
         fill = "Credible interval") +
    theme_cowplot() +
    theme(legend.position = "bottom")

# Plot marginal effects per site
predictions(mod,
            re_formula = NULL,
            newdata = new_tibble) %>% 
  posterior_draws() %>% 
  ggplot(aes(LST, draw, fill = site, color = site)) +
    stat_lineribbon(alpha = .25) +
    facet_grid(~ site) +
    scale_color_viridis_d(option = 'magma') +
    scale_fill_viridis_d(option = 'magma') +
    theme_cowplot()


# Plots ----
# ..........

## Plot all predictors ----
aoi_name <- aoi_names[[3]]

### Pair plot of predictors ----
gp <- ggpairs(model_data[c('burned_fraction','site',mod_vars)])
ggsave(gp,filename = 'figures/model/predictors_pairplot_all.png',
       height = 10, width = 10)

### small map gallery ----
predictors_unscaled <- rast(paste0("data/geodata/raster/predictors/",
                                   aoi_name,"_predictors_",
                                   window_side_length,"m.tif") )

# Define plotting parameters
parameters <- tibble(palettes = list(elevation = 'grayC',
                                     slope = "buda",
                                     aspect = "romaO",
                                     northness = 'vik',
                                     eastness = 'vik',
                                     tpi_500 = "bam",
                                     LST = "lajolla",
                                     NDVI = "cork",
                                     NDVI_sd = 'acton',
                                     burned_fraction = "bilbao"),
                     labels = list(elevation = 'Elevation (m)',
                                   slope = 'Slope (°)',
                                   aspect = 'Aspect (°)',
                                   northness = 'Northness (unitless)',
                                   eastness = 'Eastness (unitless)',
                                   tpi_500 = expression(TPI[500~m]),
                                   LST = 'LST (° C)',
                                   NDVI = expression(paste(NDVI[Landsat],
                                                           " (unitless)")), 
                                   NDVI_sd = expression(paste(sigma~NDVI[Planet],
                                                              " (unitless)")),
                                   burned_fraction = NULL),
                     titles = list(elevation = 'Arctic DEM Elevation',
                                   slope = 'Slope',
                                   aspect = 'Aspect',
                                   northness = 'Northness',
                                   eastness = 'Eastness',
                                   tpi_500 = 'Topographic potitison index',
                                   LST = 'Land surface temperature',
                                   NDVI = 'Landsat NDVI',
                                   NDVI_sd = 'Planet NDVI standard deviation \n(proxy for vegetation heterogeneity)',
                                   burned_fraction = 'Burned fraction')
)

create_ggplot <- function(layer, parameters) {
  lyr_name <- names(layer)
  
  direction <- ifelse(lyr_name == 'burned_fraction',-1,1)
  axis_labels <- ifelse(lyr_name == 'burned_fraction','10 km',"")
  
  # adjust scale limits dynamicall
  
  # see if layer has negative values, if yes apply new scale
  has_negative <- global(layer,'min',na.rm = T) < 0
  
  abs_max <- max(abs(global(layer,'min',na.rm = T)),
                 abs(global(layer,'max',na.rm = T)) )
  
  # Set color scale limits dynamically    
  if (has_negative){
    scale_limits <- c(-abs_max, abs_max)
  } else {
    scale_limits <- c(as.numeric(global(layer,'min',na.rm = T)), abs_max)
  }
  
  ggplot() +
    geom_spatraster(data = layer) +
    scale_fill_scico(palette = parameters$palettes[[lyr_name]], 
                     direction = direction,
                     limits = scale_limits) +
    theme_map(12) +
    theme(legend.position = "bottom",
          legend.box = "horizontal", 
          legend.justification = "center",
          legend.key.width = unit(0.05, "npc"),
          plot.margin = margin(t = 2, r = 0, b = 2, l = 7)) +
    guides(fill = guide_colourbar(title.position="bottom", title.hjust = 0.5)) +
    labs(title = parameters$titles[[lyr_name]],
         fill = parameters$labels[[lyr_name]])
}

ggplots <- lapply(mod_vars[-1], function(layer) {
  create_ggplot(predictors_unscaled[[layer]], parameters)
})

pred_plots <- cowplot::plot_grid(plotlist = ggplots, ncol = 3, 
                                 rel_heights = rep(.9, length(ggplots)),
                                 labels = 'auto',
                                 align = 'vh')

ggsave2(pred_plots,
        filename = sprintf('figures/model/model_predictors_%s_unscaled.png',aoi_name),
        height = 18, width = 18)

## Histograms of burned fractions in each site ----
load_burned_fraction_rasters <- function(aoi_name,window_side_length){
  bp <- vect(sprintf('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_%s.shp',aoi_name))
  df <- rast(paste0("data/geodata/raster/predictors/",aoi_name,
                    "_predictors_",window_side_length,"m.tif") ) %>% 
    select(burned_fraction) %>%
    mask(bp) %>% 
    as.data.frame() %>% 
    mutate(site = aoi_name)
  return(df)
}

histogram_data <- do.call(rbind,
                          lapply(aoi_names,FUN = load_burned_fraction_rasters,window_side_length))

# create new column indicating the 0s or 1s
histogram_data <- histogram_data %>% 
  mutate(is_zero_or_one = (burned_fraction == 0) | (burned_fraction == 1)) %>% 
  mutate(is_zero = burned_fraction == 0,
         burned_fraction = ifelse(is_zero, burned_fraction - 0.01, burned_fraction)) %>% 
  mutate(is_one = burned_fraction == 1,
         burned_fraction = ifelse(is_one & !is_zero, burned_fraction + 0.01, burned_fraction)) 

ggplot(histogram_data, aes(x = burned_fraction )) +
  geom_histogram(aes(fill = is_zero_or_one), binwidth = 0.05, 
                 boundary = 0) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 1) +
  scale_x_continuous(labels = scales::label_percent()) +
  scale_fill_viridis_d(option = "inferno", end = 0.7,
                       guide = guide_legend(reverse = TRUE)) +
  labs(x = "Burned fraction",
       y = "Count", fill = "0 or 100 % burned?") +
  facet_wrap(~site, ncol = 2) + 
  theme_cowplot() +
  theme(legend.position = "bottom")

ggsave('figures/burned_fraction_histograms_in_burnperimeter.png',bg = 'white',
       width = 10,height = 10)

## Scatterplot of observed vs. predicted burn fraction ----
data_df <- data.frame(obs = values(fraction_burned,bp),
                      pred = values(pred))

colnames(data_df) <- c('obs','pred')

ggplot(data_df, aes(x = obs, y = pred)) +
  geom_point() +
  xlab('Burned fraction (Planet)')+
  xlab('Burned fraction (RF predicted)')

## Map of differences observed - predicted burn fractions ----
ggplot() +
  geom_spatraster(data=(fraction_burned - pred) ) +
  scale_fill_scico(palette = 'vik',
                   direction = -1,
                   guide = guide_colourbar(direction = "horizontal",
                                           title.position = 'top')) +
  labs(fill = 'Burned fraction residuals Observed - Predicted') +
  theme_minimal() +
  theme(legend.position = "bottom") 

## Compare histograms of LST products ----
clrs <- c("#5A4A6F",  "#EBB261")

pl <- ggplot() +
  geom_density(data = lst_og, aes(fill = clrs[1], x = LST_og),
               color = clrs[1],alpha = .4) +
  geom_density(data = lst, aes(fill = clrs[2],x = LST),
               color = clrs[2],alpha = .4) +
  scale_fill_manual(values = clrs,name = NULL, 
                    labels = c("Level 2 product (30 m x 30 m)", "Resampled to 100 m x 100 m")) +
  labs(title = "Pre-fire Landsat LST",
       x = "Surface temperature (°C)",
       y = "Density") +
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.box = "horizontal", legend.justification = "center") 

ggsave(pl,filename = 'figures/aoi_LST_densities.png',bg = 'white')




## Plot water area products ----
blue_colors <- brewer.pal(n = 3, name = "Blues")[-2]

# Load water mask in native resolution
wa_og <- rast(wa_path) 

# Set categories of water mask
cls <- data.frame(id=1:2, cover=c('no_water','water'))
levels(wa) <- cls

p1 <- ggplot() +
  geom_spatraster(data = wa_og) +
  scale_fill_manual(values = blue_colors) +
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.box = "horizontal", legend.justification = "center") +
  labs(fill = NULL,
       title = "Water area (3 m x 3 m)")

p2 <- ggplot() +
  geom_spatraster(data = wa) +
  scale_fill_manual(values = blue_colors) +  
  theme_cowplot() +
  theme(legend.position = "bottom",
        legend.box = "horizontal", legend.justification = "center")  +
  labs(fill = NULL,
       title = "Water area (100 m x 100 m)")

combined_plot <- cowplot::plot_grid(p1, p2, ncol = 2,align = 'h');print(combined_plot)

ggsave2(combined_plot,filename = 'figures/water_mask_aggregation_aoi.png',
        bg = 'white',
        width = 12, height = 8)

## Plot band histograms of AOI 5 raster ----

# Define n of pixels to sample
n_sample <- 1e5

# Load and extract raster values
df_2019 <- rast('C:/data/8_planet/2019/cropped/20190901_LargeScarCenter_PS2_composite.tif') %>% 
  spatSample(n_sample,as.df = T, xy = T)
df_2020 <- rast('C:/data/8_planet/2020/cropped/20200911_LargeScarCenter_PS2-SD_composite.tif') %>% 
  terra::extract(df_2019[,c('x','y')],ID = F)
df_2019_unharmonized <- rast('C:/data/8_planet/2019/unharmonized/20190901/composite.tif') %>% 
  rename("Blue"=1,"Green"=2,"Red"=3,"NIR"=4) %>% 
  terra::extract(df_2019[,c('x','y')],ID = F)

# reshape dataframe
planet_dataframe <- data.frame(rbind(select(df_2019,-c(x,y)),df_2020,df_2019_unharmonized) ) %>% 
  mutate( Year = rep(c("2019_harmonized", "2020","2019_unharmonized"), each = n_sample) ) %>% 
  pivot_longer(cols = 1:4)

# Plot histograms
ggplot(data = planet_dataframe) + 
  geom_histogram(bins = 90, aes(x = value,fill = Year), alpha = 0.5) +
  facet_wrap(~name,scales = 'free') +
  labs(x = "Pixel Value", y = "Frequency", 
       title = "Comparison of Bands' Histograms for 2019 and 2020") +
  theme_cowplot()

