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

### i.component 1: ZOI ----
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

### ii. component 2: COI ----
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
   
### iii. component 3: MU ----
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

## c) Plot grid and export ----
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

