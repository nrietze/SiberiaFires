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

# Load ZOIB model results
# mod <- readRDS("zoi_beta_model_2024-03-28.rds") 
mod <- readRDS("zoi_beta_model_2024-08-26.rds")
summary(mod)

# retrieve predictor and response data as table
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
             b_LST = expression(Land~surface~temperature~(LST[Landsat])),
             b_NDVI = expression(Greenness~(NDVI[Landsat])), 
             b_NDVI_sd = expression(Greenness~heterogeneity~(sigma~NDVI[Planet]))
)

# as HTML:
mod_labs <-c(b_elevation = 'Elevation',
             b_slope = 'Slope',
             b_northness = 'Northness',
             b_eastness = 'Eastness',
             b_tpi_500 = "Landform (<i>TPI<sub>500 m</sub></i>)",
             b_LST = "Land surface temperature<br>(<i>LST<sub>Landsat</sub></i>)",
             b_NDVI = "Greenness<br>(<i>NDVI<sub>Landsat</sub></i>)", 
             b_NDVI_sd = "Greenness heterogeneity<br>(<i>&sigma; NDVI<sub>Planet</sub></i>)")
# Plot the posterior estimates, more info here: 
# https://cran.r-project.org/web/packages/tidybayes/vignettes/tidy-brms.html
font_size <- 28
subtitle_vmargin <- 100

### i.component 1: ZOI ----
p1 <- mod %>%
   gather_draws(`b_.*`, regex = TRUE) %>%          # filter estimated effects
   filter(!grepl("Intercept", `.variable`)) %>%    # remove Intercept from list
   filter(grepl("zoi", `.variable`)) %>%           # remove precision parameters
   mutate(`.variable` = gsub("zoi_", "", `.variable`) ) %>%
   mutate(`.variable` = factor(`.variable`,levels = names(mod_labs))) %>% 
   # build plot 
  ggplot(aes(y = .variable, x = .value)) +
    stat_halfeye(.width = c(0.95),size = 4,fill = "#E8CEB6") +
    geom_vline(xintercept = 0, linewidth = 0.3) +
    scale_y_discrete(labels = mod_labs) +
    labs(x = "Coefficient estimate \n(scaled)",
         y = "",
         subtitle = "<span style='color: #E8CEB6'>**a) Component 1: <br> 
         extreme or intermediate <br>
         fractions?**</span>") +
   theme_minimal_hgrid(font_size) +
   theme(
     axis.text.y.left = ggtext::element_markdown(hjust = 1),
     plot.subtitle = ggtext::element_markdown(size = font_size,
                                              margin = ggplot2::margin(0,0,subtitle_vmargin,0)),
     axis.title.x = element_text(hjust = 0.5),
     legend.position = "none");p1

### ii. component 2: COI ----
p2 <- mod %>%
   gather_draws(`b_.*`, regex = TRUE) %>%          # filter estimated effects
   filter(!grepl("Intercept", `.variable`)) %>%    # remove Intercept from list
   filter(grepl("coi", `.variable`)) %>%           # remove precision parameters
   mutate(`.variable` = gsub("coi_", "", `.variable`) ) %>%
   mutate(`.variable` = factor(`.variable`,levels = names(mod_labs))) %>% 
   # build plot 
  ggplot(aes(y = .variable, x = .value)) +
   stat_halfeye(.width = c(0.95),size = 4,fill = "#BF96AB") +
   geom_vline(xintercept = 0, linewidth = 0.3) +
   scale_y_discrete(labels = mod_labs) +
   labs(x = "Coefficient estimate \n(scaled)",
        y = "",
        subtitle = "<span style='color: #BF96AB'>**b) Component 2: <br>
        effect on full <br>
        burned fraction**</span>") +
   theme_minimal_hgrid(font_size) +
   theme(
     axis.title.x = element_text(hjust = 0.5),
     plot.subtitle = ggtext::element_markdown(size = font_size,
                                              margin = ggplot2::margin(0,0,subtitle_vmargin,0)),
     legend.position = "none") 
   
### iii. component 3: MU ----
p3 <- mod %>%
   gather_draws(`b_.*`, regex = TRUE) %>%          # filter estimated effects
   filter(!grepl("Intercept", `.variable`)) %>%    # remove Intercept from list
   filter(!grepl("phi", `.variable`)) %>%          # remove precision parameters
   filter(!grepl("zoi", `.variable`)) %>%          # remove zero-one-inflation
   filter(!grepl("coi", `.variable`)) %>%          # remove conditional one inflation
   mutate(`.variable` = factor(`.variable`,levels = names(mod_labs))) %>% 
   # build plot 
  ggplot(aes(y = .variable, x = .value)) +
   stat_halfeye(.width = c(0.95),size = 4,fill = "#b5ccb9") +
   geom_vline(xintercept = 0, linewidth = 0.3) +
   scale_y_discrete(labels = mod_labs) +
   labs(x = "Coefficient estimate \n(scaled)", 
        y = "",
        subtitle = "<span style='color: #b5ccb9'>**c) Component 3: <br> 
        effect on intermediate <br>
        fractions**</span>") +
   theme_minimal_hgrid(font_size) +
   theme(
     axis.line.y = element_line(colour = 'gray'),
     plot.subtitle = ggtext::element_markdown(size = font_size,
                                              margin = ggplot2::margin(0,0,subtitle_vmargin,0)),
     axis.title.x = element_text(hjust = 0.5),
     legend.position = "none") 

## c) Plot grid and export ----
y_hist <- .81
x_hist <- .6
width_hist <- .24
height_hist <- .08

pg <- cowplot::plot_grid(ggdraw(p1) +
                     draw_plot(hist_zoi, 
                               x_hist, y_hist, width_hist, height_hist
                               ) , 
                   ggdraw(p2 + 
                            theme(axis.text.y = element_blank(),
                                  axis.line.y = element_line(colour = 'gray'),
                                  axis.title.y= element_blank(),
                                  axis.ticks.y= element_blank())) +
                     draw_plot(hist_coi, 
                               x_hist - .3, y_hist,width_hist *2, height_hist
                               ) ,
                   ggdraw(p3 + 
                            theme(axis.text.y = element_blank(),
                                  axis.line.y = element_line(colour = 'gray'),
                                  axis.title.y= element_blank(),
                                  axis.ticks.y= element_blank())) +
                     draw_plot(hist_mu, 
                               x_hist - .3, y_hist,width_hist *2, height_hist
                               ),
                   nrow = 1,
                   rel_widths = c(1.6, .8, .8),
                   align = 'h', axis = 'tb')

# export entire plot
ggsave2(pg, filename = 'figures/Figure_3.png',
        bg = 'white',width = 18, height = 12)

# Save subplots individually
# ggsave(ggdraw(p1) +
#          draw_plot(hist_zoi, 
#                    x_hist, y_hist, width_hist, height_hist
#          ),filename = 'figures/Figure_3_Component1.png',
#        bg = 'white',width = 14, height = 10)
# ggsave(p3,filename = sprintf('figures/model/zoib_coi_model_%s.png',today()),
#        bg = 'white',width = 12, height = 8)
# ggsave(p2,filename = sprintf('figures/model/zoib_zoi_model_%s.png',today()),
#        bg = 'white',width = 12, height = 8)

# 4. Create table of all effect sizes & CIs ----
library(gt)
remove <- c("coi_","zoi_")

mod %>% 
  summarise_draws( mean = ~mean(.x),
                   q0.025 = ~quantile(.x, 0.025),
                   q0.975 = ~quantile(.x, 0.975),
                   .num_args = list(sigfig = 2, notation = "dec")) %>% 
  filter(!grepl("phi", variable),
         str_detect(variable, "^b_")) %>% 
  mutate(summary = sprintf("%.2f (%.2f, %.2f)", mean, `2.5%`, `97.5%`),
         sign = (`2.5%` > 0 & `97.5%` > 0) | (`2.5%` < 0 & `97.5%` < 0) ) %>% 
  select(c(variable,summary, sign)) %>%
  mutate(
    category = case_when(
      str_detect(variable, "_zoi") ~ "zoi",
      str_detect(variable, "_coi") ~ "coi",
      TRUE ~ "none"
      ),
    # remove prefixes
    variable = str_remove_all(variable, paste(remove, collapse = "|")) ) %>%
  # reshape table
  pivot_wider(
    names_from = category,
    values_from = c(summary, sign),
    names_glue = "{category}_{.value}"
  ) %>% 
  # rename variable labels
  mutate(variable = coalesce(mod_labs[variable], variable),
         variable = str_remove(variable, "^b_")) %>% 
  # reorder columns
  select(variable, zoi_summary, coi_summary, none_summary, 
         zoi_sign, coi_sign, none_sign) %>% 
  # reorder rows to match Figure 3
  arrange(desc(row_number())) %>%
  { .[c(2, 1, 3:nrow(.)), ] } %>% 
  # format as gtable
  gt() %>% 
  fmt_markdown(columns = variable) %>% 
  tab_style(
    style = cell_text(v_align = "top", weight = 'bold'),
    locations = cells_column_labels()) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = zoi_summary,
      rows = zoi_sign == TRUE
    )
  ) %>% 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = coi_summary,
      rows = coi_sign == TRUE
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = none_summary,
      rows = none_sign == TRUE
    )
  ) %>%
  cols_hide(
    columns = c(zoi_sign, coi_sign, none_sign)
  ) %>% 
  # rename column labels
  cols_label(
    variable = "Predictor",
    zoi_summary = "Component 1",
    coi_summary = "Component 2",
    none_summary = "Component 3",
  ) %>% 
  gtsave(filename = "tables/Table_S11_oldmod.html")
  

