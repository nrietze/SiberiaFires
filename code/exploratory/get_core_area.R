library(terra)
library(tidyterra)
library(cowplot)
library(tidyverse)
library(landscapemetrics)
library(motif)
library(scico)

ba <- rast('data/geodata/raster/burned_area/planet/Berelech_2020_burned_area.tif')
ba <- rast('data/geodata/raster/burned_area/planet/Kosukhino_2020_burned_area.tif')
check_landscape(ba)

plot(ba)

patch_area <- lsm_p_area(ba)
patch_area$area_m2 <- patch_area$value * 1e4

patch_area %>% 
  ggplot(aes(x = area_m2,color = factor(class))) +
  stat_ecdf(geom = "line", pad = FALSE, 
            aes(y = after_stat(1-y)), size = 1) +
  lims(x = c(0,400)) +
  scale_color_manual(values = c("green", "black"), labels = c("unburned", "burned")) +
  labs(x = "Patch area (m²)", y = "Inverse Cumulative Proportion",color = NULL) +
  theme_cowplot()

ba_perc <- lsp_signature(ba,type = "composition", threshold = 1, window = 33)
ba_perc2 <- lsp_restructure(ba_perc)
ba_perc_rast <- lsp_add_terra(ba_perc2)

ggplot() +
  geom_spatraster(data=ba_perc_rast['X2']) +
  scale_fill_scico(palette = 'buda') +
  theme_cowplot()

writeRaster(ba_perc_rast,
            filename = "data/geodata/raster/burned_area/planet/AOI_2_burn_pct.tif",
            overwrite = T
)
