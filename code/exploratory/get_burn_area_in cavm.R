library(sf)
library(tidyverse)
library(terra)
library(raster)

library(scico)
scico_palette_show()


# Load CAVM file
cavm <- read_sf('C:/data/6_vegetation/cavm/vegetation_zones.shp')

# Load raster data
CCI_PATH <- 'C:/data/3_fire_data/burned_area/fire_cci/2020/'
OUT_PATH <- "./output/cavm_vs_burned_area/"

ba <- raster(paste0(CCI_PATH,'20200401-ESACCI-L3S_FIRE-BA-MODIS-AREA_1-fv5.1-JD.tif'))

ba_in_cavm <- read.csv(paste0(OUT_PATH,'Arctic_BA_FireCCI51.csv'), 
                        header = TRUE, sep = ";", colClasses = c("factor", rep("numeric", 7))) 

n_zones <- 7
colors <- rev(scico(n_zones, palette = "batlowW"))

ba_in_cavm %>% 
  pivot_longer(cols = c(2:8), 
               names_to = 'Zone_name',
               values_to = 'area') %>% 
  ggplot(aes(x = Year, y = area, fill = Zone_name)) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values = colors) +
  xlab("") + ylab("Burned area (ha)") +
  labs(fill = NULL,
       title = "Fire CCI burned area within CAVM") +
  theme_classic() +
  theme(legend.position = "bottom")
# ggsave(paste0(OUT_PATH, "ba_cavm_relative_tile_17_ls.png") )


#  Plot burned fractions per CAVM type ----
preds <- rast("data/geodata/raster/predictors/Libenchik_predictors_30m.tif")
bp <- vect('data/geodata/feature_layers/burn_polygons/planet/rough_burn_perimeter_Libenchik.shp')
cavm <- rast('C:/data/6_vegetation/cavm/raster_cavm_v1_NE_sib.tif') %>% 
  crop(ext(preds)) %>% 
  as.polygons() %>% 
  crop(bp)

terra::extract(preds$burned_fraction,cavm) %>% 
  mutate(CAVM_type = factor(ID,labels = c('B1','G3','G4','S1','W3','FW'))) %>% 
  ggplot(aes(x = CAVM_type,y = burned_fraction,fill = CAVM_type)) +
  geom_violin(position = position_nudge(x = 0.2, y = 0), alpha = 0.8) +
  geom_point(aes(y = burned_fraction, color = CAVM_type),
             position = position_jitter(width = 0.15), size = 1, alpha = 0.1) +
  scale_fill_manual(values = c('B1' = '#70AD47','G3' = '#4472C4', # same colors as in Excel figure
                               'G4' = '#FFC000','S1' = '#43682B',
                               'W3' = '#997300','FW' = '#9DC3E6')) +
  scale_color_manual(values = c('B1' = '#70AD47','G3' = '#4472C4',
                                'G4' = '#FFC000','S1' = '#43682B',
                                'W3' = '#997300','FW' = '#9DC3E6')) +
  theme_cowplot() + 
  theme(legend.position = 'none')
