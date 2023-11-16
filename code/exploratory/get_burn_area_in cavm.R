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
