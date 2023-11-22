library(tidyverse)
library(scales)
library(sf)
library(terra)
library(scico)
scico_palette_show()

# Load MODIS hotspot archive data
FNAME <- 'C:/data/3_fire_data/active_fires/modis/fire_archive_modis_cavm.csv'
fire_data <- read.csv(FNAME) %>% 
  na.omit() %>% 
  st_as_sf(., coords = c("LONGITUDE","LATITUDE")) %>% 
  vect()

# Load areas of interest
aois <- vect('./data/geodata/feature_layers/aoi_wv/aois.shp') %>%
  mutate(., id = 1:nrow(.))

# Get fire points from 2020 that intersect with AOIs
fires_aoi <- fire_data %>% 
  terra::intersect(.,aois) %>% 
  mutate(ACQ_DATE = dmy_hms(ACQ_DATE)) %>% 
  group_by(id) %>% 
  mutate(first_day = min(ACQ_DATE)) %>% 
  mutate(last_day = max(ACQ_DATE))

first_day <- min(fires_aoi$ACQ_DATE)
last_day <- max(fires_aoi$ACQ_DATE)

# Plot
ggplot(fires_aoi) +
  geom_segment( aes(x=site, xend=site, y=first_day, yend=last_day), color="grey") +
  geom_point(aes(x=site, y=ACQ_DATE) , colour=rgb(0.7,0.2,0.1,1), size=3 ) +
  scale_y_datetime(labels = date_format("%b %d"),
                   limits = as.POSIXct(strptime(c("2020-06-10", "2020-07-30"), 
                                                format = "%Y-%m-%d")),
                   breaks = date_breaks("1 week")) +
  annotate("text", x = 2.2, y = first_day, label = first_day) + 
  annotate("segment", x = 2.15, xend = 2, y = first_day+days(1), yend = first_day, colour = "black") +
  annotate("text", x = 3.5, y = last_day, label = last_day) + 
  annotate("segment", x = 3.55, xend = 4, y = last_day+days(1), yend = last_day, colour = "black") +
  coord_flip()+
  theme_cowplot() +
  theme(
    legend.position = "none",
  ) +
  xlab("") +
  ylab("Acquisition date")

ggsave('./figures/fire_durations.png')


