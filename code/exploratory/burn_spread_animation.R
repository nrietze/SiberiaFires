library(tidyverse)
library(ggmap)
library(gganimate)
library(sf)
library(magick)
library(basemaps)
library(scico)
scico_palette_show()
register_stadiamaps('1e687c0d-c7c1-4dec-9d98-29829a90697d') 

FNAME <- 'C:/data/3_fire_data/active_fires/modis/fire_archive_modis_cavm.csv'
fire_data <- read.csv(FNAME) %>% na.omit()

OUT_PATH <- 'C:/Users/nrietze/Documents/1_PhD/8_CHAPTER2/output/fire_spread_animation/'

doys <- seq(min(fire_data$doy), max(fire_data$doy))

# all kytalyk fires
xmin <- 143.83
ymin <- 70.43
xmax <- 148.96
ymax <- 71.67

# large fire only
# xmin <- 143.84
# ymin <- 71.2
# xmax <- 146.45
# ymax <- 71.63

area_coords <- c(
  xmin, ymin, xmax, ymax
)

terrain_bg <- F

if (terrain_bg == TRUE){
  bg_layer <- ggmap::get_stadiamap(
    bbox = area_coords,
    zoom = 8,
    maptype = "stamen_terrain",
    color = "bw",
    force = T
  )
} else {
  # load handdrawn extent as ext_draw
  load('C:/Users/nrietze/Documents/1_PhD/8_CHAPTER2/SiberiaFires/data/geodata/feature_layers/animation_extent.rda')
  
  ext <- st_bbox(c(xmin = xmin, 
                   ymin = ymin,
                   xmax = xmax, 
                   ymax = ymax), 
                 crs = st_crs(4326)) %>% 
    st_as_sfc() 
  # %>%
  #   st_transform(crs = 3857) %>%
  #   st_bbox()
  
  map <- basemap_ggplot(ext_draw, map_service = "esri", map_type = "world_imagery",map_res = .2)
}

plot_list <- list()

i <- 1
for (doy_i in 150:244){
  if (doy_i %in% fire_data$doy) {
    data <- fire_data %>% 
      filter(doy <= doy_i) %>% 
      mutate(time_since_first = doy_i - doy) %>% 
      mutate(time_since_first = ifelse(time_since_first > 7, 8, time_since_first))  
    
    lab_date <- as.Date(paste0("2020-", doy_i), format = "%Y-%j")
    
    
    # Plot with stamen terrain background
    if (terrain_bg == TRUE){
      g <- ggmap(bg_layer) +
        geom_point(data = data, aes(x = LONGITUDE, y = LATITUDE, color = time_since_first)) +
        xlim(xmin, xmax) +  # Set x axis limits, xlim(min, max)
        ylim(ymin, ymax) +
        # color_scale +
        scico::scale_color_scico(palette = 'lajolla') +
        labs(title = lab_date) +
        theme_nothing() +
        theme(plot.title = element_text(
          size = 16, color = "grey10",
          hjust = .5
        ))
    } else {
      
    # Plot with ESRI satellite imagery background
      data <- data %>% 
        st_as_sf(coords = c('LONGITUDE','LATITUDE'),
                 crs = 4326) %>% 
        st_transform(3857)
    
      ( g <- map + 
          geom_sf(data = data, aes(color = time_since_first),stat = 'sf') +
          scico::scale_color_scico(palette = 'lajolla') +
          labs(title = lab_date) +
          theme_nothing() +
          theme(plot.title = element_text(
            size = 16, color = "grey10",
            hjust = .5
          )) + 
          xlim(xmin, xmax) +  # Set x axis limits, xlim(min, max)
          ylim(ymin, ymax) )
      }
      
    # Add plots to list for export
    plot_list[[i]] <- g
    ggsave(paste0(OUT_PATH,"plot_", i, ".png"),
           plot = g, width = 6, height = 6, dpi = 300)
    i <- i+1
  }
}

# convert to GIF
png_files <- list.files(OUT_PATH, pattern = "*.png", full.names = TRUE)
png_files <- gtools::mixedsort(sort(png_files))

# Read in PNG files as image objects
png_images <- lapply(png_files, image_read)

# Combine PNG images into a GIF
gif <- image_join(png_images)

# Set GIF delay and loop options
gif <- image_animate(gif, fps = 10, dispose = "previous", loop = 0)

# Write GIF to file
image_write(gif, paste0(OUT_PATH,"kytalyk_spread_satellite.gif"))
