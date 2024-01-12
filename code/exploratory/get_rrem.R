library(terra)
library(rrrem)

dem <- rast('C:/data/4_geodata/arcticDEM/v3/aoi_4_dem_v3.tif')
rem <- make_rem(dem)

writeRaster(rem,
            filename = 'C:/data/4_geodata/arcticDEM/v3/aoi_4_rrem_v3.tif',
            overwrite = T
)

hillshade <- terra::shade(
  terra::terrain(dem, "slope", unit = "radians"),
  terra::terrain(dem, "aspect", unit = "radians")
)

trimmed_rem <- rem
trimmed_rem[trimmed_rem > 200] <- 200

par(
  oma = c(0, 0, 0, 0),
  mar = c(0, 0, 0, 0),
  mgp = c(0, 0, 0)
)


terra::plot(
  hillshade, 
  col = grey(0:100/100), 
  legend = FALSE, 
  axes = FALSE,
  oma = NA,
  mar = NA,
  xaxs="i", 
  yaxs="i"
)


terra::plot(
  trimmed_rem,
  col = viridis::mako(50, direction = -1), 
  legend = FALSE, 
  axes = FALSE,
  add = TRUE,
  alpha = 0.75
)
