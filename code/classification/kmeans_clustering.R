library(sp)
library(terra)
library(raster)
library(ggplot2)
library(viridis)
library(rasterVis)
library(scico)

wv_aoi1_2022 <- brick('C:/data/7_worldview/2022/aoi_1/015284154010_01_P001_MUL/22JUL21020532-M2AS-015284154010_01_P001.TIF')
wv_aoi1_2022

hist(wv_aoi1_2022[[8]],
     main = "Distribution of reflectance",
     xlab = "NDVI",
     ylab= "Frequency",
     col = "aquamarine3",
     xlim = c(0, 1200),
     breaks = 30,
     xaxt = 'n')
# axis(side = 1, at = seq(-0.5,1, 0.05), labels = seq(-0.5,1, 0.05))

  
VI <- function(img, k, i) {
  bk <- img[[k]]
  bi <- img[[i]]
  vi <- (bk - bi) / (bk + bi)
  return(vi)
}

# NDVI, NIR = 8, Red = 5
ndvi <- VI(wv_aoi1_2022, 8, 5)

plot(ndvi, col = scico(10, palette = 'cork'), main = 'Worldview-3 NDVI')


hist(ndvi,
     main = "Distribution of NDVI values",
     xlab = "NDVI",
     ylab= "Frequency",
     col = "aquamarine3",
     xlim = c(-0.5, 1),
     breaks = 30,
     xaxt = 'n')
axis(side = 1, at = seq(-0.5,1, 0.05), labels = seq(-0.5,1, 0.05))

veg <- reclassify(ndvi, cbind(-Inf, 0, NA))
plot(veg, col = rev(scico(10, palette = 'bamako')), main = 'Veg cover')

nr <-getValues(ndvi)
nr_nonan <- na.omit(nr)

set.seed(99)

# create 10 clusters, allow 500 iterations, start with 5 random sets using 'Lloyd' method

kmncluster <- kmeans(nr_nonan, centers = 5, iter.max = 500,
                     nstart = 5, algorithm = "Lloyd")

str(kmncluster)
cluster_labels <- kmncluster$cluster

new_raster <- ndvi
new_raster[] <- NA

# Assign the cluster labels to the new raster layer
new_raster[!is.na(ndvi)] <- cluster_labels

par(mfrow = c(1, 2))
plot(ndvi, col =rev(scico(10, palette = 'bamako')), main = "NDVI")
plot(new_raster, main = "Kmeans", col = viridis_pal(option = "D")(5))

par(mfrow = c(1, 1))
plot(new_raster, main = "Kmeans", col = viridis_pal(option = "D")(5))