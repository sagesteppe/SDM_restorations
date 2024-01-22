setwd('/media/steppe/hdd/SDM_restorations/scripts')

library(tidyverse)
library(sf)
library(terra)

riv <- st_read('../data/processed/rivers_noram/rivers_noram_37341.shp', quiet = TRUE)|>
  select(MAJ_NAME, Strahler) |>
  filter(Strahler >= 3)

domain <- rast(nrows = 1, ncols = 1) # create a big empty raster, you can go in through sf too. 
ext(domain) <- c( -125.5, -100, 27,  50) # set the extent
crs(domain) <- "EPSG:4326"
riv <- st_crop(riv, project(domain, crs(riv))) %>% 
  st_simplify() %>% 
  st_transform(5070) %>% 
  st_simplify(dTolerance = 500) |>
  st_buffer(2000)

# import slope and identify large features... 

p <- '~/Documents/WesternPlantPredictors2/slope'
paths <- file.path(p, list.files(p))
slope <- mosaic(
  sprc(
    lapply(paths, rast)
    )
  )

# burn rivers to their own raster

rivr <- terra::rasterize(vect(riv), rast(slope), touches = TRUE, update = TRUE)
plot(rivr)

writeRaster()
