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
  st_simplify(dTolerance = 250) |>
  st_buffer(140)

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
writeRaster(rivr, '../../Geospatial_data/river_resistance/rivR.tif', overwrite = TRUE)



##############################################################
## read and watershed boundaries and convert some to raster ##
## to delineate groups (spatial structured populations, population, deme) ##
p2gdb <- '../../Geospatial_data/WBD/WBD_National_GDB.gdb'
st_layers(p2gdb)

hu10 <- st_read(p2gdb, layer = 'WBDHU10', quiet = TRUE)  %>% 
  st_make_valid() %>% 
  st_crop(., project(domain, crs(.))) %>% 
  st_simplify() %>% 
  st_transform(5070) %>% 
  st_simplify(dTolerance = 80)

hu12 <- st_read(p2gdb, layer = 'WBDHU12', quiet = TRUE)  %>% 
  st_make_valid() %>% 
  st_crop(., project(domain, crs(.))) %>% 
  st_simplify() %>% 
  st_transform(5070) %>% 
  st_simplify(dTolerance = 80)

hu10 <- select(hu10, huc10)
hu12 <- select(hu12, huc12)

hu10R <- terra::rasterize(vect(hu10), rast(slope), update = TRUE, field = 'huc10')
hu12R <- terra::rasterize(vect(hu12), rast(slope), update = TRUE, field = 'huc12')

writeRaster(hu10R, '../../Geospatial_data/WBD_rast/hu10R.tif', overwrite = TRUE)
writeRaster(hu12R, '../../Geospatial_data/WBD_rast/hu12R.tif', overwrite = TRUE)

