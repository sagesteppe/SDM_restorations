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

##################################################
### now burn lakes and oceans into the raster ###
##################################################

lakes10 <- rnaturalearth::ne_download(scale = 10, type = "lakes_north_america",
                                      category = "physical") %>% 
  st_as_sf() %>% 
  st_make_valid() %>%
  filter(featurecla == 'Lake') %>% 
  select(featurecla)

lakes10 <- lakes10[ which(st_is_valid(lakes10) == T), ]
lakes10 <- st_crop(lakes10, domain)

ocean <- rnaturalearth::ne_download(scale = 10, type = "ocean",
                                    category = "physical") 

rivr <- terra::mask(rivr, vect(ocean), touches = TRUE, update = TRUE)
rivr <- terra::mask(rivr, vect(lakes10), touches = TRUE, update = TRUE)

plot(lakes10[,1])

##############################################################
## read and watershed boundaries and convert some to raster ##
## to delineate groups (spatial structured populations, population, deme) ##
p2gdb <- '../../Geospatial_data/WBD/WBD_National_GDB.gdb'
st_layers(p2gdb)

hu10 <- st_read(p2gdb, layer = 'WBDHU10', quiet = TRUE)  %>% 
  st_make_valid() %>% 
  select(huc10) %>% 
  st_crop(., project(domain, crs(.))) %>% 
  st_simplify() %>% 
  st_transform(5070) %>% 
  st_simplify(dTolerance = 80)

hu10_line <- st_cast(hu10 , 'LINESTRING')

hu12 <- st_read(p2gdb, layer = 'WBDHU12', quiet = TRUE)  %>% 
  st_make_valid() %>% 
  st_crop(., project(domain, crs(.))) %>% 
  st_simplify() %>% 
  st_transform(5070) %>% 
  st_simplify(dTolerance = 80)

hu12_line <- st_cast(hu12 , 'LINESTRING')

hu10R <- terra::rasterize(vect(hu10), rast(slope), update = TRUE, field = 'huc10', 
                          filename = '../../Geospatial_data/WBD_rast/hu10R.tif', overwrite = T)

terra::rasterize(vect(hu10_line), rast(slope), update = TRUE, field = 'huc10', 
                 filename = '../../Geospatial_data/WBD_rast/hu10R-line.tif', overwrite = T)

terra::rasterize(vect(hu12), rast(slope), update = TRUE, field = 'huc12', 
                          filename = '../../Geospatial_data/WBD_rast/hu12R.tif', overwrite = T)

terra::rasterize(vect(hu12_line), rast(slope), update = TRUE, field = 'huc12', 
                 filename = '../../Geospatial_data/WBD_rast/hu12R-line.tif', overwrite = T)
