setwd('~/Documents/SDM_restorations/scripts')
getwd()

library(tidyverse)
library(sf)
library(terra)

domain <- rast(nrows = 1, ncols = 1) # create a big empty raster, you can go in through sf too. 
ext(domain) <- c( -125.5, -100, 27, 50) # set the extent
crs(domain) <- "EPSG:4326" # define the projection
domain <- as.polygons(domain) |>  # convert to vector data
  st_as_sf() 

rivers <- st_read('../data/processed/rivers_noram/rivers_noram_37341.shp', quiet = TRUE) %>% 
  st_crop(., domain) %>% 
  select(Strahler, SUB_NAME, SUBBAS_ID) %>% 
  filter(Strahler > 2) %>% 
  group_by(SUBBAS_ID) %>% 
  reframe(geometry = st_union(geometry)) %>% 
  st_as_sf()

rivers <- rivers[ as.numeric(st_length(rivers)) > 30000, ]
rivers <- st_simplify(rivers) %>% 
  st_transform(5070) %>% 
  st_buffer(1000) %>% 
  st_transform(4326)

lakes10 <- rnaturalearth::ne_download(scale = 10, type = "lakes_north_america",
                       category = "physical") %>% 
  st_make_valid()

lakes10 <- lakes10[ which(st_is_valid(lakes10) == T), ]
lakes10 <- st_crop(lakes10, domain)


ocean <- rnaturalearth::ne_download(scale = 10, type = "ocean",
                                    category = "physical") 

ggplot() + 
  geom_sf(data = lakes10)  + 
  geom_sf(data = rivers) +
  theme_bw()

domain <- rast(nrows = 25000, ncols = 25000) 
# create a big empty raster, you can go in through sf too. 
ext(domain) <- c( -125.5, -100, 27, 50) # set the extent
crs(domain) <- "EPSG:4326"
values(domain) <- 1

domain <- terra::mask(domain, vect(ocean), inverse = T)
domain <- terra::mask(domain, vect(lakes10), inverse = T)
domain <- terra::mask(domain, vect(rivers), inverse = T, touches = TRUE)

plot(domain)
