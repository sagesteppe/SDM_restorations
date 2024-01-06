setwd('/home/sagesteppe/Documents/SDM_restorations/scripts')

library(sf)
library(tidyverse)
library(terra)

boundaries <- read_sf('../data/raw/geomorph90_Boundary.kmz') %>% 
  st_transform(4326)

domain <- rast()
domain <- rast(nrows = 1, ncols = 1) # create a big empty raster, you can go in through sf too. 
ext(domain) <- c( -125.5, -100, 27,  50) # set the extent
crs(domain) <- "EPSG:4326"
domain <- as.polygons(domain) %>% 
  st_as_sf()

weird <- data.frame(
  ID = c(rep(LETTERS[1:2], each = 5)),
  X = c(-105, -100, -100, -105, -105, -105, -100, -100, -105, -105),
  Y = c(40, 40, 35, 35, 40, 35, 35, 30, 30, 35)
) %>% 
  st_as_sf(coords = c('X', 'Y'), crs = 4326) %>% 
  group_by(ID) %>%  
  dplyr::summarise() %>%
  st_cast("POLYGON")

poly <- bind_rows(
  boundaries[ lengths( st_overlaps(boundaries, domain)) > 0, ],
  boundaries[ unlist(st_contains(domain, boundaries)), ],
  weird
) %>% 
  select(geometry)  %>% 
  mutate(ID = 1:n())

ggplot() +
  geom_sf(data = poly) +
  geom_sf(data = domain, fill = NA, color = 'red')

