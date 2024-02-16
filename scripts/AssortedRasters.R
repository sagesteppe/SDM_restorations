setwd('/media/steppe/hdd/SDM_restorations/scripts')

library(tidyverse)
library(terra)
library(sf)

mrlc <- "/media/steppe/hdd/Geospatial_data/MRLC/"
blm <- vect('/media/steppe/hdd/Geospatial_data/Public_land_ownership/BLM_West_simplified/BLM_West_simplified.shp')
list.files(mrlc, recursive = T)

change <- rast(paste0(mrlc, 'total_change_intensity_index/rcmap_total_change_intensity_index.tif'))
terra::aggregate(change, fact = 5, fun = 'mean', cores = parallel::detectCores(), overwrite = T, 
                 filename = paste0(mrlc, 'total_change_intensity_index/TCII_3arc.tif'))
change <- rast(paste0(mrlc, 'total_change_intensity_index/TCII_3arc.tif'))
change <- mask(change, blm)
writeRaster(change, filename = paste0(mrlc, 'total_change_intensity_index/TCII_3arc.tif'), overwrite = T)

rm(change)

may_inv <- rast(
  paste0(mrlc, 'ExoticAnnualGrass_May152023_PercentCover/ExoticAnnualGrass_May152023_PercentCover.img'))
may_inv <- ifel(may_inv > 100, NA, may_inv)

june_inv <- rast(
  paste0(mrlc, 'ExoticAnnualGrass_June192023_PercentCover/ExoticAnnualGrass_June192023_PercentCover.img'))
june_inv <- ifel(june_inv > 100, NA, june_inv)

plot(june_inv)

inv_grass <- c(may_inv, june_inv)
inv_grass$twentyThreeMean <- app(inv_grass, max, na.rm = F)
twentyThreeMean <- inv_grass$twentyThreeMean

terra::aggregate(twentyThreeMean, fact = 5, fun = 'mean', cores = parallel::detectCores(), overwrite = T, 
                 filename = paste0(mrlc, '2023_InvGrassCover.tif'))

inv_g <- rast(paste0(mrlc, '2023_InvGrassCover.tif'))
writeRaster(inv_g, filename = paste0(mrlc, '2023_InvGrassCover.tif'), overwrite = T)

inv_g <- mask(inv_g, blm)

rm(june_inv, may_inv, inv_grass, twentyThreeMean, mrlc)

###############################################
# Create Roads Data set

library(tigris)

states <- states() %>% 
  filter(STUSPS %in% c('NV', 'UT', 'CA', 'AZ', 'CO')) %>% 
  select(STUSPS, ST_NAME = NAME, STATEFP) %>% 
  st_drop_geometry()

county <- counties(states$STUSPS) %>% 
  st_transform(5070)
county <- county[ lengths( st_intersects(county, blm)) > 0, ]

county_rds <- roads(county$STATEFP[1], county$COUNTYFP[1]) %>% 
  st_transform(st_crs(blm))

p2 <- '/media/steppe/hdd/Geospatial_data/blmRoadsSOS24/raw/'

road_setter <- function(x){
  
  statefp <- x$STATEFP[1]
  county_rds <- lapply(X = x$STATEFP, county = x$COUNTYFP, FUN = tigris::roads) %>% 
    dplyr::bind_rows() %>% 
    sf::st_transform(5070)
  blm_roads <- sf::st_intersection(county_rds, blm)
  
  st_write(blm_roads, paste0(p2, statefp, '.shp'))
}

rs <- split(county, f = county$STATEFP)
lapply(rs, road_setter)

roads <- lapply(paste0(p2, list.files(p2, pattern = 'shp$')), sf::st_read, quiet = TRUE) |>
  dplyr::bind_rows() %>% 
  select(LINEARID)

roads_simp <- st_simplify(roads, dTolerance = 10)
object.size(roads_simp)/10^9 # ~55% of original size. 
st_write(roads_simp, '/media/steppe/hdd/Geospatial_data/blmRoadsSOS24/processed/BLM_roads.shp', append = F)

rm(states, county, p2, roads, roads_simp, rs, road_setter)



## simplify boundaries of allotments

allottments <- st_read(paste0(p, 
  'Public_land_ownership/BLM_Natl_Grazing_Allotment_Polygons/BLM_Natl_Grazing_Allotment_Polygons.shp')) %>% 
  filter(ADMIN_ST != 'AK')

allottments15.5 <- rmapshaper::ms_simplify(allottments, keep = 0.2, weighting = 0.6, keep_shapes = TRUE, sys = TRUE)
allottments15.5 <- st_make_valid(allottments15.5)
st_write(allottments15.5, paste0(p, 'Public_land_ownership/BLM_Natl_Grazing_Allotment_Polygons/Allotment-simp.shp'), append = F)
#### CREATE boundaries of field offices. 

p <- '/media/steppe/hdd/Geospatial_data/'
p2gdb <- '/Public_land_ownership/BLM_National_Administrative_Unit_B.gdb'
bounds <- st_read(paste0(p, p2gdb), layer = 'blm_natl_admu_field_poly_webpub')  %>% 
  filter(ADMIN_ST != 'AK') %>% 
  select(ADMU_NAME, PARENT_NAME) 

bounds <- bounds %>% 
  mutate(
    across(ADMU_NAME:PARENT_NAME, ~ str_to_upper(.x)),
    ADMU_NAME = paste(str_remove(ADMU_NAME, ' FIELD.*$'), 'FIELD OFFICE'), 
    PARENT_NAME = paste(str_remove(PARENT_NAME, ' DISTRICT.*$'), 'DISTRICT OFFICE'),
    ADMU_NAME = str_replace(ADMU_NAME, 'ROSEBURG DISTRICT SWIFTWATER FO FIELD OFFICE', 
                            'ROSEBURG DISTRICT SWIFTWATER FIELD OFFICE'), 
    ADMU_NAME = str_replace(ADMU_NAME, 'ROSEBURG DISTRICT SOUTH RIVER FO FIELD OFFICE', 
                            'ROSEBURG DISTRICT SOUTH RIVER FIELD OFFICE'), 
    ADMU_NAME = str_replace(ADMU_NAME, 'PALM SPRINGS/S. COAST FIELD OFFICE', 
                            'PALM SPRINGS/SOUTH COAST FIELD OFFICE'), 
    PARENT_NAME = str_replace(PARENT_NAME, 'DAKS', 'DAKOTAS')
    )

st_write(bounds, paste0(p, 'Public_land_ownership/BLM_ADMU/BLM_ADMU_BOUNDARIES.shp'))

ggplot() +
  geom_sf(data = bounds)

head(bounds)
