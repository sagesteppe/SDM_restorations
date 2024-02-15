setwd('/media/steppe/hdd/SDM_restorations/scripts')

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


domain <- rast(nrows = 1, ncols = 1) # create a big empty raster, you can go in through sf too. 
ext(domain) <- c( -125.5, -100, 27,  50) # set the extent
crs(domain) <- "EPSG:4326"


p <- '/media/steppe/hdd/Geospatial_data/Transportation_National_GDB/Transportation_National_GDB.gdb'
roads <- st_read(p, layer = 'Trans_RoadSegment', quiet = TRUE) 

st_crop(roads, domain) # restrict to the west
# now restrict to BLM 

  
# mask away urban areas. 


