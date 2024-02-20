setwd('/media/steppe/hdd/SDM_restorations/scripts')

library(terra)
source('functions.R')

p <- '../results/suitability_maps'
f <- file.path(p, list.files(p, pattern = '1k-'))

template <- rast('../data/processed/gddPCA.tif') |>
  project(crs(rast(f[1])))
noGrow <- rast('../data/processed/nlcd.tif') |>
  resample(template) |>
  project(crs(template))

phenAggregator <- function(x){
  
  taxon <- basename(gsub('1k-.*$', '', x))
  r <- terra::rast(x)
  r <- terra::ifel(r < 0.5, NA, 1)
  r <- terra::resample(r, template, threads = parallel::detectCores())
  terra::mask(r, terra::crop(noGrow, terra::ext(r)), inverse = T, overwrite = T, 
                   filename = file.path('../results/PhenPredSurfaces', paste0(taxon, '.tif')))
  
  terra::tmpFiles(current = FALSE, orphan = TRUE, old = TRUE, remove = TRUE)
}

lapply(f, phenAggregator)
