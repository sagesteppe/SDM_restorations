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

r1 <- rast(f[1])
r1 <- terra::ifel(r1 < 0.5, NA, r1)
ng1 <- terra::crop(noGrow, r1)
plot(r1)

terra::mask(r1, ng1, maskvalues = 1)

phenAggregator <- function(x){
  
  taxon <- basename(x)
  r <- terra::rast(x)
  r <- terra::ifel(r < 0.5, NA, r)
  r <- terra::resample(r, template, threads = parallel::detectCores())
  terra::mask(r, terra::crop(noGrow, terra::ext(r)), inverse = T, 
                   filename = file.path('../results/PhenPredSurfaces', taxon))
  
}

lapply(f, phenAggregator)
