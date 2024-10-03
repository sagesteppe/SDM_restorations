setwd('~/Documents/SDM_restorations/scripts')
library(tidyverse)

p <- '../data/raw/AbsenceSampling'
f <- file.path(p, list.files(p))

PresAbs <- do.call("rbind", sapply(f, read.csv, simplify = FALSE, skip = 4)) |>
  rownames_to_column(var = 'SampleEvent') |>
  mutate(
    SampleEvent = gsub('[.].*$','', basename(SampleEvent)),
    Census = if_else(is.na(Census), 0, Census),
    Pres.Abs = case_when(
      Pres.Abs == 'P' ~ 1, 
      Pres.Abs == 'A' ~ 0, 
      .default = as.numeric(Pres.Abs)
    )) 

ob <- do.call("rbind", sapply(f, read.csv, simplify = FALSE))[,1:3]

readheadR <- function(x){
  
  dat <- read.csv(x, header = F)[1:3,1:2] 
  sampleEvent <- gsub('[A-z]|[.]', '', basename(x))
  
  colnames(dat) <- c('Variable', 'Response')
  Crew <- gsub('CO150', 'UFO',
       gsub('CO130', 'GJFO',
          gsub('UT070|UT70', 'PFO',
              gsub('UT080|UT80', 'VFO', dat$Response))))[3]
  sampleEvent <- paste0(Crew, sampleEvent)
  
  dat <- tidyr::pivot_wider(dat, names_from = 'Variable', values_from = 'Response') 
  
  ob <- data.frame(
    cbind(
      matrix(
        sapply(dat[c('Latitude', 'Longitude')], as.numeric), nrow = 1)
    ), sampleEvent) |>
    setNames(c('Latitude', 'Longitude', 'SampleEvent'))
  
  return(ob)
  
}

SiteData <- lapply(f, readheadR) |>
  bind_rows() |>
  drop_na(Longitude) |>
  sf::st_as_sf(coords = c('Longitude', 'Latitude'))

left_join(PresAbs, SiteData, by = 'SampleEvent') 
