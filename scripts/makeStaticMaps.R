library(tidyverse)
library(sf)
library(lubridate)
library(terra)
library(tigris)
library(ggthemes)

domain <- rast(nrows = 1, ncols = 1) # create a big empty raster, you can go in through sf too. 
ext(domain) <- c( -125.0, -100, 27, 50) # set the extent
crs(domain) <- "EPSG:4326" # define the projection
domain <- project(domain, "EPSG:5070") |>
  as.polygons() |>
  st_as_sf()
#domain <- st_bbox(domain)

setwd('/media/steppe/hdd/SDM_restorations/results')
p <- '/media/steppe/hdd/SDM_restorations/results/suitability_maps'
f <- file.path(p, list.files(p))
f <- f[grep('1k-2024', f)]

f1 <- rast(f[1])
st_bbox(f1)

st <- tigris::states() %>% 
  filter(REGION == '4' | STUSPS %in% c('SD', 'ND', 'OK', 'NE', 'TX')) %>% 
  filter(! STUSPS %in% c('AK', 'HI')) %>% 
  st_transform(st_crs(f1)) %>% 
  st_intersection(., domain)

ocean <- rnaturalearth::ne_download(scale = 10, type = "ocean", category = "physical") %>% 
  st_as_sf() %>% 
  st_transform(st_crs(f1)) %>% 
  st_intersection(., domain)

occ_p <- '/media/steppe/hdd/SDM_restorations/data/raw/occurrence/jittered_occurrences'
occurrences <- st_read(file.path(occ_p, 'jittered_occurrences.shp'), quiet = T)

dem <- mosaic(
  sprc(
    lapply(file.path('/media/steppe/hdd/Geospatial_data/GMTED',
                     list.files('/media/steppe/hdd/Geospatial_data/GMTED',
                                recursive = T, pattern = 'mea300.tif'))
           , rast))
  ) |> 
  project(crs(domain)) 

dem <- crop(dem, domain)
dem <- mask(dem, vect(ocean), inverse = T)
slope <- terrain(dem, "slope", unit="radians")
aspect <- terrain(dem, "aspect", unit="radians")
hill <- shade(slope, aspect, 30, 315)

rm(slope, aspect)

st_bbox(st)
#' plot many same sdm maps
#' 
#' @param x a raster stack of species to plot
#' @param y an sf tibble of the same species occurrence data to plot
#' @param path_dir directory to save the outfiles
map_maker_fn <- function(x, y, path_dir){
  
  spp <- rast(x)
  spp <- aggregate(spp, 5, fun = 'max', na.rm=TRUE)
  spp <- ifel(spp < 0.5, NA, spp)
  
  dem_sub <- crop(dem, spp)
  hill_sub <- crop(hill, spp)
  dem <- as.data.frame(dem_sub, xy = T)
  colnames(dem) <- c('x', 'y', 'elevation')
  hillshade <- as.data.frame(hill_sub, xy = T)
  
  focal_taxon <- gsub('1k.*$', '', gsub('_', ' ', basename(x)))

  occurrence <- filter(y, taxon == focal_taxon) %>% 
    st_crop(st_bbox(spp))
  pred <- as.data.frame(spp, xy= T)
  
  st_sub <- st_crop(st, st_bbox(spp))
  ocean_sub <- st_crop(ocean, st_bbox(spp))
  
  ggplot() + #3C6997
    geom_sf(data = ocean_sub, fill = '#3C6997') + # #97DFFC
    geom_tile(data = hillshade, aes(x = x, y = y, fill = hillshade))  +
    scale_fill_gradient(low = "grey50", high = "grey100") +
    guides(fill = 'none') +
    ggnewscale::new_scale_fill() + 
    
    geom_tile(data = dem, aes(x = x, y = y, fill = elevation)) + 
    scale_fill_viridis_c(direction = -1, alpha = 0.3) + 
    guides(fill = 'none') +
    ggnewscale::new_scale_fill() + 
    
    geom_sf(data = occurrence, shape = 5, alpha = 0.8) + 
    geom_tile(data = pred, aes(x = x, y = y, fill = predicted_suitability), alpha =0.9) + 
    scale_fill_distiller('Probability of suitable habitat      ', 
                        palette = "RdPu", direction = 1, 
                         breaks=c(0.625, 0.75, 0.875), limits = c(0.5, 1.0),
                         labels=c("Mild", "Medium", "High")) + 
    
    geom_sf(data = st_sub, fill = NA, lwd = 1) + 
    
    theme_void() + 
    
    theme(legend.position = 'bottom', 
          plot.background = element_rect(
            fill = "transparent", colour = "#FFC759", linewidth = 1),
          panel.background = element_rect(fill = 'white'),
          legend.box.background = element_rect(
            color='#FFC759', fill = "transparent", linewidth = 1),
          legend.box.margin = margin(2, 40, 2, 40), 
    )
  
  ggsave(filename = file.path(path_dir, paste0(focal_taxon, '.png')), plot = last_plot(), 
         width = 420, height = 420, units = "px",
         dpi = 72, device = 'png', bg = 'transparent')
  
}

lapply(f[7], FUN = map_maker_fn, 
       y = occurrences,
       path_dir = '/media/steppe/hdd/SDM_maps_CBG_guidebooks2024'
       )



breaks <- c(121, 152, 182, 213, 244)
labels <- c('May', 'June', 'July', 'Aug.', 'Sept.')

spliecies <- split(records, records$species)

phen_maker_fn <- function(x){
  
  species <- gsub('_', ' ', x$species[1])
  
  ggplot(x) +
    geom_density(aes(x = DOY), fill = '#FF7B9C', color = '#FFC759', alpha = 0.5) +
    theme_bw() + 
    labs(title = 'Estimated Flowering') + 
    theme(aspect.ratio = 6/16, 
          plot.title = element_text(
            hjust = 0.5, colour = "black", size = 7, face = "bold"),
          axis.text.x= element_text(
            family = "Tahoma", face = "bold", colour = "black", size=5),
          panel.background = element_rect(fill='#607196'),
          plot.background = element_rect(fill='#607196'),
          panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.x = element_blank(), 
          panel.grid.minor.y = element_blank()) +
    scale_x_continuous(breaks = breaks, labels=labels, limits = c(120, 245) )
  
  ggsave(filename = paste0('../results/phen/', species, '.png'), plot = last_plot(), 
         dpi = 150, width = 260, height = 150,  units = "px",  bg = 'transparent')
}