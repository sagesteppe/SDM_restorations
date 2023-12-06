absence_drawer <- function(x, bg_abs){
  
  
  # create the observed hull, the buffered hull for random sample selection, 
  # and a hull to form the basis for a bounding box to generate regular absences in
  
  observed_hull <- sf::st_union(x) |>
    sf::st_convex_hull() |>
    sf::st_sfc() 
  
  obs_area <- as.numeric(st_area(observed_hull)) * 0.15 / 2
  obs_area <- if(obs_area > 160934){obs_area = 160934} # cap area with 100 mi buffer. 
  
  buffer_hull <- st_buffer(observed_hull, obs_area) 
  
  buf_area <- as.numeric(sf::st_area(buffer_hull) *0.1) / 2
  obs_area <- if(buf_area > 40233){buf_area = 40233} # cap area with 25 mi buffer.
  
  box_hull <- sf::st_buffer(buffer_hull, dist = obs_area)
  
  bbox <- sf::st_bbox(box_hull) |>
    sf::st_as_sfc() 
  
  box_hull <- sf::st_difference(bbox, buffer_hull)
  
  # perform the outline background absence sample
  reg_pts <- sf::st_sample(box_hull, size = round(nrow(x) * bg_abs, 0), type = 'regular')
  
  rm(box_hull, bbox, observed_hull)
  

  # identify AIM plots and generate absences. 
  
  blm_surf_area <- as.numeric(sf::st_area(blm_surf))
  prcnt_blm <- (blm_surf_area / buf_area) * 100
  aim_recs <- (1 - bg_abs) * round(prcnt_blm, 0) * 0.1
  
  AIM_points <- AIM_points[lengths(st_intersects(AIM_points, buffer_hull)) > 0,]
  AIM_points <- AIM_points |> 
    filter(! PrimaryKey %in% x$PrimaryKey ) # these are presences
  
  # generate AIM plot absences
  
  
#  st_intersection(aim_plots, blm_admin_surf)
  
  
  # generate Random absences
  
#  rand_recs <- ((1 - bg_abs) * (1 - prcnt_blm)) * 0.1
#  buffered_presences <- sf::st_buffer(x, dist = 10000) |>
#    sf::st_union()
  
#  non_blm_area <- sf::st_difference(blm_admin_surf, buffer_hull) 
#  non_blm_area <- sf::st_difference(non_blm_area, buffered_presences)
  
#  random_samples <- sf::st_sample(non_blm_area, prcnt_blm)
  
}

