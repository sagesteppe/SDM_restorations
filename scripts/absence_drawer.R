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
  reg_pts <- sf::st_sample(box_hull, size = round(nrow(x) * bg_abs, 0), type = 'regular') |>
    sf::st_as_sf() |>
    dplyr::mutate(Occurrence = 0, PtType = 'Background Absence') |>
    dplyr::rename(geometry = x)
  
  rm(bbox, observed_hull)
  
  # identify AIM plots and generate absences. 
  
  blm_surf_area <- as.numeric(sf::st_area(blm_surf))
  prcnt_blm <- (blm_surf_area / as.numeric(sf::st_area(buffer_hull))) 
  aim_recs <- nrow(x) * (1 - bg_abs) * prcnt_blm
  
  AIM_points <- AIM_points[lengths(st_intersects(AIM_points, buffer_hull)) > 0,]
  AIM_points <- AIM_points |> 
    filter(! PrimaryKey %in% x$PrimaryKey ) # these are presences
  
  if(aim_recs >= nrow(AIM_points)){aim_recs <- nrow(AIM_points)} # in case species v. abundant. 
  
  AIM_absence <- sample_n(AIM_points, size = round(aim_recs * 1.1, 0), replace = FALSE) |> 
    dplyr::mutate(Occurrence = 0, PtType = 'AIM Absence') |>
    dplyr::select(Occurrence, PtType, geometry) |>
    sf::st_transform(sf::st_crs(x))

  # generate Random absences
  
  rand_recs <- nrow(x) * ((1 - bg_abs) * (1 - prcnt_blm))
  buffered_presences <- sf::st_buffer(x, dist = 10000) |>
    sf::st_union() |>
    sf::st_transform(sf::st_crs(blm_surf))
  
  buffer_hull <- sf::st_transform(buffer_hull, sf::st_crs(blm_surf)) |>
    sf::st_make_valid()
  non_blm_area <- sf::st_difference(buffer_hull, blm_surf) 
  non_blm_area <- sf::st_difference(non_blm_area, buffered_presences)
  
  random_samples <- sf::st_sample(non_blm_area, rand_recs, 
                                  size = round(rand_recs * 1.1, 0),
                                  type = 'random') |> 
    sf::st_as_sf() |>
    dplyr::mutate(Occurrence = 0, PtType = 'Random Absence') |>
    sf::st_transform(sf::st_crs(x)) |>
    dplyr::rename(geometry = x)
  
  absences <- bind_rows(reg_pts, AIM_absence, random_samples)
  
  return(absences)
  
}


