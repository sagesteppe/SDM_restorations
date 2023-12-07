idig_records <- function(query){
  
  result <- ridigbio::idig_search_records(rq = list(scientificname = query)) |> 
    dplyr::select(scientificname, datecollected, collector, geopoint.lon, geopoint.lat) |> 
    tidyr::drop_na(geopoint.lon, geopoint.lat) |> 
    dplyr::mutate(datecollected = as.Date(datecollected)) |> 
    dplyr::filter(datecollected > '1950-01-01')
  
  if(nrow(result) > 0) {result$status = 'found'} else {
    result <- data.frame(
      'scientificname' = query, 
      'datecollected' = NA, 'collector' = NA, 'geopoint.lon' = NA, 'geopoint.lat' = NA, 
      'status' = 'not found'
    )
  }
  return(result)
}

ONEperCELL <- function(x){
  
  records <- sf::st_transform(x, 5070) |>
    sf::st_buffer( dist = 90)
  
  intersects <- sf::st_intersects(records)
  intersect_sets <- unique(intersects) # identify intersected points
  
  records_list <- vector(mode = 'list', length = length(intersect_sets))
  for (i in seq_along(intersect_sets)){
    records_list[[i]] <- records[unlist((intersect_sets[i])),]
  }
  
  records_list <- bind_rows(
    lapply(records_list, function(x) x[which.max(x$date),])
  )
  
  return(records_list)
  
}

multiplesPERcell <- function(x){
  
  records <- sf::st_transform(x, 5070) |>
    sf::st_buffer( dist = 90)
  
  intersects <- sf::st_intersects(records)
  intersect_sets <- unique(intersects[ lengths(intersects) > 1])
  
  records_list <- vector(mode = 'list', length = length(intersect_sets))
  for (i in seq_along(intersect_sets)){
    records_list[[i]] <- records[unlist((intersect_sets[i])),]
  }
  
  min <- lapply(records_list, function(x) x[which.min(x$date), ])
  max <- lapply(records_list, function(x) x[which.max(x$date), ])
  dates_min <- dplyr::bind_rows(min); dates_max <- dplyr::bind_rows(max)
  subsets <- dates_min$date != dates_max$date
  
  min <- min[subsets] ; max <- max[subsets]
  same_sites <- dplyr::bind_rows(min, max) 
  if(nrow(same_sites) > 1){
  	same_sites <- dplyr::mutate(same_sites, SitePair = rep(1:length(min), times = 2), .before = geometry) |>
      dplyr::arrange(SitePair)} else {
    same_sites <- tibble(taxon = NA, date = NA, SitePair = NA, geometry = NA)}
  
  return(same_sites)
  
}

