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
