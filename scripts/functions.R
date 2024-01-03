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
  ) |>
    sf::st_centroid()
  
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
      dplyr::arrange(SitePair)|>
  	  sf::st_centroid() } else {
    same_sites <- tibble(taxon = NA, date = NA, SitePair = NA, geometry = NA)}
  
  return(same_sites)
  
}

absence_drawer <- function(x, bg_abs, terrestrial){
  
  
  # create the observed hull, the buffered hull for random sample selection, 
  # and a hull to form the basis for a bounding box to generate regular absences in
  
  observed_hull <- sf::st_union(x) |> # this is the true hull around the observations
    sf::st_convex_hull() |>
    sf::st_sfc() 
  
  obs_area <- as.numeric(st_area(observed_hull)) * 0.15 / 2 # this hull should go at most 100 miles around the inner hull
  obs_area <- if(obs_area > 160934){obs_area = 160934} # cap area with 100 mi buffer. 
  
  buffer_hull <- st_buffer(observed_hull, obs_area) 
  
  buf_area <- as.numeric(sf::st_area(buffer_hull) *0.1) / 2 # this final hull will be used to create a box around the buffer hull
  obs_area <- if(buf_area > 40233){buf_area = 40233} # cap area with 25 mi buffer.
  
  box_hull <- sf::st_buffer(buffer_hull, dist = obs_area)
  bbox <- sf::st_bbox(box_hull) |>
    sf::st_as_sfc() 
  
  box_hull <- sf::st_difference(bbox, buffer_hull)
  
  # remove large oceans, or lakes, or what have you. 
  observed_hull <- sf::st_intersection(observed_hull, terrestrial) |> sf::st_union()
  buffer_hull <- sf::st_intersection(buffer_hull, terrestrial) |> sf::st_union()
  box_hull <- sf::st_intersection(box_hull, terrestrial) |> sf::st_union()
  bbox <- sf::st_intersection(bbox, terrestrial) |> sf::st_union()
  
  # perform the outline background absence sample
  reg_pts <- sf::st_sample(box_hull, size = round(nrow(x) * bg_abs, 0), type = 'regular') |>
    sf::st_as_sf() |>
    dplyr::mutate(Occurrence = 0, PtType = 'Background Absence') |>
    dplyr::rename(geometry = x)
  
  rm(bbox, observed_hull)
  
  # identify AIM plots and generate absences. 
  
  blm_surf <- sf::st_intersection(
    sf::st_transform(buffer_hull, sf::st_crs(blm_surf)), blm_surf)
  blm_surf_area <- as.numeric(sf::st_area(blm_surf))
  prcnt_blm <- (blm_surf_area / as.numeric(sf::st_area(buffer_hull))) 
  aim_recs <- round( nrow(x) * (1 - bg_abs) * prcnt_blm, 0)
  
  AIM_points <- AIM_points[lengths(st_intersects(AIM_points, buffer_hull)) > 0,]
  AIM_points <- AIM_points |> 
    filter(! PrimaryKey %in% x$PrmryKy ) # these are presences
  
  if(aim_recs >= nrow(AIM_points)){aim_recs <- nrow(AIM_points)}
  
  AIM_absence <- sample_n(AIM_points, size = aim_recs, replace = FALSE) |> 
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
  
  random_samples <- sf::st_sample(non_blm_area, size = round(rand_recs * 1.25, 0),
                                  type = 'random') |> 
    sf::st_as_sf() |>
    dplyr::mutate(Occurrence = 0, PtType = 'Random Absence') |>
    sf::st_transform(sf::st_crs(x)) |>
    dplyr::rename(geometry = x)
  
  absences <- dplyr::bind_rows(reg_pts, AIM_absence, random_samples)|>
    dplyr::mutate(Taxon = unique(x$taxon), .before = 'geometry') 
  
  print(unique(x$taxon))
  return(absences)
  
}





quantile_flagger <- function(x, quant, flag_threshold){
  
  if(missing(quant)){quant <- 0.025}
  if(missing(flag_threshold)){flag_threshold <- 0.2}
  
  # perform calculations# 
  quant_horse <- function(x, quant){
    lwr <- quantile(x, quant, na.rm = TRUE)
    upr <- quantile(x, 1 - quant, na.rm = TRUE)
    
    lwr_flag <- which(x <= lwr)
    upr_flag <- which(x >= upr)
    flags <- c(lwr_flag, upr_flag)
    
    return(list(flags))
  }
  x_ng <- sf::st_drop_geometry(x)
  flagged_vals <- apply(x_ng[, grepl(pattern = 'bio*', colnames(x_ng))], MARGIN = 2, 
                        FUN = quant_horse, quant)
  
  nvars <- ncol(x[, grepl(pattern = 'bio*', colnames(x))])
  colsReq2flag <- round(flag_threshold * nvars, 0)
  if(colsReq2flag < 2){colsReq2flag <- 2}
  
  counts <- table(unlist(flagged_vals)) >= colsReq2flag
  flagged <- as.numeric(names(counts[counts == TRUE]))
  
  x['EnvFlag'] <- FALSE
  x[flagged,'EnvFlag'] <- TRUE
  
  # now flag geographically farflung records
  
  dist <- as.numeric(
    sf::st_distance(x, x[sf::st_nearest_feature(x),], by_element = TRUE)
  )
  
  ub_geo <- quantile(dist, (1 - quant) ) 
  x_out <- dplyr::mutate(x, distance = dist, 
                         DistFlag = dplyr::if_else(dist > ub_geo, T, F),
                         OutlierFlag = dplyr::case_when(
                           DistFlag == T & EnvFlag == T ~ 'Dist-Env',
                           DistFlag == T & EnvFlag == F ~ 'Dist',
                           DistFlag == F & EnvFlag == T ~ 'Env',
                           .default = NA,
                         ),
                         .before = geometry) |>
    dplyr::select(-EnvFlag, -DistFlag)
  
  return(x_out)
  
}


colon_balloon <- function(x){
  
  if(stringr::str_detect(x, ':')){
    
    suppressWarnings({
      ranges <- as.vector(str_extract_all(x, '[:digit:]{1,}:[:digit:]{1,}'))
      lower <- stringr::str_extract_all(ranges, '[:digit:]{1,}(?=:)') # lower end of the sequence
      upper <- stringr::str_extract_all(ranges, '(?<=:)[:digit:]{1,}') # upper end of the sequence
      lower <- unlist(lapply(lower, as.numeric)) ; upper <- unlist(lapply(upper, as.numeric))
      numeric_ranges <- as.numeric(unlist( mapply(seq, from = lower, to = upper)))
    })
    
    return(numeric_ranges)
  } else {return(NA)}
  
}






# to determine whether the pseudo-random absences are in areas which are possibly suitable
# we will use linear discrimination analysis (LDA). LDA uses a label, e.g. presence or absence
# as a response which can be predicted by independent variables which are linear combinations 
# generated by reducing a higher dimensional data set to two axis. LDA will be trained on a data set 
# of AIM Plots with true absences, geographic absences, and a subset of the Psuedo-absences. 
# Only up to 750:750 presences:absences per species will be used to train a model. The maximum
# imbalance between presences to absences will be 750:250, or 3:1. 

lda_PA_dropper <- function(x, path, col_names){
  
  recs_w_vars <- terra::extract(bio, x)
  
  recs_w_vars_scaled <- data.frame(apply(recs_w_vars[,2:ncol(recs_w_vars)], MARGIN = 2, FUN = scales::rescale, to = c(0,1)))
  recs_w_vars_scaled$ID <- x$ID
  recs_w_vars_scaled<- tidyr::drop_na(recs_w_vars_scaled)
  x1 <- dplyr::right_join(sf::st_drop_geometry(x), recs_w_vars_scaled, by = 'ID') 
  
  
  pres <- x1[x1$Occurrence == 1, ]
  real_abs <- x1[x1$PtType == 'AIM Absence', ]
  geo_abs <- x1[x1$PtType == 'Background Absence', ]
  model_abs <- dplyr::bind_rows(real_abs, geo_abs)
  
  if(nrow(pres) > 750){pres <- pres[sample(1:nrow(pres), size = 750, replace = FALSE),]}
  if(nrow(model_abs) > 750){model_abs <- model_abs[sample(1:nrow(model_abs), size = 750, replace = FALSE),]}
  
  training <- dplyr::bind_rows(pres, model_abs) |>
    sf::st_drop_geometry()
  training <- training[, c('Occurrence', 'ID', col_names)]
  
  testing <- dplyr::filter(x1, ! ID %in% training$ID) |>
    sf::st_drop_geometry()
  
  LDA_model <- MASS::lda(Occurrence ~ ., training) # generate model
  
  # Here the accuracy of the trained modeled is evaluated
  LDA_train <- predict(LDA_model, training)$class
  tab <- table(Predicted = LDA_train, Actual = training$Occurrence)
  ACC_train <- sum(diag(tab))/sum(tab) # accuracy
  
  # Here the model is applied to unclassified data
  LDA_test <- predict(LDA_model, testing)$class
  testing$Predicted <- predict(LDA_model, testing)$class
  tab1 <- table(Predicted = LDA_test, Actual = testing$Occurrence)
  ACC_test <- sum(diag(tab1))/sum(tab1) # accuracy. 
  
  n_pres <- nrow(x1[x1$Occurrence == 1,])
  n_abs <- nrow(x1[x1$Occurrence == 0,])
  conflicted <- dplyr::filter(testing, Predicted != Occurrence) 
  max_rcs_2_remove <- n_abs - n_pres
  
  if(max_rcs_2_remove < 1 ){
    max_rcs_2_remove <- round(n_pres * 0.1, 0)} # if absences are fewer than presence, allow a 0.9:1 ratio
  
  if(nrow(conflicted) == 0){ # identify the records we want to remove from the test data. 
    removals <- NA
  } else if (nrow(conflicted) > max_rcs_2_remove) {
    removals <- conflicted[sample(1:nrow(conflicted), size = max_rcs_2_remove, replace = FALSE), 'ID']|> dplyr::pull(ID)
  } else {
    removals <- conflicted$ID
  }
  
  if(!is.na(removals)){data_out <- dplyr::filter(x, ! ID %in% removals)}
  
  metrics <- data.frame(
    Train_Presences = nrow(pres), 
    Train_Absences = nrow(model_abs),
    Pres_in = n_pres, 
    Abs_in = n_abs,
    Abs_out = (n_abs - length(removals)),
    Max_Recs_Removal = max_rcs_2_remove,
    ACC_train = round(ACC_train, 5), 
    ACC_test = round(ACC_test, 5),
    Flagged = nrow(conflicted), 
    Removed = length(removals)
  )
  
  species <- x1$taxon[1]
  message(length(removals), ' records dropped from ', species)
  species <- gsub(' ', '_', species)
  write.csv(metrics, file.path(path, paste0(species, '.csv')), row.names = FALSE)
  
  sf::st_write(data_out, file.path('../data/raw/occurrence/combined_records', paste0(species, '.shp')), quiet = TRUE, append = FALSE)
  
}



Boruta_var_selector <- function(x){
  
  BorutaRes <- Boruta::Boruta(Occurrence ~ ., data = x, num.threads = cores, doTrace = 1)
  importance <- Boruta::attStats(BorutaRes)
  rn <- rownames(importance)
  important_vars <- Boruta::getSelectedAttributes(BorutaRes, withTentative = F)
  
  
  ### Remove highly autocorrelated variables ###
  co_pairs <- c('cfvo', 'phh2o', 'gdd', 'ngd', 'gdgfgd')
  pair_removals <- unlist(lapply(co_pairs, function(x){ rn[min(grep(x, rn))]}))
  
  triplets <- c('sand', 'silt', 'clay')
  triplet_removals <- unlist(lapply(triplets, function(x){ rn[grep(x, rn)[1:2]]}))
  
  # reduce two the maximum of three soil textures
  sand_imp <- mean(importance [ grep('sand', rn), 'meanImp'])
  silt_imp <- mean(importance [ grep('silt', rn), 'meanImp'])
  clay_imp <- mean(importance [ grep('clay', rn), 'meanImp'])
  
  texture2remove <- c('sand', 'silt', 'clay')[which.min(c(sand_imp, silt_imp, clay_imp))]
  textures2remove <- rn[ grep(texture2remove, rn)]
  
  removals <- unique(c(pair_removals, triplet_removals, textures2remove))
  
  x_sub <- x[ , !names(x) %in% removals]
  return(x_sub)
}



#' perform random forest modelling and save outputs to locations
#' test_data output from Boruta_var_selector
randomForests <- function(train_data, test_data, species){
  
  # identify the appropriate mtry
  opt_mtry <- data.frame(
    randomForest::tuneRF(train_data[,-1], train_data[,1], stepFactor=1.5))
  opt_mtry <- opt_mtry[ which.min(opt_mtry$OOBError), 'mtry'] # here we subset the lowest OOB. 
  
  # model and predict
  rf_model <- randomForest::randomForest(Occurrence ~ ., data = train_data, mtry = opt_mtry)
  prediction <- predict(rf_model, test_data)
  cmRestrat <- caret::confusionMatrix(prediction, test_data$Occurrence)
  
  # plot and save results
  species <- gsub(' ', '_', species)
  saveRDS(rf_model,
          file = paste0('../results/rf_models/', species, '-', gsub(' ', '_', Sys.time()), '.rds'))
  
  vip::vip(rf_model) +
    theme_bw() +
    labs(title = species) + 
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0('../results/vip_plots/', species, '.png'), device = 'png')
  
  prediction_prob <- predict(rf_model, test_data, type = 'prob') # use this for auc/roc
  result.roc <- pROC::roc(test_data$Occurrence, prediction_prob[,1]) # Draw ROC curve.
  
  png(file = paste0('../results/roc_plots/', species, '.png'))
  plot(
    result.roc,
    print.thres="best",
    print.thres.best.method="closest.topleft")
  dev.off()
  
  summary <- setNames(
    rbind(stack(cmRestrat$byClass), stack(cmRestrat$overall)), 
    c('Value', 'Metric')
  ) 
  samples <- data.frame(
    Value = c(
      nrow(train_data[train_data$Occurrence==0,]), 
      nrow(train_data[train_data$Occurrence==1,]), 
      nrow(test_data[test_data$Occurrence==0,]), 
      nrow(test_data[test_data$Occurrence==1,]) 
    ),
    Metric = c('nTrainPres', 'nTrainAbs', 'nTestPres', 'nTestAbs')
  )
  summary <- rbind(
    setNames(data.frame(as.numeric(result.roc$auc), 'AUC'), c('Value', 'Metric')), 
    summary, samples)
  write.csv(summary, paste0('../results/summary/', species, '.csv'), row.names = F)
  
}

modeller <- function(x){
  
  species <- sf::st_drop_geometry(x) %>% 
    dplyr::pull(taxon)
  species <- species[1]
  x1 <- terra::extract(layers, x)
  
  x_dat <- cbind(x$Occurrence, x1) %>% 
    dplyr::mutate(Occurrence = as.factor(`x$Occurrence`)) %>% 
    dplyr::select(-ID, -`x$Occurrence`) %>% 
    dplyr::relocate(Occurrence, .before = 1) %>% 
    dplyr::filter(if_all(2:ncol(.), ~ !is.na(.)))
  
  trainIndex <- caret::createDataPartition(x_dat$Occurrence, p = .8, 
                                           list = FALSE, times = 1)
  TRAIN <- x_dat[ trainIndex,]
  TEST  <- x_dat[-trainIndex,]
  
  TRAIN_DATA <- Boruta_var_selector(TRAIN)
  randomForests(train_data = TRAIN_DATA, test_data = TEST, species)
  message(species, ' complete')
  
}


predict_wrapper <- function(x){
  
  sdm <- readRDS(x$path2model)
  occ_path <- file.path(occpath, x$occurrence_data)
  
  # read in occurrences and subset prediction to 100 miles from border
  bound <- sf::st_read(occ_path, quiet = TRUE) |>
    dplyr::filter(Occurrence == 1) |>
    sf::st_union() |>
    sf::st_convex_hull() |>
    sf::st_buffer(160934) |>
    sf::st_bbox()
  
  pout <- '../results/suitability_maps'
  layer_subset <- terra::crop(layers, bound)
  terra::predict(
    type = 'prob', cores = 16,
    layer_subset, sdm, cpkgs="randomForest",
    filename = file.path(pout, paste0(x$species, '.tif')),
    wopt = c(names = 'predicted_suitability'))
  
}


term_grab <- function(x){
  
  pieces <- unlist(
    strsplit(
      as.character(x[['terms']][[3]]),
      " ")
  )
  
  terms <- pieces[grepl('[A-z]', pieces)]
  return(terms)
}
