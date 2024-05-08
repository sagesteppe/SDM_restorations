idig_records <- function(query){
  
  
  fields2getrecord <- c(
    "data.dwc:scientificName",  "data.dwc:decimalLatitude", 
    "data.dwc:decimalLongitude", "collector", "uuid",
    "data.dwc:year", "data.dwc:month", "data.dwc:day")
  names(fields2getrecord) <- c(
    'scientificname', 'geopoint.lat', 'geopoint.lon', 'collector', 'uuid', 'year', 'month', 'day')
  
  result <- ridigbio::idig_search_records(rq = 
                                            list(scientificname = query), fields = fields2getrecord) |> 
    dplyr::rename(all_of(fields2getrecord)) |>
    dplyr::select(scientificname, year, collector, geopoint.lon, geopoint.lat) |> 
    tidyr::drop_na(geopoint.lon, geopoint.lat, year) |> 
    dplyr::filter(year >= 1950)
  
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

#' @param x sf/tibble/dataframe containing cleaned presences
#' @param bg_abs the proportion of records which should be background absences, i..e come from outside the species known range
#' @param terrestrial an object of terrestrial surfaces
#' @param rand-dist the distance offset between the presences and random absences. 
absence_drawer <- function(x, bg_abs, terrestrial, rand_dist){
  
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
  buffered_presences <- sf::st_buffer(x, dist = rand_dist) |>
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
  
  recs_w_vars <- terra::extract(layers, x)
  
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
  
  
  if(!is.na(removals[1])){data_out <- dplyr::filter(x, ! ID %in% removals)}
  
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
randomForests <- function(train_data, test_data, species, suffix){
  
  if(missing(suffix)){suffix <- '1k'}
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
          file = paste0('../results/rf_models/', species, suffix, '-', gsub(' ', '_', Sys.time()), '.rds'))
  
  vip::vip(rf_model) +
    theme_bw() +
    labs(title = species) + 
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(paste0('../results/vip_plots/', species, suffix, '.png'), device = 'png')
  
  prediction_prob <- predict(rf_model, test_data, type = 'prob') # use this for auc/roc
  result.roc <- pROC::roc(test_data$Occurrence, prediction_prob[,1]) # Draw ROC curve.
  
  png(file = paste0('../results/roc_plots/', species, suffix, '.png'))
  plot(
    result.roc,
    print.thres="best",
    print.thres.best.method="closest.topleft")
  dev.off(); graphics.off()
  
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
    setNames(data.frame((ncol(train_data)-1), 'Independent_Variables'), c('Value', 'Metric')), 
    summary, samples)
  write.csv(summary, paste0('../results/summary/', species, suffix, '.csv'), row.names = F)
  
}


modeller <- function(x, suffix){
  
  if(missing(suffix)){suffix <- '1k'}
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
  TRAIN_DATA <- rfe_var_selector(TRAIN_DATA)
  randomForests(train_data = TRAIN_DATA, test_data = TEST, species, suffix)
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


rfe_var_selector <- function(x){
  
  ctrl <- caret::rfeControl(functions = caret::rfFuncs, method = "repeatedcv", number = 10,
                            repeats = 5, rerank = TRUE, allowParallel = TRUE)
  
  subsets = seq(from = 10, to = (floor((ncol(x)-1)*0.1)*10), by = 5)
  
  cl <- parallel::makeCluster(cores, type='PSOCK')
  doParallel::registerDoParallel(cl)
  rfProfile <- caret::rfe(x[,2:ncol(x)], x$Occurrence,
                          sizes = subsets, rfeControl = ctrl)
  ParallelLogger::stopCluster(cl)
  
  # accept a percent change with less than -1.5% 
  results <- rfProfile[['results']] %>% 
    dplyr::mutate(pct_change = (Accuracy/lead(Accuracy) - 1) * 100)
  
  v_no <- if(any(results$pct_change > 0 & results$Variables < 20)){ # here if performance increases, select the
    A <- results[results$pct_change  > 0 & results$Variables < 20,] # increased value with the fewest number of terms. Do not allow the 
    v_no <- A[which.min(A$Variables), 'Variables'] # 30 vars!! that is essentially where we started.
  } else { # if perfomance only moderately decreases across all models with fewer vars, select that with 
    B <- results[results$pct_change > -1.5,] # the fewest terms, but with less than a 1.5% drop in performance
    v_no <- B[which.min(B$Variables), 'Variables']
  }
  
  rfe_tab <- rfProfile[["variables"]] %>% 
    dplyr::filter(Variables == v_no) %>% 
    dplyr::group_by(var) %>% 
    dplyr::summarise(imp = mean(Overall)) %>% 
    dplyr::arrange(-imp) 
  
  x_sub <- x[ , names(x) %in% c('Occurrence', rfe_tab$var)]
  return(x_sub)
  
}




# this bash script was made to write files to the sdd, and copy to an hdd periodically 

#!/bin/bash
# cd /home/steppe/Documents/suitability_maps
#maps=(*) # save all file names into object
#noFiles="$((${#maps[@]}-1))" # the number of completed files. 
#if [ "$noFiles" -gt 1 ]; then
#maps2move=${maps[@]::${#maps[@]}-1} # subset the completed files
#  for f in $maps2move; do 
#  mv "$f" "/media/steppe/hdd/SDM_restorations/results/suitability_maps/"
#  done
#  else
#    :
#    fi
# we set up a cron job to run every24 hours to move the files
# 0 */3 * * * sdd2hdd.sh




patcheR <- function(x){
  
  # read in file
  r <- rast(x)
  taxon <- gsub('1k.*$', '', basename(x))
  pout <- '../results/patches'
  fsize <- file.info(x)$size/1e9
  
  # mask raster values < 0.8
  r <- terra::mask(r, ifel(r < 0.80, NA, r))
  
  # burn away rivers, and burn patches to hydrologic units ~ populations
  r <- terra::mask(r, terra::crop(hydr_bound, r), inverse = TRUE)
  r <- terra::mask(r, terra::crop(hydr_ftr, r), inverse = TRUE)
  
  # detect patches - this scales to the extent of analysis. #
  # Areas need to be aggregateed for computatations to be run in meaningful amount of
  # time
  
  if(fsize < 0.3){r1 <- terra::aggregate(r, 2, fun = 'mean')} else
    if(fsize < 0.5){r1 <- terra::aggregate(r, 3, fun = 'mean')} else 
      if(fsize < 0.7){r1 <- terra::aggregate(r, 4, fun = 'mean')} else
        if(fsize > 0.7){r1 <- terra::aggregate(r, 5, fun = 'mean')}
  
  pat <- terra::patches(r1,  directions = 4, allowGaps = FALSE)
  
  # remove patches < 5 acres
  sizes <- terra::zonal(cellSize(pat, unit="ha"), pat, sum, as.raster=TRUE)
  size_mask <- terra::ifel(sizes < 4.0470, NA, sizes)
  r1 <- mask(r1, size_mask)
  
  # resample the patches to the original raster resolution
  terra::resample(r1, r, threads = 16, method = 'near',
                  filename = file.path(pout, paste0(taxon, '.tif')))
  
  terra::tmpFiles(current = FALSE, orphan = TRUE, old = TRUE, remove = TRUE)
  
}



#' identify patches which are known to have populations of a target species in them
patchTaggR <- function(x){
  
  r <- terra::rast(x)
  r_taxon <- gsub('_', ' ', gsub('-patchIDS[.]tif', '', basename(x)))
  occs <- dplyr::filter(recs, taxon == r_taxon) |>
    dplyr::mutate(ID = dplyr::row_number()) |>
    terra::vect()
  
  occ_pts <- terra::extract(r, occs, bind = TRUE) |>
    sf::st_as_sf()
  
  # if pts were NA, buffer to 90m, many occupied areas were 'cut' away by the 
  # hydrologic mapping. 
  occ_poly <- dplyr::filter(occ_pts, is.na(lyr.1)) |>
    sf::st_buffer(250) # we cut out 90, but let's go a bit further for corner cases.
  
  occ_poly <- terra::extract(r, occ_poly, weights = TRUE) %>% 
    dplyr::group_by(ID) |> # now make sure we get the patch with the 'most' overlap with the polygon
    tidyr::drop_na(lyr.1) |> # these values actually NA. 
    dplyr::slice_max(weight, with_ties = F) |>
    dplyr::left_join(filter(occ_pts) |>
                       select(-lyr.1), by = join_by('ID')) |>
    dplyr::select(-weight)
  
  occ_pops <- dplyr::bind_rows( # a data frame of all raster patch ID's with populations known from within
    # the last few decades in them. 
    tidyr::drop_na(occ_pts, lyr.1), occ_poly) |>
    dplyr::rename(PatchID = lyr.1)
  
  # write out summary values of the patch matching
  
  fp <- '../results/occupiedPatches'
  data.frame(
    noPres = nrow(occs), 
    PresTopSuitability = nrow(occ_pops), # the number of presences located in top suitability patches
    UniquePatchID = nrow(dplyr::distinct(occ_pops, PatchID)), # the number of uniquely identified patches
    NoPatchesTopSuitability = terra::minmax(r)[2] # the total number of most suitable patches from the RASTER
  ) |>
    write.csv(row.names = F, 
              file.path(fp, paste0(gsub(' ', '_', r_taxon), '-summary', '.csv')))
  
  dplyr::distinct(occ_pops, PatchID, .keep_all = TRUE) |> # this is multiple pts per patch - not necessary for us. 
    sf::st_drop_geometry() |>
    dplyr::select(date, PatchID) |>
    dplyr::arrange(PatchID) |>
    write.csv(row.names = F,
              file.path(fp, paste0(gsub(' ', '_', r_taxon), '-topSuitability', '.csv')))
  
  # for points which were not patched to a most suitable patch, put them into a hydrologic basin
  marginal_pops <- sf::st_as_sf(occs)|>
    dplyr::filter(!ID %in% occ_pops$ID)
  
  # identify the drainage's within the HU12 database
  drain_nos <- terra::extract(hu12, marginal_pops, bind = TRUE) |>
    sf::st_as_sf()
  r_drain <- terra::ifel(hu12 %in% drain_nos$huc12, 1, NA)
  r_drain_sub <- terra::crop(hu12, r_drain, mask = TRUE)
  
  # mask original suitability map to these areas, and filter > 0.5 suitable habitat
  f1 <- list.files('../results/suitability_maps')
  # f1 <- f1[grepl('1k', f1)]# this for the first records
  f1 <- f1[ grepl(pattern = gsub(' ', '_', r_taxon), gsub('1k.*$', '', f1)) ][1]
  
  suitable <- terra::rast(paste0('../results/suitability_maps/', f1))
  r_drain_sub <- terra::resample(r_drain_sub, suitable, threads = 16, method = 'near')
  suitable <- terra::crop(suitable, r_drain_sub, mask = TRUE)
  msk <- terra::ifel(suitable > 0.5 & suitable < 0.8, 1, NA) # let's make sure these areas do not include existing clusters. 
  suitable <- terra::mask(suitable, terra::crop(hydr_bound, suitable), inverse = TRUE)
  terra::mask(suitable, msk, overwrite = TRUE,
              filename = paste0('../results/marginal_habitat/', gsub(' ', '_', r_taxon), '.tif'))
  # this is now a second raster for occupied marginal populations
  
  terra::tmpFiles(current = FALSE, orphan = TRUE, old = TRUE, remove = TRUE)
}


combine_pop_patches <- function(x){
  
  pIDs <- read.csv(paste0('../results/occupiedPatches/', x, '-topSuitability.csv'))
  top_patches <- terra::rast(paste0('../results/fine_patches/', x, '-patchIDS.tif'))
  highest_no <- minmax(top_patches)[2]
  top_patches <- terra::subst(top_patches, pIDs$PatchID,  pIDs$PatchID, others = NA)
  names(top_patches)<- 'patches'
  pts <- st_read(paste0('../data/raw/occurrence/combined_records/', x, '.shp'), quiet = TRUE) |>
    dplyr::filter(Occurrence == 1)
  marg_patch <- terra::rast(paste0('../results/marginal_habitat/', x, '.tif'))
  
  fsize <- file.info(paste0('../results/marginal_habitat/', x, '.tif'))$size/1e9
  if(fsize < 0.025){marg_patch <- terra::patches(marg_patch, directions = 4)} else
    if(fsize > 0.025){marge1 <- terra::patches(terra::aggregate(marg_patch, 2, fun = 'mean'), directions = 4)} else
      if(fsize > 0.04){marge1 <- terra::patches(terra::aggregate(marg_patch, 3, fun = 'mean'), directions = 4)} else
        if(fsize > 0.05){marge1 <- terra::patches(terra::aggregate(marg_patch, 4, fun = 'mean'), directions = 4)} else
          if(fsize > 0.06){marge1 <- terra::patches(terra::aggregate(marg_patch, 5, fun = 'mean'), directions = 4)} else
            if(fsize > 0.07){marge1 <- terra::patches(terra::aggregate(marg_patch, 6, fun = 'mean'), directions = 4)}
  if(exists('marge1')){marg_patch <- terra::resample(marge1, marg_patch, threads = 16, method = 'near')} 
  
  occ_pts <- terra::extract(marg_patch, terra::vect(pts), bind = TRUE) |>
    sf::st_as_sf()
  occ_poly <- dplyr::filter(occ_pts, is.na(patches)) |>
    sf::st_buffer(250) # we cut out 90, but let's go a bit further for corner cases.
  occ_poly <- terra::extract(marg_patch, terra::vect(occ_poly), weights = TRUE) %>% 
    dplyr::group_by(ID) |> # now make sure we get the patch with the 'most' overlap with the polygon
    tidyr::drop_na(patches) |> # these values actually NA. 
    dplyr::slice_max(weight, with_ties = F)
  
  occ_patches <- c(
    dplyr::filter(occ_pts, ! is.na(patches)) |>
      dplyr::pull(patches), 
    occ_poly$patches
  )
  
  if(length(occ_patches) > 0){ 
    
    marg_patch <- terra::subst(marg_patch, occ_patches,  occ_patches, others = NA)
    old_ID <- unique(marg_patch)$patches # gather and relabel patches 
    lkp <- data.frame(
      old_ID = old_ID,
      new_ID = highest_no + 1:length(old_ID)
    )
    marg_patch <- terra::subst(marg_patch, lkp$old_ID, lkp$new_ID)
    final_patches <- terra::cover(top_patches, marg_patch)
    
    terra::writeRaster(final_patches, paste0('../results/allOccupiedPatches/', x, '.tif'))
  } else {
    terra::writeRaster(top_patches, paste0('../results/allOccupiedPatches/', x, '.tif'))
  }
}


landscapR <- function(x, pout, pRout){
  
  taxon <- gsub('[.]tif', '', basename(x))
  x <- terra::rast(file.path('../results/patches', x))
  x <- terra::ifel(x > 0, 1, NA)
  x <- terra::mask(x, terra::crop(hydr_bound, x), inverse = TRUE)
  landscape <- landscapemetrics::get_patches(x, directions = 4)
  calcs <- landscapemetrics::calculate_lsm(x, what = c('lsm_c_enn_mn', 'lsm_c_enn_cv',
                                                       "lsm_p_enn", "lsm_p_cai", "lsm_p_para", "lsm_p_frac"), 
                                           neighbourhood = 4, directions = 4)
  
  p <- '../results/'
  write.csv(
    calcs[calcs$level == 'patch', c('id', 'metric', 'value')],
    file = paste0(p, pout, '/', taxon, '-', 'patch_vals.csv'),  row.names = F
  )
  
  write.csv(
    calcs[calcs$level == 'class', c('metric', 'value')],
    file = paste0(p, pout, '/', taxon, '-', 'class_vals.csv'), row.names = F
  )
  
  writeRaster(landscape[["layer_1"]][["class_1"]], 
              filename = paste0(p, pRout, '/', taxon, '-', 'patchIDS.tif')
  )
  
}




#' identify contiguous neighbors and count them. 
#' @param x 
neigh_type <- function(x){
  
  # assemble the data sets. 
  taxon <- gsub('.shp$', '', basename(x))
  pm <- read.csv(paste0('../results/patch_metrics/', taxon, '-patch_vals.csv')) %>% 
    dplyr::filter(metric == 'cai' & value > 0) %>% 
    dplyr::select(id)
  
  occ <- sf::st_read(paste0('../results/allOccupiedPatchesVector/', taxon, '.shp'), quiet = T) |> 
    dplyr::mutate(Presence = 'Known')
  pred <- sf::st_read(paste0('../results/fine_patchesVector/', taxon, '-patchIDS.shp'), quiet = T) |> 
    dplyr::mutate(Presence = 'Predicted') |> 
    dplyr::filter(lyr_1 %in% pm$id) # filter for core area indexes
  
  pred1 <- pred[ lengths(sf::st_covers(pred, occ)) == 0, ] # make this dataset distinct. 
  pred2 <- sf::st_buffer(pred1, dist = 80) # for up north. 
  pred2 <- sf::st_difference(pred2)
  
  occ <- sf::st_buffer(occ, dist = 80) # for up north. 
  occ <- sf::st_difference(occ)
  
  patches <- dplyr::bind_rows(occ, pred2) %>% 
    sf::st_difference()
  
  # identify contiguous neighbors 
  nblags <- spdep::nblag(
    neighbours = spdep::poly2nb(patches, queen = TRUE), maxlag = 5)
  
  # rows 1: here are presences. 
  known_r <- nrow(patches[patches$Presence=='Known',])
  
  # count the contiguous neighbors
  nblags <- rapply(nblags, function(x){x <- length(as.numeric(x[x <= known_r & x != 0]))}, how = 'replace')
  nblags <- data.table::rbindlist(nblags) %>% 
    t() %>% 
    as.data.frame() %>% 
    setNames(., c('FirstOrder', 'SecondOrder', 'ThirdOrder', 'FourthOrder', 'FifthOrder')) 
  rownames(nblags) <- 1:nrow(nblags)
  
  # count the proximal neighbors
  nbdist <- spdep::dnearneigh(x = st_centroid(patches), d1 = 0, d2 = 5000) # neighbors within 5k
  nbdist <- lapply(nbdist, function(x){x <- length(as.numeric(x[x <= known_r & x != 0]))})
  
  nbdist <- setNames(
    data.frame(as.matrix(nbdist)), 'OccupiedNeighborsIn5k')
  nbdist$OccupiedNeighborsIn5k <- as.numeric(nbdist$OccupiedNeighborsIn5k)
  
  dat <- cbind(nbdist, nblags) %>% 
    tibble::rownames_to_column('ROW') %>% 
    tidyr::pivot_longer(FirstOrder:FifthOrder, names_to = 'OrderConnection', values_to = 'NoConnections') %>% 
    dplyr::filter(ROW %in% 1:known_r | OccupiedNeighborsIn5k > 0 | NoConnections > 0)
  
  ranks <- dplyr::mutate(dat, 
                         Rank = dplyr::case_when(
                           OrderConnection == 'FirstOrder' & NoConnections >= 2 ~ 2, 
                           OrderConnection == 'FirstOrder' & NoConnections == 1 ~ 3, 
                           OrderConnection == 'SecondOrder' & NoConnections >= 3 ~ 4, 
                           OrderConnection == 'SecondOrder' & NoConnections <= 2 ~ 5, 
                           OrderConnection == 'ThirdOrder' & NoConnections >= 4 ~ 6, 
                           OrderConnection == 'ThirdOrder' & NoConnections <= 3 ~ 7, 
                           OccupiedNeighborsIn5k >= 3 ~ 5, 
                           OccupiedNeighborsIn5k <= 2 ~ 6, 
                           .default = as.numeric(7)
                         )) |>
    dplyr::slice_min(Rank, n = 1, by = 'ROW') |>
    dplyr::select(ROW, Rank) |>
    dplyr::mutate(Rank = dplyr::if_else(ROW %in% 1:known_r, 1, Rank))
  
  x1 <- patches |>
    tidyr::unite('PatchIDCombined', c('patches', 'lyr_1'), na.rm = T) |>
    tibble::rownames_to_column('ROW') |>
    dplyr::right_join( ranks, by = 'ROW') |>
    dplyr::select(-ROW, -Presence, PchIDComb = PatchIDCombined)
  
  write.csv(dat, paste0('../results/connections/', taxon, '.csv'), row.names = F)
  sf::st_write(x1, paste0('../results/rankedPatches/', taxon, '.shp'), quiet = T, append = F)
  
}


############################################################################
# write out portions of 

#' Write out all spatial data for crews to a directory
#'
#' @param project_areas split sf data set, with each crews field office(s), and name
#' @param admu administrative boundaries
#' @param intype one of 'crew' or 'senior', defaults to crew
#' @param target_species a dataframe containing possible target species
#' @param crew_id column in project_areas holding the crews identifier info 'e.g. UFO'
#' @param blm_surf sf data set of surface management as multipolygon at most
#' @param fs_surf sf data set of surface USFS land management as multipolygon at most
#' @param fire sf data set of historic fires as multipopylgon at most
#' @param invasive sf dataset of invasive species relative cover estimates
#' @param historic_SOS sf dataset with historic SoS endeavors and relevant info
#' @param roads sf dataset of roads with relevant attritbutes
#' @param seed_transfer sf dataset of seed transfer zones 
#' @param total_change Total Change Intensity Index from the MRLC
#' @param drought6 netcdf of drought dataset from SPEI website, we recommend using 6 and 12 month
#' @param drought12 12 netcdf of drought dataset from SPEI website, we recommend using 6 and 12 month
project_maker_original <- function(x, target_species, intype,
                          admu, 
                          blm_surf, fs_surf, # ownership stuff
                          fire,  invasive,  occurrences, historic_SOS, 
                          total_change, 
                          roads, seed_transfer#, 
                     #     drought6, drought12
                     ){
  
  defaultW <- getOption("warn")
  options(warn = -1)
  if(missing(intype)){intype <- 'crew'}
  
  #### Initiate Project Directories ####
  # create directory to hold all contents
  
  if(intype == 'crew'){
    crew_dir <- paste0('/media/steppe/hdd/2024SOS_CrewGeospatial/', x['Contract'][[1]][1])
  } else {crew_dir <- paste0('/media/steppe/hdd/2024SOS_CrewGeospatial/', x['Senior'][[1]][1])}
  
  ifelse(!dir.exists(file.path('../Crews/')), dir.create(file.path('../Crews/')), FALSE)
  ifelse(!dir.exists(file.path(crew_dir)), dir.create(file.path(crew_dir)), FALSE)
  ifelse(!dir.exists(file.path(crew_dir, 'Geodata')), 
         dir.create(file.path(crew_dir, 'Geodata')), FALSE)
  ifelse(!dir.exists(file.path(crew_dir, 'Data')), 
         dir.create(file.path(crew_dir, 'Data')), FALSE)
  
  # write out BLM and Forest Service
  dir.create( file.path(crew_dir, 'Geodata/Admin'), showWarnings = F) 
  dir.create( file.path(crew_dir, 'Geodata/Admin/Boundaries'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Admin/Allotments'), showWarnings = F)
  dir.create( file.path(crew_dir, 'Geodata/Admin/Surface'), showWarnings = F ) # both BLM and Forest Service go in here
  
  # fire and invasive species data
  dir.create( file.path(crew_dir, 'Geodata/Disturb'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Disturb/Fire'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Disturb/Invasive'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Disturb/TotalChangeIntensity'), showWarnings = F ) 
  
  # target species information
  dir.create( file.path(crew_dir, 'Geodata/Species'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Species/SDM'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Species/SDM-raw'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Species/Occurrences'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Species/Historic_SoS'), showWarnings = F)
  
  # Roads
  dir.create( file.path(crew_dir, 'Geodata/Roads')) 
  # Seed transfer zones
  dir.create( file.path(crew_dir, 'Geodata/STZ')) 
#  # drought 
#  dir.create( file.path(crew_dir, 'Geodata/Drought'))
  
  #### Process geographic data to a mask of the field office ####
  focal_bbox <- dplyr::filter(admu, ADMU_NA %in% stringr::str_to_upper(x$FieldOffice) ) %>%  
    sf::st_union() %>% 
    sf::st_transform(5070)
  
  focal_vect  <- focal_bbox %>%  
    terra::vect()
  focal_bbox <- focal_bbox 
  
  # write out ownership details
  filter(admu, ADMU_NA %in% stringr::str_to_upper(x$FieldOffice)) |>
    sf::st_cast("MULTILINESTRING") |> 
    sf::st_write( dsn = file.path(crew_dir, 'Geodata/Admin/Boundaries', 'Field_Office_Boundaries.shp' ), quiet = T, append = F)
  blm_surf_sub <- sf::st_intersection(blm_surf, focal_bbox) 
  sf::st_write(blm_surf_sub, dsn = file.path(crew_dir, 'Geodata/Admin/Surface', 'BLM_Surface.shp'), quiet = T, append = F)
  sf::st_intersection(allotments, focal_bbox) |> 
    sf::st_write(dsn = file.path(crew_dir, 'Geodata/Admin/Allotments', 'Allotments.shp'), quiet = T, append = F)
  sf::st_intersection(fs_surf, focal_bbox) |>
    sf::st_write(dsn = file.path(crew_dir, 'Geodata/Admin/Surface', 'USFS_Surface.shp'), quiet = T, append = F)
  
  # write out invasive species and fire
  fire <- sf::st_intersection(fire, blm_surf_sub)
  fire |> 
    dplyr::mutate(FIRE_YEAR = as.numeric(FIRE_YEAR)) |>  
    sf::st_write(dsn =  file.path(crew_dir, 'Geodata/Disturb/Fire', 'Fire.shp'), quiet = T, append = F)
  
  terra::crop(invasives, focal_vect, mask = T, threads = T, overwrite = T, filename =
         file.path(crew_dir, 'Geodata/Disturb/Invasive', 'Invasive.tif'))
  
  # write out the change_intensity
  terra::crop(total_change, focal_vect, mask = T, threads = T, overwrite = T, filename = 
         file.path(crew_dir, 'Geodata/Disturb/TotalChangeIntensity', 'TCII.tif'))
  
  # write out assorted data 
  roads[lengths(sf::st_intersects(roads, focal_bbox)) > 0,] |>
    sf::st_write( dsn = file.path(crew_dir, 'Geodata/Roads', 'roads.shp'), quiet = T, append = F)
  sf::st_intersection(seed_transfer, focal_bbox) |>  
    sf::st_write( dsn = file.path(crew_dir, 'Geodata/STZ', 'STZ.shp'), quiet = T, append = F)
  
  # drought
  #  crop(drought6, focal_vect, mask = T, threads = T, filename =
  #         file.path(crew_dir, 'Geodata/Drought', 'drought-6.tif'))
  #  crop(drought12, focal_vect, mask = T, threads = T, filename =
  #         file.path(crew_dir, 'Geodata/Drought', 'drought-12.tif'))
  
  
  # write out species occurrence data
  if(intype == 'crew'){ContractInfo <- x$Contract[1]} else{ContractInfo <- unique(x$Contract)}

  t_spp <- target_species |> 
    dplyr::filter(Requisition %in% ContractInfo) |> 
    dplyr::pull(binomial)
  
  if(intype == 'senior'){t_spp <- unique(unlist(t_spp))}
  
  occurrences_sub <- dplyr::filter(occurrences, taxon %in% gsub('_', ' ', t_spp))
  occurrences_sub <- sf::st_intersection(occurrences_sub, blm_surf_sub)
  occurrences_list <- split(occurrences_sub, f = occurrences_sub$taxon)
  
  occ_writer <- function(x){
    
    binomial <- paste0(gsub(' ', '_', sf::st_drop_geometry(x$taxon[1])), '.shp')
    st_write(x,  dsn = file.path(crew_dir, 'Geodata/Species/Occurrences', binomial),
             quiet = T, append = F)
  } 
  
  lapply(occurrences_list, occ_writer)
  
  ### identify target species and load their patches data set
  patches_p <- '/media/steppe/hdd/SDM_restorations/results/PatchesClippedBLM'
  sdm_files <- file.path(patches_p, list.files(patches_p, pattern = 'shp$'))
  sdm_files <- sdm_files[ grep(paste0(gsub(' ', '_', t_spp), collapse = '.shp|'), sdm_files)]
  
  sdm_writer <- function(x){
    f <- sf::st_read(x, quiet = TRUE)
    f_crop <- sf::st_crop(f, focal_bbox)
    f_crop <- f_crop[ lengths(st_intersects(f_crop, blm_surf_sub)) > 0, ] %>% 
      dplyr::filter(!st_is_empty(.)) %>% 
      sf::st_as_sf() 
    f_crop <- st_intersection(f_crop, blm_surf_sub) %>% 
      select(-Unit_Nm,-Loc_Nm)
    f_crop <- f_crop[sf::st_is(f_crop, c('POLYGON', 'MULTIPOLYGON')),]
    
    binomial <- paste0(gsub(' ', '_', basename(x)))
    if(nrow(f_crop) > 1){
      sf::st_write(f_crop, 
                 file.path(crew_dir, 'Geodata/Species/SDM', binomial), append = F, quiet = TRUE)
    }
  }
  
  lapply(sdm_files, sdm_writer)
  
  ## also write out the raw SDM output
  raw_P <- '/media/steppe/hdd/SDM_restorations/results/suitability_maps'
  sdm_files <- file.path(raw_P, list.files(raw_P, pattern = 'tif$'))
  sdm_files <- sdm_files[grepl('1k-', sdm_files)]
  sdm_files1 <- gsub('1k.*$', '', sdm_files)
  sdm_files <- sdm_files[ grep(paste0(gsub(' ', '_', t_spp), collapse = '|'), basename(sdm_files1))]
  
  sdm_writer2 <- function(x){
    f <- terra::rast(x)
    f_crop <- terra::crop(f, focal_bbox)
    f_crop <- terra::mask(f_crop, blm_surf_sub)
    
    f_crop <- terra::mask(f_crop, ifel(f_crop < 0.70, NA, f_crop)) # if suitability < 60% remove
    # now aggregate rasters to half resolutions
    f_crop <- terra::aggregate(f_crop, fact = 2, fun = 'mean', na.rm = TRUE)
    
    binomial <- gsub('1k-.*$', '.tif', basename(x))
    terra::writeRaster(f_crop, 
                   file.path(crew_dir, 'Geodata/Species/SDM-raw', binomial), overwrite = TRUE)
  }
  
  lapply(sdm_files, sdm_writer2)
  
  # now copy the qgis project template. 
  
  fs::dir_copy('/media/steppe/hdd/2024SOS_CrewGeospatial/qml_styles',
               paste0(crew_dir, '/Geodata/qml_styles'), overwrite = TRUE)
  fs::file_copy('/media/steppe/hdd/2024SOS_CrewGeospatial/HitchPlanner.qgz',
                paste0(crew_dir, '/Geodata/HitchPlanner.qgz'), overwrite = TRUE)
 # write.csv(t_spp, file = file.path(crew_dir, 'Data', 'Target-species.csv'), row.names = F)
 # sub <- sdms[[str_remove(names(sdms), '_[0-9].*$') %in% t]]
  
  # st_write(., dsn = file.path(crew_dir, 'Geodata/Species/Occurrences', 'Occurrences.shp'), quiet = T)
 # st_intersection(historic_SOS, focal_bbox) %>% 
#    st_write(., dsn = file.path(crew_dir, 'Geodata/Species/Historic_SoS', 'Historic_SoS.shp'), quiet = T, append = F)
  
  # write out original species information
#  sdm_fo <- crop(sub, focal_vect, mask = T) # need to make a vect of this... 
#  fnames <- paste0(crew_dir, '/Geodata/Species/SDM/', str_remove(names(sdm_fo),
#                                                                 '_[0-9].*$'), ".tif")
#  writeRaster(sdm_fo, fnames)
  
  options(warn = defaultW)
}




























































































############################################################################
# write out portions of 

#' Write out all spatial data for crews to a directory
#' @param x a project species list
#' @param project_areas split sf data set, with each crews field office(s), and name
#' @param admu administrative boundaries
#' @param intype one of 'crew' or 'senior', defaults to crew
#' @param target_species a dataframe containing possible target species
#' @param crew_id column in project_areas holding the crews identifier info 'e.g. UFO'
#' @param blm_surf sf data set of surface management as multipolygon at most
#' @param fs_surf sf data set of surface USFS land management as multipolygon at most
#' @param fire sf data set of historic fires as multipopylgon at most
#' @param invasive sf dataset of invasive species relative cover estimates
#' @param historic_SOS sf dataset with historic SoS endeavors and relevant info
#' @param roads sf dataset of roads with relevant attritbutes
#' @param seed_transfer sf dataset of seed transfer zones 
#' @param total_change Total Change Intensity Index from the MRLC
#' @param drought6 netcdf of drought dataset from SPEI website, we recommend using 6 and 12 month
#' @param drought12 12 netcdf of drought dataset from SPEI website, we recommend using 6 and 12 month
project_maker <- function(x, target_species, intype,
                          admu, 
                          blm_surf, fs_surf, # ownership stuff
                          fire,  invasive,  occurrences, historic_SOS, 
                          total_change, 
                          roads, seed_transfer#, 
                          #     drought6, drought12
){
  
  defaultW <- getOption("warn")
  options(warn = -1)
  if(missing(intype)){intype <- 'crew'}
  
  #### Initiate Project Directories ####
  # create directory to hold all contents
  
  if(intype == 'crew'){
    crew_dir <- paste0('/media/steppe/hdd/2024SOS_CrewGeospatial/', x['Contract'][[1]][1])
  } else {crew_dir <- paste0('/media/steppe/hdd/2024SOS_CrewGeospatial/', x['Senior'][[1]][1])}
  
  ifelse(!dir.exists(file.path('../Crews/')), dir.create(file.path('../Crews/')), FALSE)
  ifelse(!dir.exists(file.path(crew_dir)), dir.create(file.path(crew_dir)), FALSE)
  ifelse(!dir.exists(file.path(crew_dir, 'Geodata')), 
         dir.create(file.path(crew_dir, 'Geodata')), FALSE)
  ifelse(!dir.exists(file.path(crew_dir, 'Data')), 
         dir.create(file.path(crew_dir, 'Data')), FALSE)
  
  # write out BLM and Forest Service
  dir.create( file.path(crew_dir, 'Geodata/Admin'), showWarnings = F) 
  dir.create( file.path(crew_dir, 'Geodata/Admin/Boundaries'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Admin/Allotments'), showWarnings = F)
  dir.create( file.path(crew_dir, 'Geodata/Admin/Surface'), showWarnings = F ) # both BLM and Forest Service go in here
  
  # fire and invasive species data
  dir.create( file.path(crew_dir, 'Geodata/Disturb'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Disturb/Fire'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Disturb/Invasive'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Disturb/TotalChangeIntensity'), showWarnings = F ) 
  
  # target species information
  dir.create( file.path(crew_dir, 'Geodata/Species'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Species/SDM'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Species/SDM-raw'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Species/Occurrences'), showWarnings = F ) 
  dir.create( file.path(crew_dir, 'Geodata/Species/Historic_SoS'), showWarnings = F)
  
  # Roads
  dir.create( file.path(crew_dir, 'Geodata/Roads')) 
  # Seed transfer zones
  dir.create( file.path(crew_dir, 'Geodata/STZ')) 
  #  # drought 
  #  dir.create( file.path(crew_dir, 'Geodata/Drought'))
  
  #### Process geographic data to a mask of the field office ####
  focal_bbox <- dplyr::filter(admu, ADMU_NA %in% stringr::str_to_upper(x$FieldOffice) ) %>%  
    sf::st_union() %>% 
    sf::st_transform(5070)
  
  focal_vect  <- focal_bbox %>%  
    terra::vect()
  focal_bbox <- focal_bbox 
  
  # write out ownership details
  filter(admu, ADMU_NA %in% stringr::str_to_upper(x$FieldOffice)) |>
    sf::st_cast("MULTILINESTRING") |> 
    sf::st_write(dsn = file.path(crew_dir, 'Geodata/Admin/Boundaries', 'Field_Office_Boundaries.shp' ), 
                  quiet = T, append = F)
  blm_surf_sub <- sf::st_intersection(blm_surf, focal_bbox) 

  
  # write out species occurrence data
  if(intype == 'crew'){ContractInfo <- x$Contract[1]} else{ContractInfo <- unique(x$Contract)}
  
  t_spp <- target_species |> 
    dplyr::filter(Requisition %in% ContractInfo) |> 
    dplyr::pull(binomial)
  
  
  if(intype == 'senior'){t_spp <- unique(unlist(t_spp))}
  
  

  occurrences_sub <- dplyr::filter(occurrences, taxon %in% gsub('_', ' ', t_spp))
  occurrences_list <- split(occurrences_sub, f = occurrences_sub$taxon)
  
  occurrences_list <- lapply(
    X = occurrences_list, FUN = subset2stz, target_stz = seed_zones,
    prov = prov_stz, requisition = x$Contract, focal_bbox = focal_vect, 
    blm_surf_sub = blm_surf_sub)
  
  occurrences_list <- Filter(function(x) nrow(x)[1] > 0, occurrences_list)
  occurrences_list <- Filter(Negate(is.null), occurrences_list)

  occ_writer <- function(x){
    
    binomial <- paste0(gsub(' ', '_', sf::st_drop_geometry(x$taxon[1])), '.shp')
    st_write(x,  dsn = file.path(crew_dir, 'Geodata/Species/Occurrences', binomial),
             quiet = T, append = F)
  } 
  
  lapply(occurrences_list, occ_writer)
  
  ### identify target species and load their patches data set ###
  patches_p <- '/media/steppe/hdd/SDM_restorations/results/PatchesClippedBLM'
  sdm_files <- file.path(patches_p, list.files(patches_p, pattern = 'shp$'))
  sdm_files <- sdm_files[ grep(paste0(gsub(' ', '_', t_spp), collapse = '.shp|'), sdm_files)]
  
  sdm_writer <- function(input, x){
    f <- sf::st_read(input, quiet = TRUE)
    f_crop <- sf::st_crop(f, focal_bbox)
    f_crop <- f_crop[ lengths(st_intersects(f_crop, blm_surf_sub)) > 0, ] %>% 
      dplyr::filter(!st_is_empty(.)) %>% 
      sf::st_as_sf() 
    f_crop <- st_intersection(f_crop, blm_surf_sub) %>% 
      select(-Unit_Nm,-Loc_Nm)
    f_crop <- f_crop[sf::st_is(f_crop, c('POLYGON', 'MULTIPOLYGON')),]
    
    lapply(
      X = occurrences_list, FUN = subset2stz, target_stz = seed_zones,
      prov = prov_stz, requisition = x$Contract, focal_bbox = focal_vect, 
      blm_surf_sub = blm_surf_sub)
    
    binomial <- paste0(gsub(' ', '_', basename(x)))
    if(nrow(f_crop) > 1){
      sf::st_write(f_crop, 
                   file.path(crew_dir, 'Geodata/Species/SDM', binomial), append = F, quiet = TRUE)
    }
  }
  
  lapply(sdm_files, sdm_writer, x = x)
  
  ## also write out the raw SDM output
  raw_P <- '/media/steppe/hdd/SDM_restorations/results/suitability_maps'
  sdm_files <- file.path(raw_P, list.files(raw_P, pattern = 'tif$'))
  sdm_files <- sdm_files[grepl('1k-', sdm_files)]
  sdm_files1 <- gsub('1k.*$', '', sdm_files)
  sdm_files <- sdm_files[ grep(paste0(gsub(' ', '_', t_spp), collapse = '|'), basename(sdm_files1))]
  
  sdm_writer2 <- function(x){
    f <- terra::rast(x)
    f_crop <- terra::crop(f, focal_bbox)
    f_crop <- terra::mask(f_crop, blm_surf_sub)
    
    f_crop <- terra::mask(f_crop, ifel(f_crop < 0.70, NA, f_crop)) # if suitability < 60% remove
    # now aggregate rasters to half resolutions
    f_crop <- terra::aggregate(f_crop, fact = 2, fun = 'mean', na.rm = TRUE)
    
    binomial <- gsub('1k-.*$', '.tif', basename(x))
    terra::writeRaster(f_crop, 
                       file.path(crew_dir, 'Geodata/Species/SDM-raw', binomial), overwrite = TRUE)
  }
  
  lapply(sdm_files, sdm_writer2)
  
  # now copy the qgis project template. 
  
  fs::dir_copy('/media/steppe/hdd/2024SOS_CrewGeospatial/qml_styles',
               paste0(crew_dir, '/Geodata/qml_styles'), overwrite = TRUE)
  fs::file_copy('/media/steppe/hdd/2024SOS_CrewGeospatial/HitchPlanner.qgz',
                paste0(crew_dir, '/Geodata/HitchPlanner.qgz'), overwrite = TRUE)
  # write.csv(t_spp, file = file.path(crew_dir, 'Data', 'Target-species.csv'), row.names = F)
  # sub <- sdms[[str_remove(names(sdms), '_[0-9].*$') %in% t]]
  
  # st_write(., dsn = file.path(crew_dir, 'Geodata/Species/Occurrences', 'Occurrences.shp'), quiet = T)
  # st_intersection(historic_SOS, focal_bbox) %>% 
  #    st_write(., dsn = file.path(crew_dir, 'Geodata/Species/Historic_SoS', 'Historic_SoS.shp'), quiet = T, append = F)
  
  # write out original species information
  #  sdm_fo <- crop(sub, focal_vect, mask = T) # need to make a vect of this... 
  #  fnames <- paste0(crew_dir, '/Geodata/Species/SDM/', str_remove(names(sdm_fo),
  #                                                                 '_[0-9].*$'), ".tif")
  #  writeRaster(sdm_fo, fnames)
  
  options(warn = defaultW)
}


################################################################################
#' reduce the occurrences, known basins, and raw SDM outputs to desired seed transfer zones
#' on a per project basis
#' @param x the input species list
#' @param target_stz a data frame with the Taxons name, Contract, desired STZs, and STZ product related to them. 
#' @param prov A STACK OF PROVISIONAL RASTERS, THIS SHOULD BE CROPPED TO THE AREA OF ANALYSIS IN THE FUNCTION AROUD THIS ONE. 
#' @param focal_bbox focal bounding box for study area, supply as a vect() so can be projected if necessary. 
#' @param blm_surf_sub subsetted to blm ownership. 
subset2stz <- function(x, target_stz, prov, requisition, focal_bbox, blm_surf_sub){
  
  # determine whether we are using empirical or provisional seed transfer zones
  
  #### Part 1 using the provisional STZ's e.g. Bower, Desert SW. ####
  # identify all species which can be processed from the Provisional STZ's
  # we have created a local copy of the stz's as the layers of a raster stack 
  # to speed up this process. 
  prov_taxa <- target_stz[target_stz$stz %in% c('Bower2013', 'DesertSW','GreatBasin_Bower', 'Mojave') &
                            target_stz$Requisition %in% c(requisition),]
  prov <- prov[[which(names(prov) == unique(prov_taxa$stz))]]  # subset the raster stack to the relevant layer
  
  # These species will be processed using their empirical seed zones. 
  need_empirical <- target_stz[! target_stz$stz %in% c('Bower2013', 'DesertSW','GreatBasin_Bower', 'Mojave') &
                            target_stz$Requisition %in% c(requisition),]
  
  # now subset the provisional seed zones to the BLM Office 
  prov_crop <- terra::crop(prov, terra::ext(terra::project(focal_bbox, terra::crs(prov))))
  prov_crop <- terra::mask(prov_crop, 
                           terra::project(terra::vect(blm_surf_sub), terra::crs(prov)))
  
  # We will now import both the raw species distribution model output, and the
  # re-sampled 
  
  # subset to species 
  prov_taxon <- prov_taxa[gsub('_', ' ', prov_taxa$binomial) == x$taxon[1], ]
  if(nrow(prov_taxon) > 0){
  
   # need to determine if raster type is categorical or numeric. to match the 
    # properties of that layer of the spatraster
    if(
      all(
        grepl('[A-z]', unique(prov_taxon$st_zone)))){
         msk_vals <- unique(prov_taxon$st_zone)} else {
          msk_vals <- as.numeric(unique(prov_taxon$st_zone))
         }
  
    msk <- terra::ifel(prov_crop %in% msk_vals, prov_crop, NA) 
    prov_crop <- terra::mask(prov_crop, msk) 
    
    
    
    
    
    
    
    # load the sdms here. 
    
    
    
    
    
    
    
    
    
    
    
    
    
  
    # need to be able to dispatch to both raster and vector data, and both point and
    # polygon data. 
    if(class(x)[1] == 'SpatRaster'){

      # we can simply mask non target seed zones in our stz rasters, and then mask
      # our species raster to that raster.
    
      targeted_st_zones <- terra::mask(prov_crop, msk, inverse = TRUE)
      raw_sdm_masked_to_target_stzs <- terra::mask(x, targeted_st_zones)
      return(raw_sdm_masked_to_target_stzs)
    
    } else if(class(x)[1] == 'sf'){ # determine type of vector data and dispatch. 
    
      if(unique(sf::st_geometry_type(x[1,])) == 'POINT'){
    
        targeted <- terra::mask(prov_crop, msk, inverse = TRUE)
        taxon_in_target_stz <- x[which( # identify point records overlapping target stz. 
          !is.na(
            terra::extract(prov_crop, terra::project(terra::vect(x), crs(prov_crop)), 
                           method = 'simple', ID = F))), ]
        # species by species subset to the target zones
        return(taxon_in_target_stz)
      
      } else { # now essentially intersect the areas of polygons. 
      
        targeted_v <- sf::st_as_sf(
          terra::as.polygons(targeted, dissolve = FALSE, values = FALSE))
      
        drainages_in_target_stz <- st_intersection(x, targeted_v)
      
      }
     } else {
        message('This class not expected by this small function\n 
              we take sf and SpatRaster objects. Small overhaul needed?')
      
       } 
    } else {
      # run empirical seed zone species down here
      
      
      
    }
}











# lapply through the species which are needed. 

empiricalSTZ_writer <- function(spp_needed, seed_zones, pSTZs){
  
  sz <- seed_zones[seed_zones$binomial == spp_needed$binomial & 
                     seed_zones$Requisition == spp_needed$Requisition,]
  if(all(
      grepl('[A-z]', unique(sz$st_zone)))){
    msk_vals <- unique(sz$st_zone)} else {
      msk_vals <- as.numeric(unique(sz$st_zone))
    }
  
  # mask to field office # 
  msk <- terra::ifel(pSTZs %in% msk_vals, pSTZs, NA) 
  msk <- terra::crop(pSTZs, msk, mask = TRUE) 
  
  ############# first write out the occurrence points ################
  occ <- sf::st_read(spp_needed$occurrence_path, quiet = T) |>
    dplyr::mutate(taxon = str_replace(taxon, ' ', '_')) |>
    dplyr::filter(taxon == spp_needed$binomial)
  
  taxon_in_target_stz <- occ[which( # identify point records overlapping target stz. 
    !is.na(
      terra::extract(msk, terra::project(terra::vect(occ), crs(pSTZs)), 
                     method = 'simple', ID = F))), ]
    
  binomial <- gsub(' ', '_', sf::st_drop_geometry(spp_needed$binomial[1]))
  sf::st_write(taxon_in_target_stz, dsn = file.path(
    crew_dir, 'Geodata/Species/Occurrences', paste0(binomial,'.shp')),
             quiet = T, append = F)
  
  ############ second write out the raw SDM ################
  raw_sdm <- terra::rast(spp_needed$path_raw_sdm)
  raw_sdm <- terra::crop(raw_sdm, msk, mask = TRUE)
  
  # if suitability < 70% remove
  raw_sdm <- terra::mask(raw_sdm, ifel(raw_sdm < 0.70, NA, raw_sdm)) 
  # now aggregate rasters to half resolutions
  raw_sdm <- terra::aggregate(raw_sdm, fact = 2, fun = 'mean', na.rm = TRUE)
  raw_sdm <- terra::trim(raw_sdm)
  
  terra::writeRaster(raw_sdm, filename = file.path(
    crew_dir, 'Geodata/Species/SDM-raw', paste0(binomial,'.tif')),
    overwrite = T)
  
  ############ third write out the basins SDM ################
  basins <- terra::vect(spp_needed$path_basins)
  drainages_in_target_stz <- terra::crop(basins, raw_sdm)
  msk <- terra::as.polygons(msk > -Inf)
  drainages_in_target_stz <- sf::st_as_sf(
    terra::mask(drainages_in_target_stz, msk))

  sf::st_write(drainages_in_target_stz, dsn = file.path(
    crew_dir, 'Geodata/Species/SDM', paste0(binomial,'.shp')),
    quiet = T, append = F)

}












































############################################################################
# write out whether a model says a species is flowering or not.
#' @param x a list of data.frames of crews which need the estimates.
#' @param path a path to the raster phenology estimate summaries from `spat_summarize` in sagesteppe::SeedPhenology
#' @param admu 
#' @param target_species 
phen_tabulator <- function(x, path, project_areas, admu, target_species){
  
  dates <- date_maker()
  # subset to field offices
  focal_area <- dplyr::filter(admu, ADMU_NA %in% stringr::str_to_upper(x$FieldOffice) ) |> 
    sf::st_union() |> 
    sf::st_transform(5070) |> 
    vect() 
  
  # determine which species we want estimates for.
  t_spp <- target_species |> 
    dplyr::filter(Requisition == x$Contract[1]) |> 
    dplyr::pull(taxon) 

  # determine which species we have estimates for
  f <- list.files(path, pattern = '.tif$', recursive = T) 
  f <- f[grep('summary_doys', f)] 
  f_name <- gsub('[.]tif', '', basename(f)) 
  f <- file.path(path, na.omit(f[match(t_spp, f_name)])) 
  missing <- data.frame( 
    taxon = setdiff(t_spp, f_name)) 
  
  # gather the simple density estimates. 

  if(length(missing) >0){
    
  path <- '/media/steppe/hdd/SeedPhenology/data/processed/lowest_priority_sheets'
  dat <- sf::st_read(file.path(path, 'lowest_priority_sheets.shp'), quiet = TRUE) |>
    dplyr::filter(scntfcnm %in% gsub('_', ' ', missing$taxon)) |>
    dplyr::group_by(scntfcnm) |>
    tidyr::drop_na(doy) |>
    dplyr::mutate(
      q95 = quantile(doy, 0.95), # trim the last %5 
      q25 = quantile(doy, 0.025)) |> # trim only the first 2.5%
    dplyr::filter(doy > q25 & doy < q95) |>
    dplyr::select(-q25, -q95) |>
    sf::st_transform(sf::st_crs(focal_area))
  
  dat_sub <- sf::st_buffer(sf::st_as_sf(focal_area), 120000) |> # ~75 miles
    sf::st_intersection(dat)
  dat_sub <- split(dat_sub, f = dat_sub$scntfcnm)
  
  simple_summarizer <- function(x){
    if(nrow(x) > 5){d <- density(x$doy)} else {d <- median(x$doy)}
    if(length(d)==1){mid <- d} else{mid <- d$x[which.max(d$y)]}
    
    summy <- data.frame(
      taxon = x$scntfcnm[1], 
      doy = round(c(min(x$doy), mid, max(x$doy)), 0),
      event = c('Initiation', 'Peak', 'Cessation')
    ) |> dplyr::mutate( # peak cannot be the same as cessation or initiation
      doy = dplyr::if_else( 
        dplyr::lead(doy, n = 1) <= (doy - 27), doy-28, doy, missing = doy ), 
      doy = dplyr::if_else( 
        dplyr::lag(doy, n = 1) >= (doy - 27), doy+28, doy, missing = doy ), 
      doy = dplyr::if_else(
        (le = dplyr::lead(doy, n = 1)) - (la =dplyr::lag(doy, n = 1)) <= 57, ((le+la)/2), doy, missing = doy)
    )
    return(summy)
  }
  
  if(length(dat_sub)>0){
    s_sum <- dplyr::bind_rows(lapply(dat_sub, simple_summarizer))
  
    s_sum <- dplyr::left_join(s_sum, dates[[1]], by = 'doy') |>
      dplyr::left_join(
        dplyr::select(dates[[2]], -doy), by = c('reports_to' = 'reporting_event'))

  } else {summy <- setNames(data.frame(
    matrix(ncol = 3, nrow = 0)), c('Initiation', 'Peak', 'Cessation'))}
  
  } else {summy <- setNames(data.frame(
    matrix(ncol = 3, nrow = 0)), c('Initiation', 'Peak', 'Cessation'))}
  
  # basically lapply this through the taxa
  summarizer <- function(x){
    focal_r <- terra::rast(x)
    taxon <- gsub('[.]tif', '', basename(x))
    out_v <- terra::extract(focal_r, focal_area, method = 'simple', fun = table, na.rm = T, ID = FALSE)
    pos <- unlist(lapply(out_v, function(x){x1 <- colnames(x)[which.max(x)]})) # most common doy for event
    summy <- data.frame(
      taxon = taxon,
      doy = as.numeric(pos), 
      event = c('Initiation', 'Peak', 'Cessation'))
    
  } 
  
  test <- dplyr::bind_rows(lapply(f, summarizer))
  
  dates <- date_maker()
  # now prepare the date in a format for reporting.
  tdf <- test |>
    dplyr::arrange(taxon) |>
    dplyr::distinct() |>
    dplyr::group_by(taxon) |>
    dplyr::mutate( # peak cannot be the same as cessation or initiation
      doy = dplyr::if_else( 
        dplyr::lead(doy, n = 1) <= (doy - 27), doy-28, doy, missing = doy ), 
      doy = dplyr::if_else( 
        dplyr::lag(doy, n = 1) >= (doy - 27), doy+28, doy, missing = doy ), 
      doy = dplyr::if_else(
        (le = dplyr::lead(doy, n = 1)) - (la =dplyr::lag(doy, n = 1)) <= 57, ((le+la)/2), doy, missing = doy)
      )

  tdf <- dplyr::left_join(tdf, dates[[1]], by = 'doy') |>
    dplyr::left_join(
      dplyr::select(dates[[2]], -doy), by = c('reports_to' = 'reporting_event')) |>
    dplyr::group_by(taxon) 
  
  tdf_peak <- filter(tdf, event == 'Peak') |>
    dplyr::select(taxon, doy) |>
    dplyr::slice(rep(1:dplyr::n(), each = 3)) |>
    dplyr::mutate(
      event = c('Early.Peak', 'Peak', 'Late.Peak'),
      doy = case_when(
        event == 'Early.Peak' ~ doy - 15, # one extra day to ensure we stretch outside
        event == 'Late.Peak' ~ doy + 15, # the current range. 
        .default = doy)
    ) |>
    dplyr::left_join( dates[[1]], by = 'doy') |>
    dplyr::left_join(
      dplyr::select(dates[[2]], -doy), by = c('reports_to' = 'reporting_event')) 
  
  if(exists('s_sum')){
    tdf1 <- dplyr::bind_rows(
      dplyr::filter(tdf, event != 'Peak'),
      s_sum,
      tdf_peak
    )
  } else {
    tdf1 <- dplyr::bind_rows(
      dplyr::filter(tdf, event != 'Peak'),
      tdf_peak
    )
  }
  
  tdf1 <- tdf1 |>
    dplyr::arrange(taxon, doy, event) |>
    dplyr::group_by(taxon, reports_to) |>
    dplyr::mutate(multiples = dplyr::n()) |>
    dplyr::filter(! ( multiples>=2 & str_detect(event, 'Late[.]|Early[.]'))) |>
    dplyr::mutate(multiples = dplyr::n()) |>
    dplyr::filter(! ( multiples>=2 & str_detect(event, 'Peak'))) |>
    dplyr::select(-multiples) |>
    dplyr::arrange(doy) |>
    tidyr::pivot_wider(names_from = reports_to, values_from = event, id_cols = 'taxon',
                       values_fill = NA)
 
  report_pts_rep <- sort(as.numeric(colnames(tdf1)[grep('[0-9]', colnames(tdf1))]))
  miss_wk <- setdiff(min(report_pts_rep):max(report_pts_rep), report_pts_rep) # see if some 

  # weeks are missing by default. 
  if(length(miss_wk)>0){
    
    empty_cols <- setNames(
      data.frame(
        matrix(
          NA, nrow = nrow(tdf1), ncol = length(miss_wk))), 
      miss_wk)
    
    tdf2 <- cbind(tdf1, empty_cols)
    rm(empty_cols)
  } else {tdf2 <- tdf1}
  
  
  tdf2 <- tdf2[, c(
    colnames(tdf2)[grep('[A-z]', colnames(tdf2))], 
    sort(as.numeric(colnames(tdf2)[ grep('[0-9]', colnames(tdf2)) ]))) ]
  
  # now update the column names
  colnames(tdf2)[ grep('[0-9]', colnames(tdf2))] <-  paste0(
    dates[[2]]$month, '-', dates[[2]]$day)[
      as.numeric(colnames(tdf2)[ grep('[0-9]', colnames(tdf2))])
    ]
  
  tdf2 <- cbind(tdf2, setNames(
    data.frame(
      matrix(NA, 
           nrow = nrow(tdf2), 
           ncol = 1)), 'Reliability')
  ) |> 
    dplyr::relocate(Reliability, .after = taxon)
  
  positions <- apply(X = tdf2, MARGIN = 1, function(x)as.numeric(grep("^Peak", x)))
  if(any(is.na(as.numeric(positions)))) {
    positions[which(lapply(positions, length) == 0)] <- median(as.numeric(positions), na.rm = T)
  }
  
  positions <- sort(as.numeric(positions), index.return = T)
  tdf2 <- tdf2[positions$ix,]

  # write out to this location
  crew_dir <- paste0('/media/steppe/hdd/2024SOS_CrewGeospatial/', x['Contract'][[1]][1])
  ifelse(!dir.exists(file.path(crew_dir, 'Data', 'Phenology')), 
         dir.create(file.path(crew_dir, 'Data', 'Phenology')), FALSE)
  crew_dir <- paste0(crew_dir, '/Data/Phenology')
  
  write.csv(tdf2, paste0(crew_dir, '/Phenology_Estimates.csv'), na = "", row.names = F)


}

#' create reporting event dates for helping crews assess phenology
date_maker <- function(){
  
  # subset to the period which the crew is active until-from
  #' rolling fill from using nearest known value. from someone on stackoverflow
  f1 <- function(dat) {
    N <- length(dat)
    na.pos <- which(is.na(dat))
    if (length(na.pos) %in% c(0, N)) {
      return(dat)
    }
    non.na.pos <- which(!is.na(dat))
    intervals  <- findInterval(na.pos, non.na.pos,
                               all.inside = TRUE)
    left.pos   <- non.na.pos[pmax(1, intervals)]
    right.pos  <- non.na.pos[pmin(N, intervals+1)]
    left.dist  <- na.pos - left.pos
    right.dist <- right.pos - na.pos
    
    dat[na.pos] <- ifelse(left.dist <= right.dist,
                          dat[left.pos], dat[right.pos])
    return(dat)
  }
  
  dates <- as.Date(0:365, origin = '2024-01-01')
  dates <- data.frame(
    month = as.numeric(gsub('2024-|-[0-9]{2}', '', dates)),
    day = as.numeric(gsub('2024-[0-9]{2}-', '', dates)),
    doy = 0:365,
    dates
  )
  
  reporting_dates <- dates |>
    dplyr::filter(day %in% c(1,15)) |>
    dplyr::mutate(reporting_event = 1:dplyr::n())
  
  dates <- dplyr::left_join(dates, reporting_dates) |>
    dplyr::mutate(reports_to = f1(reporting_event)) |>
    dplyr::select(doy, reports_to)
  
  reporting_dates <- reporting_dates |>
    dplyr::mutate(
      month = format(as.Date(doy, origin = '2024-01-01'), format = '%B'))
  
  return(list(dates, reporting_dates))
}
