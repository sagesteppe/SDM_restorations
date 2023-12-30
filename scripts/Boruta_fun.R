cores <- parallel::detectCores() # run this outside fn to avoid repeat. 

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
  
  x_sub <- x[ , !names(x) %in% removals]# remove
  x_sub <- x_sub[, names(x_sub) %in% c('Occurrence', important_vars)] # keep informative
  
  return(x_sub)
}



