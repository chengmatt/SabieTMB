# Purpose: Function to run francis reweighting

francis_rwgt <- function(data, model) {
  
# Fishery Ages ------------------------------------------------------------

  # Get weights with data 
  new_fish_age_wts <- matrix(nrow = nrow(data$Wt_FishAgeComps), ncol = ncol(data$Wt_FishAgeComps)) # matrix to store new weights
  fish_age_indices <- which(!is.na(data$Wt_FishAgeComps), arr.ind = TRUE)
  
  for(i in 1:nrow(fish_age_indices)) { # loop through data sources for fishery ages with data
    # Filter to years and fleets with data
    tmp_obs <- data$ObsFishAgeComps[,,fish_age_indices[i,1],fish_age_indices[i,2]]
    tmp_obs <- na.omit(tmp_obs, drop = FALSE)
    tmp_iss_obs <- data$ISS_FishAgeComps[,fish_age_indices[i,1],fish_age_indices[i,2], drop = FALSE]
    tmp_iss_obs <- tmp_iss_obs[rownames(tmp_iss_obs) %in% rownames(tmp_obs)]
    
    # Get expected values from model
    rownames(model$rep$CAA) <- data$yrs # input row names
    tmp_exp <- model$rep$CAA[rownames(model$rep$CAA) %in% rownames(tmp_obs),,fish_age_indices[i,1],fish_age_indices[i,2]] # index to correct expected dimensions
    tmp_exp <- t(apply(X = tmp_exp, MARGIN = 1, FUN = function(x) x / sum(x))) # normalize
    data_yrs <- nrow(tmp_obs)
    
    # Set up reweighting vectors
    exp_bar <- vector() # mean expected
    obs_bar <- vector() # mean observed
    v_y <- vector() # variance
    w_denom <- vector() # weight factor in denominator
    
    for(y in 1:data_yrs) {
      exp_bar[y] <- sum(data$ages * tmp_exp[y,]) # get mean pred comps
      obs_bar[y] <- sum(data$ages * tmp_obs[y,]) # get mean obs comps
      v_y[y]<-sum(data$ages^2*tmp_exp[y,])-exp_bar[y]^2 # get variance
      w_denom[y]<-(obs_bar[y]-exp_bar[y])/sqrt(v_y[y]/tmp_iss_obs[y]) # get weights
    } # end y loop
    
    new_fish_age_wts[fish_age_indices[i,1],fish_age_indices[i,2]] <- 1 / var(w_denom) # return weights
    
  } # end i loop
  
  # Fishery Lengths ------------------------------------------------------------
  
  # Get weights with data 
  new_fish_len_wts <- matrix(nrow = nrow(data$Wt_FishLenComps), ncol = ncol(data$Wt_FishLenComps)) # matrix to store new weights
  fish_len_indices <- which(!is.na(data$Wt_FishLenComps), arr.ind = TRUE)
  
  for(i in 1:nrow(fish_len_indices)) { # loop through data sources for fishery ages with data
    # Filter to years and fleets with data
    tmp_obs <- data$ObsFishLenComps[,,fish_len_indices[i,1],fish_len_indices[i,2]]
    tmp_obs <- na.omit(tmp_obs, drop = FALSE)
    tmp_iss_obs <- data$ISS_FishLenComps[,fish_len_indices[i,1],fish_len_indices[i,2], drop = FALSE]
    tmp_iss_obs <- tmp_iss_obs[rownames(tmp_iss_obs) %in% rownames(tmp_obs)]
    
    # Get expected values from model
    rownames(model$rep$CAL) <- data$yrs # input row names
    tmp_exp <- model$rep$CAL[rownames(model$rep$CAL) %in% rownames(tmp_obs),,fish_len_indices[i,1],fish_len_indices[i,2]] # index to correct expected dimensions
    tmp_exp <- t(apply(X = tmp_exp, MARGIN = 1, FUN = function(x) x / sum(x))) # normalize
    data_yrs <- nrow(tmp_obs)
    
    # Set up reweighting vectors
    exp_bar <- vector() # mean expected
    obs_bar <- vector() # mean observed
    v_y <- vector() # variance
    w_denom <- vector() # weight factor in denominator
    
    for(y in 1:data_yrs) {
      exp_bar[y] <- sum(data$lens * tmp_exp[y,]) # get mean pred comps
      obs_bar[y] <- sum(data$lens * tmp_obs[y,]) # get mean obs comps
      v_y[y]<-sum(data$lens^2*tmp_exp[y,])-exp_bar[y]^2 # get variance
      w_denom[y]<-(obs_bar[y]-exp_bar[y])/sqrt(v_y[y]/tmp_iss_obs[y]) # get weights
    } # end y loop
    
    new_fish_len_wts[fish_len_indices[i,1],fish_len_indices[i,2]] <- 1 / var(w_denom) # return weights
    
  } # end i loop
  
  # Survey Ages ------------------------------------------------------------
  
  # Get weights with data 
  new_srv_age_wts <- matrix(nrow = nrow(data$Wt_SrvAgeComps), ncol = ncol(data$Wt_SrvAgeComps)) # matrix to store new weights
  srv_age_indices <- which(!is.na(data$Wt_SrvAgeComps), arr.ind = TRUE)
  
  for(i in 1:nrow(srv_age_indices)) { # loop through data sources for fishery ages with data
    # Filter to years and fleets with data
    tmp_obs <- data$ObsSrvAgeComps[,,srv_age_indices[i,1],srv_age_indices[i,2]]
    tmp_obs <- na.omit(tmp_obs, drop = FALSE)
    tmp_iss_obs <- data$ISS_SrvAgeComps[,srv_age_indices[i,1],srv_age_indices[i,2], drop = FALSE]
    tmp_iss_obs <- tmp_iss_obs[rownames(tmp_iss_obs) %in% rownames(tmp_obs)]
    
    # Get expected values from model
    rownames(model$rep$SrvIAA) <- data$yrs # input row names
    tmp_exp <- model$rep$SrvIAA[rownames(model$rep$SrvIAA) %in% rownames(tmp_obs),,srv_age_indices[i,1],srv_age_indices[i,2]] # index to correct expected dimensions
    tmp_exp <- t(apply(X = tmp_exp, MARGIN = 1, FUN = function(x) x / sum(x))) # normalize
    data_yrs <- nrow(tmp_obs)
    
    # Set up reweighting vectors
    exp_bar <- vector() # mean expected
    obs_bar <- vector() # mean observed
    v_y <- vector() # variance
    w_denom <- vector() # weight factor in denominator
    
    for(y in 1:data_yrs) {
      exp_bar[y] <- sum(data$ages * tmp_exp[y,]) # get mean pred comps
      obs_bar[y] <- sum(data$ages * tmp_obs[y,]) # get mean obs comps
      v_y[y]<-sum(data$ages^2*tmp_exp[y,])-exp_bar[y]^2 # get variance
      w_denom[y]<-(obs_bar[y]-exp_bar[y])/sqrt(v_y[y]/tmp_iss_obs[y]) # get weights
    } # end y loop
    
    new_srv_age_wts[srv_age_indices[i,1],srv_age_indices[i,2]] <- 1 / var(w_denom) # return weights
    
  } # end i loop
  
  # Survey Lengths ------------------------------------------------------------
  
  # Get weights with data 
  new_srv_len_wts <- matrix(nrow = nrow(data$Wt_SrvLenComps), ncol = ncol(data$Wt_SrvLenComps)) # matrix to store new weights
  srv_len_indices <- which(!is.na(data$Wt_SrvLenComps), arr.ind = TRUE)
  
  for(i in 1:nrow(srv_len_indices)) { # loop through data sources for fishery ages with data
    # Filter to years and fleets with data
    tmp_obs <- data$ObsSrvLenComps[,,srv_len_indices[i,1],srv_len_indices[i,2]]
    tmp_obs <- na.omit(tmp_obs, drop = FALSE)
    tmp_iss_obs <- data$ISS_SrvLenComps[,srv_len_indices[i,1],srv_len_indices[i,2], drop = FALSE]
    tmp_iss_obs <- tmp_iss_obs[rownames(tmp_iss_obs) %in% rownames(tmp_obs)]
    
    # Get expected values from model
    rownames(model$rep$SrvIAL) <- data$yrs # input row names
    tmp_exp <- model$rep$SrvIAL[rownames(model$rep$SrvIAL) %in% rownames(tmp_obs),,srv_len_indices[i,1],srv_len_indices[i,2]] # index to correct expected dimensions
    tmp_exp <- t(apply(X = tmp_exp, MARGIN = 1, FUN = function(x) x / sum(x))) # normalize
    data_yrs <- nrow(tmp_obs)
    
    # Set up reweighting vectors
    exp_bar <- vector() # mean expected
    obs_bar <- vector() # mean observed
    v_y <- vector() # variance
    w_denom <- vector() # weight factor in denominator
    
    for(y in 1:data_yrs) {
      exp_bar[y] <- sum(data$lens * tmp_exp[y,]) # get mean pred comps
      obs_bar[y] <- sum(data$lens * tmp_obs[y,]) # get mean obs comps
      v_y[y]<-sum(data$lens^2*tmp_exp[y,])-exp_bar[y]^2 # get variance
      w_denom[y]<-(obs_bar[y]-exp_bar[y])/sqrt(v_y[y]/tmp_iss_obs[y]) # get weights
    } # end y loop
    
    new_srv_len_wts[srv_len_indices[i,1],srv_len_indices[i,2]] <- 1 / var(w_denom) # return weights
    
  } # end i loop
  
  return(list(new_fish_age_wts = new_fish_age_wts, new_fish_len_wts = new_fish_len_wts,
              new_srv_age_wts = new_srv_age_wts, new_srv_len_wts = new_srv_len_wts))
  
} # end function
