# Purpose: Function to run Francis reweighting

#' Title
#'
#' @param data List of data inputs
#' @param model List of model outputs (includes report file)
#'
#' @return
#' @export
#'
#' @examples
francis_rwgt <- function(data, model) {
  
  # Set up
  require(here)
  source(here("R", "functions", "Restrc_Comps.R"))
  
  # Get indexing
  n_fish_fleets <- data$n_fish_fleets
  n_srv_fleets <- data$n_srv_fleets
  n_sexes <- data$n_sexes

# Fishery Ages ------------------------------------------------------------

  # Get weights with data 
  new_fish_age_wts <- matrix(nrow = nrow(data$Wt_FishAgeComps), ncol = ncol(data$Wt_FishAgeComps)) # matrix to store new weights
  
  for(f in 1:n_fish_fleets) {
    data_yrs <- which(data$UseFishAgeComps[,f] == 1) # get years with data for a given fleet
    # Set up reweighting vectors
    exp_bar <- matrix(NA, length(data_yrs), n_sexes) # mean expected
    obs_bar <- matrix(NA, length(data_yrs), n_sexes)  # mean observed
    v_y <- matrix(NA, length(data_yrs), n_sexes)  # variance
    w_denom <- matrix(NA, length(data_yrs), n_sexes)  # weight factor in denominator
    
    for(y in data_yrs) {
      # Get temporary comps list
      comps_tmp_list = Restrc_Comps_Francis(Exp = model$rep$CAA[y,,,f], Obs = data$ObsFishAgeComps[y,,,f], 
                                            ISS = data$ISS_FishAgeComps[y,,f], Wt_Mltnml = data$Wt_FishAgeComps[,f], 
                                            Comp_Type = data$FishAgeComps_Type[f], n_sexes = n_sexes, age_or_len = 0, AgeingError = data$AgeingError)
      
      # Extract out temporary variables
      tmp_iss_obs = data$ISS_FishAgeComps[y,,f] # temporary ISS
      tmp_exp = comps_tmp_list$Exp # temporary expected values
      tmp_obs = comps_tmp_list$Obs # temporary observed values
      yr_alt_idx = which(data_yrs == y) # get indexing to start from 1
      
      # Aggregated comps fitting
      if(data$FishAgeComps_Type[f] == 0) {
        exp_bar[yr_alt_idx,1] = sum(data$ages * as.vector(tmp_exp[,1])) # get mean pred comps
        obs_bar[yr_alt_idx,1] = sum(data$ages * as.vector(tmp_obs[,1])) # get mean obs comps
        v_y[yr_alt_idx,1] = sum(data$ages^2*tmp_exp[,1])-exp_bar[yr_alt_idx,1]^2 # get variance
        w_denom[yr_alt_idx,1] = (obs_bar[yr_alt_idx,1]-exp_bar[yr_alt_idx,1])/sqrt(v_y[yr_alt_idx,1]/tmp_iss_obs[1]) # get weights      
        } # end if for 'Aggregated"
      
      # Split comps fitting
      if(data$FishAgeComps_Type[f] == 1) {
        for(s in 1:n_sexes) {
          exp_bar[yr_alt_idx,s] = sum(data$ages * as.vector(tmp_exp[,s])) # get mean pred comps
          obs_bar[yr_alt_idx,s] = sum(data$ages * as.vector(tmp_obs[,s])) # get mean obs comps
          v_y[yr_alt_idx,s] = sum(data$ages^2*tmp_exp[,s])-exp_bar[yr_alt_idx,s]^2 # get variance
          w_denom[yr_alt_idx,s] = (obs_bar[yr_alt_idx,s]-exp_bar[yr_alt_idx,s])/sqrt(v_y[yr_alt_idx,s]/tmp_iss_obs[s]) # get weights     
        } # end s loop
      } # end if for 'Split'
        
      # Joint comps fitting
      if(data$FishAgeComps_Type[f] == 2) {
        exp_bar[yr_alt_idx,1] = sum(rep(data$ages, n_sexes) * as.vector(tmp_exp)) # get mean pred comps
        obs_bar[yr_alt_idx,1] = sum(rep(data$ages, n_sexes) * as.vector(tmp_obs)) # get mean obs comps
        v_y[yr_alt_idx,1] = sum(rep(data$ages, n_sexes)^2*as.vector(tmp_exp))-exp_bar[yr_alt_idx,1]^2 # get variance
        w_denom[yr_alt_idx,1] = (obs_bar[yr_alt_idx,1]-exp_bar[yr_alt_idx,1])/sqrt(v_y[yr_alt_idx,1]/tmp_iss_obs[1]) # get weights          
      } # end if for 'Joint'
    } # end y loop
    
    # Calculate weights
    # Aggregated comps fitting
    if(data$FishAgeComps_Type[f] == 0) new_fish_age_wts[1,f] <- 1 / var(w_denom[,1]) # return weights
    # Split comps fitting
    if(data$FishAgeComps_Type[f] == 1) for(s in 1:n_sexes) new_fish_age_wts[s,f] <- 1/var(w_denom[,s])
    # Joint comps fitting
    if(data$FishAgeComps_Type[f] == 2) new_fish_age_wts[1,f] <- 1 / var(w_denom[,1]) # return weights
  } # end f loop

# Fishery Lengths ------------------------------------------------------------
  
  # Get weights with data 
  new_fish_len_wts <- matrix(nrow = nrow(data$Wt_FishLenComps), ncol = ncol(data$Wt_FishLenComps)) # matrix to store new weights
  
  for(f in 1:n_fish_fleets) {
    data_yrs <- which(data$UseFishLenComps[,f] == 1) # get years with data for a given fleet
    # Set up reweighting vectors
    exp_bar <- matrix(NA, length(data_yrs), n_sexes) # mean expected
    obs_bar <- matrix(NA, length(data_yrs), n_sexes)  # mean observed
    v_y <- matrix(NA, length(data_yrs), n_sexes)  # variance
    w_denom <- matrix(NA, length(data_yrs), n_sexes)  # weight factor in denominator
    
    for(y in data_yrs) {
      # Get temporary comps list
      comps_tmp_list = Restrc_Comps_Francis(Exp = model$rep$CAL[y,,,f], Obs = data$ObsFishLenComps[y,,,f], 
                                            ISS = data$ISS_FishLenComps[y,,f], Wt_Mltnml = data$Wt_FishLenComps[,f], 
                                            Comp_Type = data$FishLenComps_Type[f], n_sexes = n_sexes, age_or_len = 1, AgeingError = data$AgeingError)
      
      # Extract out temporary variables
      tmp_iss_obs = data$ISS_FishLenComps[y,,f] # temporary ISS
      tmp_exp = comps_tmp_list$Exp # temporary expected values
      tmp_obs = comps_tmp_list$Obs # temporary observed values
      yr_alt_idx = which(data_yrs == y) # get indexing to start from 1
      
      # Aggregated comps fitting
      if(data$FishLenComps_Type[f] == 0) {
        exp_bar[yr_alt_idx,1] = sum(data$lens * as.vector(tmp_exp[,1])) # get mean pred comps
        obs_bar[yr_alt_idx,1] = sum(data$lens * as.vector(tmp_obs[,1])) # get mean obs comps
        v_y[yr_alt_idx,1] = sum(data$lens^2*tmp_exp[,1])-exp_bar[yr_alt_idx,1]^2 # get variance
        w_denom[yr_alt_idx,1] = (obs_bar[yr_alt_idx,1]-exp_bar[yr_alt_idx,1])/sqrt(v_y[yr_alt_idx,1]/tmp_iss_obs[1]) # get weights      
      } # end if for 'Aggregated"
      
      # Split comps fitting
      if(data$FishLenComps_Type[f] == 1) {
        for(s in 1:n_sexes) {
          exp_bar[yr_alt_idx,s] = sum(data$lens * as.vector(tmp_exp[,s])) # get mean pred comps
          obs_bar[yr_alt_idx,s] = sum(data$lens * as.vector(tmp_obs[,s])) # get mean obs comps
          v_y[yr_alt_idx,s] = sum(data$lens^2*tmp_exp[,s])-exp_bar[yr_alt_idx,s]^2 # get variance
          w_denom[yr_alt_idx,s] = (obs_bar[yr_alt_idx,s]-exp_bar[yr_alt_idx,s])/sqrt(v_y[yr_alt_idx,s]/tmp_iss_obs[s]) # get weights     
        } # end s loop
      } # end if for 'Split'
      
      # Joint comps fitting
      if(data$FishLenComps_Type[f] == 2) {
        exp_bar[yr_alt_idx,1] = sum(rep(data$lens, n_sexes) * as.vector(tmp_exp)) # get mean pred comps
        obs_bar[yr_alt_idx,1] = sum(rep(data$lens, n_sexes) * as.vector(tmp_obs)) # get mean obs comps
        v_y[yr_alt_idx,1] = sum(rep(data$lens, n_sexes)^2*as.vector(tmp_exp))-exp_bar[yr_alt_idx,1]^2 # get variance
        w_denom[yr_alt_idx,1] = (obs_bar[yr_alt_idx,1]-exp_bar[yr_alt_idx,1])/sqrt(v_y[yr_alt_idx,1]/tmp_iss_obs[1]) # get weights          
      } # end if for 'Joint'
    } # end y loop
    
    # Calculate weights
    # Aggregated comps fitting
    if(data$FishLenComps_Type[f] == 0) new_fish_len_wts[1,f] <- 1 / var(w_denom[,1]) # return weights
    # Split comps fitting
    if(data$FishLenComps_Type[f] == 1) for(s in 1:n_sexes) new_fish_len_wts[s,f] <- 1/var(w_denom[,s])
    # Joint comps fitting
    if(data$FishLenComps_Type[f] == 2) new_fish_len_wts[1,f] <- 1 / var(w_denom[,1]) # return weights
  } # end f loop
  
# Survey Ages ------------------------------------------------------------
  # Get weights with data 
  new_srv_age_wts <- matrix(nrow = nrow(data$Wt_SrvAgeComps), ncol = ncol(data$Wt_SrvAgeComps)) # matrix to store new weights
  
  for(sf in 1:n_srv_fleets) {
    data_yrs <- which(data$UseSrvAgeComps[,sf] == 1) # get years with data for a given fleet
    # Set up reweighting vectors
    exp_bar <- matrix(NA, length(data_yrs), n_sexes) # mean expected
    obs_bar <- matrix(NA, length(data_yrs), n_sexes)  # mean observed
    v_y <- matrix(NA, length(data_yrs), n_sexes)  # variance
    w_denom <- matrix(NA, length(data_yrs), n_sexes)  # weight factor in denominator
    
    for(y in data_yrs) {
      # Get temporary comps list
      comps_tmp_list = Restrc_Comps_Francis(Exp = model$rep$SrvIAA[y,,,sf], Obs = data$ObsSrvAgeComps[y,,,sf], 
                                            ISS = data$ISS_SrvAgeComps[y,,sf], Wt_Mltnml = data$Wt_SrvAgeComps[,sf], 
                                            Comp_Type = data$SrvAgeComps_Type[sf], n_sexes = n_sexes, age_or_len = 0, AgeingError = data$AgeingError)
      
      # Extract out temporary variables
      tmp_iss_obs = data$ISS_SrvAgeComps[y,,sf] # temporary ISS
      tmp_exp = comps_tmp_list$Exp # temporary expected values
      tmp_obs = comps_tmp_list$Obs # temporary observed values
      yr_alt_idx = which(data_yrs == y) # get indexing to start from 1
      
      # Aggregated comps fitting
      if(data$SrvAgeComps_Type[sf] == 0) {
        exp_bar[yr_alt_idx,1] = sum(data$ages * as.vector(tmp_exp[,1])) # get mean pred comps
        obs_bar[yr_alt_idx,1] = sum(data$ages * as.vector(tmp_obs[,1])) # get mean obs comps
        v_y[yr_alt_idx,1] = sum(data$ages^2*tmp_exp[,1])-exp_bar[yr_alt_idx,1]^2 # get variance
        w_denom[yr_alt_idx,1] = (obs_bar[yr_alt_idx,1]-exp_bar[yr_alt_idx,1])/sqrt(v_y[yr_alt_idx,1]/tmp_iss_obs[1]) # get weights      
      } # end if for 'Aggregated"
      
      # Split comps fitting
      if(data$SrvAgeComps_Type[sf] == 1) {
        for(s in 1:n_sexes) {
          exp_bar[yr_alt_idx,s] = sum(data$ages * as.vector(tmp_exp[,s])) # get mean pred comps
          obs_bar[yr_alt_idx,s] = sum(data$ages * as.vector(tmp_obs[,s])) # get mean obs comps
          v_y[yr_alt_idx,s] = sum(data$ages^2*tmp_exp[,s])-exp_bar[yr_alt_idx,s]^2 # get variance
          w_denom[yr_alt_idx,s] = (obs_bar[yr_alt_idx,s]-exp_bar[yr_alt_idx,s])/sqrt(v_y[yr_alt_idx,s]/tmp_iss_obs[s]) # get weights     
        } # end s loop
      } # end if for 'Split'
      
      # Joint comps fitting
      if(data$SrvAgeComps_Type[sf] == 2) {
        exp_bar[yr_alt_idx,1] = sum(rep(data$ages, n_sexes) * as.vector(tmp_exp)) # get mean pred comps
        obs_bar[yr_alt_idx,1] = sum(rep(data$ages, n_sexes) * as.vector(tmp_obs)) # get mean obs comps
        v_y[yr_alt_idx,1] = sum(rep(data$ages, n_sexes)^2*as.vector(tmp_exp))-exp_bar[yr_alt_idx,1]^2 # get variance
        w_denom[yr_alt_idx,1] = (obs_bar[yr_alt_idx,1]-exp_bar[yr_alt_idx,1])/sqrt(v_y[yr_alt_idx,1]/tmp_iss_obs[1]) # get weights          
      } # end if for 'Joint'
    } # end y loop
    
    # Calculate weights
    # Aggregated comps fitting
    if(data$SrvAgeComps_Type[sf] == 0) new_srv_age_wts[1,sf] <- 1 / var(w_denom[,1]) # return weights
    # Split comps fitting
    if(data$SrvAgeComps_Type[sf] == 1) for(s in 1:n_sexes) new_srv_age_wts[s,sf] <- 1/var(w_denom[,s])
    # Joint comps fitting
    if(data$SrvAgeComps_Type[sf] == 2) new_srv_age_wts[1,sf] <- 1 / var(w_denom[,1]) # return weights
  } # end f loop
  
# Survey Lengths ------------------------------------------------------------
  # Get weights with data 
  new_srv_len_wts <- matrix(nrow = nrow(data$Wt_SrvLenComps), ncol = ncol(data$Wt_SrvLenComps)) # matrix to store new weights
  
  for(sf in 1:n_srv_fleets) {
    data_yrs <- which(data$UseSrvLenComps[,sf] == 1) # get years with data for a given fleet
    # Set up reweighting vectors
    exp_bar <- matrix(NA, length(data_yrs), n_sexes) # mean expected
    obs_bar <- matrix(NA, length(data_yrs), n_sexes)  # mean observed
    v_y <- matrix(NA, length(data_yrs), n_sexes)  # variance
    w_denom <- matrix(NA, length(data_yrs), n_sexes)  # weight factor in denominator
    
    for(y in data_yrs) {
      # Get temporary comps list
      comps_tmp_list = Restrc_Comps_Francis(Exp = model$rep$SrvIAL[y,,,sf], Obs = data$ObsSrvLenComps[y,,,sf], 
                                            ISS = data$ISS_SrvLenComps[y,,sf], Wt_Mltnml = data$Wt_SrvLenComps[,sf], 
                                            Comp_Type = data$SrvLenComps_Type[sf], n_sexes = n_sexes, age_or_len = 1, AgeingError = data$AgeingError)
      
      # Extract out temporary variables
      tmp_iss_obs = data$ISS_SrvLenComps[y,,sf] # temporary ISS
      tmp_exp = comps_tmp_list$Exp # temporary expected values
      tmp_obs = comps_tmp_list$Obs # temporary observed values
      yr_alt_idx = which(data_yrs == y) # get indexing to start from 1
      
      # Aggregated comps fitting
      if(data$SrvLenComps_Type[sf] == 0) {
        exp_bar[yr_alt_idx,1] = sum(data$lens * as.vector(tmp_exp[,1])) # get mean pred comps
        obs_bar[yr_alt_idx,1] = sum(data$lens * as.vector(tmp_obs[,1])) # get mean obs comps
        v_y[yr_alt_idx,1] = sum(data$lens^2*tmp_exp[,1])-exp_bar[yr_alt_idx,1]^2 # get variance
        w_denom[yr_alt_idx,1] = (obs_bar[yr_alt_idx,1]-exp_bar[yr_alt_idx,1])/sqrt(v_y[yr_alt_idx,1]/tmp_iss_obs[1]) # get weights      
      } # end if for 'Aggregated"
      
      # Split comps fitting
      if(data$SrvLenComps_Type[sf] == 1) {
        for(s in 1:n_sexes) {
          exp_bar[yr_alt_idx,s] = sum(data$lens * as.vector(tmp_exp[,s])) # get mean pred comps
          obs_bar[yr_alt_idx,s] = sum(data$lens * as.vector(tmp_obs[,s])) # get mean obs comps
          v_y[yr_alt_idx,s] = sum(data$lens^2*tmp_exp[,s])-exp_bar[yr_alt_idx,s]^2 # get variance
          w_denom[yr_alt_idx,s] = (obs_bar[yr_alt_idx,s]-exp_bar[yr_alt_idx,s])/sqrt(v_y[yr_alt_idx,s]/tmp_iss_obs[s]) # get weights     
        } # end s loop
      } # end if for 'Split'
      
      # Joint comps fitting
      if(data$SrvLenComps_Type[sf] == 2) {
        exp_bar[yr_alt_idx,1] = sum(rep(data$lens, n_sexes) * as.vector(tmp_exp)) # get mean pred comps
        obs_bar[yr_alt_idx,1] = sum(rep(data$lens, n_sexes) * as.vector(tmp_obs)) # get mean obs comps
        v_y[yr_alt_idx,1] = sum(rep(data$lens, n_sexes)^2*as.vector(tmp_exp))-exp_bar[yr_alt_idx,1]^2 # get variance
        w_denom[yr_alt_idx,1] = (obs_bar[yr_alt_idx,1]-exp_bar[yr_alt_idx,1])/sqrt(v_y[yr_alt_idx,1]/tmp_iss_obs[1]) # get weights          
      } # end if for 'Joint'
    } # end y loop
    
    # Calculate weights
    # Aggregated comps fitting
    if(data$SrvLenComps_Type[sf] == 0) new_srv_len_wts[1,sf] <- 1 / var(w_denom[,1]) # return weights
    # Split comps fitting
    if(data$SrvLenComps_Type[sf] == 1) for(s in 1:n_sexes) new_srv_len_wts[s,sf] <- 1/var(w_denom[,s])
    # Joint comps fitting
    if(data$SrvLenComps_Type[sf] == 2) new_srv_len_wts[1,sf] <- 1 / var(w_denom[,1]) # return weights
  } # end sf loop
  
  return(list(new_fish_age_wts = new_fish_age_wts, new_fish_len_wts = new_fish_len_wts,
              new_srv_age_wts = new_srv_age_wts, new_srv_len_wts = new_srv_len_wts))
  
} # end function
