sabie_RTMB = function(pars) {
  
  require(RTMB)
  RTMB::getAll(pars, data) # load in starting values and data

# Get_Selex Function ------------------------------------------------------
Get_Selex = function(Selex_Model, ln_Pars, Age) {
  selex = rep(0, length(Age)) # Temporary container vector
  if(Selex_Model == 0) { # logistic selectivity
    # Extract out and exponentiate the parameters here
    a50 = exp(ln_Pars[1]); # a50
    k = exp(ln_Pars[2]); # slope
    selex = 1 / (1 + exp(-k * (Age - a50))) 
  }
  
  if(Selex_Model == 1) { # gamma dome-shaped selectivity 
    # Extract out and exponentiate the parameters here
    amax = exp(ln_Pars[1]) # age at max selex
    delta = exp(ln_Pars[2]) # slope parameter
    # Now, calculate/derive power parameter + selex values
    p = 0.5 * (sqrt( amax^2 + (4 * delta^2)) - amax)
    selex = (Age / amax)^(amax/p) * exp( (amax - Age) / p ) 
  }
  
  if(Selex_Model == 2) { # power function selectivity
    # Extract out and exponentiate the parameters here
    power = exp(ln_Pars[1]); # power parameter
    selex = 1 / Age^power
  }
  
  return(selex)
} # end function
  

# Dirichlet-Multinomial Likelihood ----------------------------------------
# From https://github.com/James-Thorson/CCSRA/blob/main/inst/executables/CCSRA_v9.cpp
ddirmult <- function(obs, pred, Ntotal, ln_theta, give_log = TRUE) {
  # Set up function variables
  n_c = length(obs) # number of categories
  p_exp = pred # expected values container
  p_obs = obs # observed values container
  dirichlet_Parm = exp(ln_theta) * Ntotal # Dirichlet alpha parameters
  
  # set up pdf
  logres = lgamma(Ntotal + 1)
  for(c in 1:n_c) logres = logres - lgamma(Ntotal*p_obs[c]+1) # integration constant
  logres = logres + lgamma(dirichlet_Parm) - lgamma(Ntotal+dirichlet_Parm) # 2nd term in formula
  
  # Summation in 3rd term in formula
  for(c in 1:n_c) {
    logres = logres + lgamma(Ntotal*p_obs[c] + dirichlet_Parm*p_exp[c])
    logres = logres - lgamma(dirichlet_Parm * p_exp[c])
  } # end c
  
  if(give_log == TRUE) return(logres)
  else return(exp(logres))
} # end function


# Model Set Up (Containers) -----------------------------------------------
  n_ages = length(ages) # number of ages
  n_yrs = length(years) # number of years
  n_lens = length(lens) # number of lengths
  
  # Recruitment stuff
  Rec = rep(0, n_yrs) # Recruitment
  sigmaR2_early = exp(ln_sigmaR_early)^2 # recruitment variability for early period
  sigmaR2_late = exp(ln_sigmaR_late)^2 # recruitment variability for late period
  
  # Population Dynamics 
  NAA = array(data = 0, dim = c(n_yrs + 1, n_ages, n_sexes)) # Numbers at age
  ZAA = array(data = 0, dim = c(n_yrs, n_ages, n_sexes)) # Total mortality at age
  SAA = array(data = 0, dim = c(n_yrs, n_ages, n_sexes)) # Survival at age
  SAA_mid = array(data = 0, dim = c(n_yrs, n_ages, n_sexes)) # Survival at age (midpoint of the year)
  natmort = array(data = 0, dim = c(n_yrs, n_ages, n_sexes)) # natural mortaltity at age
  SSB = rep(0, n_yrs) # Spawning stock biomass
  
  # Fishery Processes
  init_F = init_F_prop * exp(ln_F_mean[1]) # initial F for age structure
  Fmort = matrix(data = 0, nrow = n_yrs, ncol = n_fish_fleets) # Fishing mortality scalar
  FAA = array(data = 0, dim = c(n_yrs, n_ages, n_sexes, n_fish_fleets)) # Fishing mortality at age
  CAA = array(data = 0, dim = c(n_yrs, n_ages, n_sexes, n_fish_fleets)) # Catch at age
  CAL = array(data = 0, dim = c(n_yrs, n_lens, n_sexes, n_fish_fleets)) # Catch at length
  PredCatch = matrix(data = 0, nrow = n_yrs, ncol = n_fish_fleets) # Predicted catch in weight
  PredFishIdx = matrix(data = 0, nrow = n_yrs, ncol = n_fish_fleets) # Predicted fishery index 
  fish_sel = array(data = 0, dim = c(n_yrs, n_ages, n_sexes, n_fish_fleets)) # Fishery selectivity
  fish_q = matrix(data = 0, nrow = n_yrs, ncol = n_fish_fleets) # Fishery catchability
  ESS_FishAgeComps = array(data = 0, dim = c(n_yrs, n_sexes, n_fish_fleets)) # Effective sample size for fishery ages (wts * ISS)
  ESS_FishLenComps = array(data = 0, dim = c(n_yrs, n_sexes, n_fish_fleets)) # Effective sample size for fishery lengths (wts * ISS)

  # Survey Processes
  SrvIAA = array(data = 0, dim = c(n_yrs, n_ages, n_sexes, n_srv_fleets)) # Survey index at age
  SrvIAL = array(data = 0, dim = c(n_yrs, n_lens, n_sexes, n_srv_fleets)) # Survey index at length
  PredSrvIdx = matrix(data = 0, nrow = n_yrs, ncol = n_srv_fleets) # Predicted survey index 
  srv_sel = array(data = 0, dim = c(n_yrs, n_ages, n_sexes, n_srv_fleets)) # Survey selectivity
  srv_q = matrix(data = 0, nrow = n_yrs, ncol = n_srv_fleets) # Survey catchability
  ESS_SrvAgeComps = array(data = 0, dim = c(n_yrs, n_sexes, n_srv_fleets)) # Effective sample size for survey ages (wts * ISS)
  ESS_SrvLenComps = array(data = 0, dim = c(n_yrs, n_sexes, n_srv_fleets)) # Effective sample size for survey lengths (wts * ISS)

  # Likelihoods
  Catch_nLL = matrix(data = 0, nrow = n_yrs, ncol = n_fish_fleets) # Fishery Catch Likelihoods
  FishIdx_nLL = matrix(data = 0, nrow = n_yrs, ncol = n_fish_fleets) # Fishery Index Likelihoods
  FishAgeComps_nLL = array(data = 0, dim = c(n_yrs, n_sexes, n_fish_fleets)) # Fishery Age Comps Likelihoods
  FishLenComps_nLL = array(data = 0, dim = c(n_yrs, n_sexes, n_fish_fleets)) # Fishery Length Comps Likelihoods
  FishAgeComps_offset_nLL = array(data = 0, dim = c(n_yrs, n_sexes, n_fish_fleets)) # Fishery Age Comps Likelihoods multinomial offset
  FishLenComps_offset_nLL = array(data = 0, dim = c(n_yrs, n_sexes, n_fish_fleets)) # Fishery Length Comps Likelihoods multinomial offset
  SrvIdx_nLL = matrix(data = 0, nrow = n_yrs, ncol = n_srv_fleets) # Survey Index Likelihoods
  SrvAgeComps_nLL = array(data = 0, dim = c(n_yrs, n_sexes, n_srv_fleets)) # Survey Age Comps Likelihoods
  SrvLenComps_nLL = array(data = 0, dim = c(n_yrs, n_sexes, n_srv_fleets)) # Survey Length Comps Likelihoods
  SrvAgeComps_offset_nLL = array(data = 0, dim = c(n_yrs, n_sexes, n_srv_fleets)) # Survey Age Comps Likelihoods multinomial offset
  SrvLenComps_offset_nLL = array(data = 0, dim = c(n_yrs, n_sexes, n_srv_fleets)) # Survey Length Comps Likelihoods multinomial offset
  
  # Penalties and Priors
  Fmort_Pen = rep(0, n_fish_fleets) # Fishing Mortality Deviation penalty
  Rec_nLL = rep(0, n_yrs - 1) # Recruitment penalty
  Init_Rec_nLL = rep(0, n_ages - 2) # Initial Recruitment penalty
  bias_ramp = rep(0, n_yrs) # bias ramp from Methot and Taylor 2011
  M_Pen = 0 # Penalty/Prior for natural mortality
  jnLL = 0 # Joint negative log likelihood

# Model Process Equations -------------------------------------------------
  ## Fishery Selectivity -----------------------------------------------------
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {
      fish_sel_blk_idx = fish_sel_blocks[y,f] # Get fishery selectivity block index
      for(s in 1:n_sexes) {
        tmp_fish_sel_vec = ln_fish_fixed_sel_pars[,fish_sel_blk_idx,s,f] # extract temporary selectivity parameters
        fish_sel[y,,s,f] = Get_Selex(Selex_Model = fish_sel_model[y,f],
                                       ln_Pars = tmp_fish_sel_vec,
                                       Age = ages) # Calculate selectivity
      } # end s loop
    } # end f loop
  } # end y loop
  
  # TESTING: Fixing selectivity at true selectivity (has a max call)
  # fish_sel[] = fish_sel_dat[]
  
  ## Survey Selectivity ------------------------------------------------------
  for(y in 1:n_yrs) {
    for(sf in 1:n_srv_fleets) {
      srv_sel_blk_idx = srv_sel_blocks[y,sf] # Get survey selectivity block index
      for(s in 1:n_sexes) {
        tmp_srv_sel_vec = ln_srv_fixed_sel_pars[,srv_sel_blk_idx,s,sf] # extract temporary selectivity parameters
        srv_sel[y,,s,sf] = Get_Selex(Selex_Model = srv_sel_model[y,sf],
                                     ln_Pars = tmp_srv_sel_vec,
                                     Age = ages) # Calculate selectivity
      } # end s loop
    } # end sf loop
  } # end y loop


  ## Mortality ---------------------------------------------------------------
  for(y in 1:n_yrs) {
    for(a in 1:n_ages) {
      for(s in 1:n_sexes) {

        # Fishing Mortality at Age calculations
        for(f in 1:n_fish_fleets) {
          if(is.na(ObsCatch[y,f])) {
            Fmort[y,f] = 0 # Set F to zero when no catch data
            FAA[y,a,s,f] = 0
          } else {
            Fmort[y,f] = exp(ln_F_mean[f] + ln_F_devs[y,f]) # Fully selected F
            FAA[y,a,s,f] = Fmort[y,f] * fish_sel[y,a,s,f] # Fishing mortality at age
          }
        } # f loop

        # Population Mortality and Survival
        if(s == 1) natmort[y,a,s] = exp(ln_M) # get natural mortality (females or single-sex)
        if(s == 2) natmort[y,a,s] = exp(ln_M) + M_offset # natural mortality with offset (males)
        ZAA[y,a,s] = sum(FAA[y,a,s,]) + natmort[y,a,s] # Total Mortality at age
        SAA[y,a,s] = exp(-1 * ZAA[y,a,s]) # Survival at age
        SAA_mid[y,a,s] = exp(-0.5 * ZAA[y,a,s]) # Survival at age at midpoint of year

      } # s loop
    } # a loop
  } # y loop

  ## Recruitment Bias Ramp (Methot and Taylor) -------------------------------
  for(y in 1:n_yrs) {
    if(do_rec_bias_ramp == 0) bias_ramp[y] = 1 # don't do bias ramp correction
    if(do_rec_bias_ramp == 1) {
      if(y < bias_year[1] || y == bias_year[4]) bias_ramp[y] = 0 # no bias correction during poor data
      if(y >= bias_year[1] && y < bias_year[2]) bias_ramp[y] = 1 * ((y - bias_year[1]) / (bias_year[2] - bias_year[1])) # ascending limb
      if(y >= bias_year[2] && y < bias_year[3]) bias_ramp[y] = 1 # full bias correction
      if(y >= bias_year[3] && y < bias_year[4]) bias_ramp[y] = 1 * (1 - ((y - bias_year[3]) / (bias_year[4] - bias_year[3]))) #descending limb
    } # if we want to do bias ramp
  } # end y loop


  ## Initial Age Structure ---------------------------------------------------
  init_age_idx = 1:(n_ages - 2) # Get initial age indexing
  for(s in 1:n_sexes) {
    NAA[1,init_age_idx + 1,s] = exp(ln_R0 + ln_InitDevs[init_age_idx] -
                                      (init_age_idx * (natmort[1, init_age_idx + 1, s] +
                                      (init_F * fish_sel[1, init_age_idx + 1, s, 1])))) * sexratio[s] # not plus group
    # Plus group calculations
    NAA[1,n_ages,s] = exp(ln_R0 - ((n_ages - 1) * (natmort[1, n_ages, s] + (init_F * fish_sel[1, n_ages, s, 1]))) ) /
                      (1 - exp(-(natmort[1, n_ages, s] + (init_F * fish_sel[1, n_ages, s, 1])))) * sexratio[s]

  } # end s loop

  ## Annual Recruitment ------------------------------------------------------
  for(y in 1:n_yrs) {
    for(s in 1:n_sexes) {
      if(y < sigmaR_switch) NAA[y,1,s] = exp(ln_R0 + ln_RecDevs[y] - (bias_ramp[y] * sigmaR2_early/2)) * sexratio[s] # early period recruitment
      if(y >= sigmaR_switch && y < n_yrs) NAA[y,1,s] = exp(ln_R0 + ln_RecDevs[y] - (bias_ramp[y] * sigmaR2_late/2)) * sexratio[s] # late period recruitment
      if(y == n_yrs) NAA[y,1,s] = exp(ln_R0) * sexratio[s] # mean recruitment in terminal year
    } # end s loop
    Rec[y] = sum(NAA[y,1,]) # get annual recruitment here
  } # end y loop


  ## Population Projection ---------------------------------------------------
  for(y in 1:n_yrs) {
    for(a in 1:n_ages) {
      for(s in 1:n_sexes) {
        # Project ages and years forward
        if(a < n_ages) NAA[y+1,a+1,s] = NAA[y,a,s] * SAA[y,a,s] # not plus group
        if(a == n_ages) NAA[y+1,n_ages,s] = sum(NAA[y,n_ages,s] * SAA[y,n_ages,s], NAA[y+1,n_ages,s]) # plus group (decrement previous year and add sum projected year)
      } # end s loop
    } # end a loop
    SSB[y] = sum(NAA[y,,1] * WAA[y,,1] * MatAA[y,,1]) # Spawning stock biomass calculation
  } # end y loop
  

  ## Fishery Observation Model -----------------------------------------------
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {

      fish_q_blk_idx = fish_q_blocks[y,f] # get time-block catchability index
      fish_q[y,f] = exp(ln_fish_q[fish_q_blk_idx,f]) # Input into fishery catchability container

      for(s in 1:n_sexes) {
        CAA[y,,s,f] = FAA[y,,s,f] / ZAA[y,,s] * NAA[y,,s] * (1 - exp(-ZAA[y,,s])) # Catch at age (Baranov's)
        CAL[y,,s,f] = SizeAgeTrans[y,,,s] %*% CAA[y,,s,f] # Catch at length
      } # end s loop

      PredCatch[y,f] = sum(CAA[y,,,f] * WAA[y,,]) # get total catch

      # Get fishery index
      if(fish_idx_type[f] == 0) PredFishIdx[y,f] = fish_q[y,f] * sum(NAA[y,,] * SAA_mid[y,,] * fish_sel[y,,,f]) # abundance
      if(fish_idx_type[f] == 1) {
        # Sablefish specific - for fitting Japanese LL RPW cpue fishery as is implemented in the ADMB assessment
        # (i.e., only using female selex to calculate this)
        if(fish_q_blk_idx == 1 && share_sel == 0) PredFishIdx[y,f] = fish_q[y,f] * sum(NAA[y,,] * SAA_mid[y,,] * fish_sel[y,,1,f] * WAA[y,,]) # first time block
        else PredFishIdx[y,f] = fish_q[y,f] * sum(NAA[y,,] * SAA_mid[y,,] * fish_sel[y,,,f] * WAA[y,,]) # for not first time block
      } # weight

    } # end f loop
  } # end y loop


  ## Survey Observation Model ------------------------------------------------
  for(y in 1:n_yrs) {
    for(sf in 1:n_srv_fleets) {

      srv_q_blk_idx = srv_q_blocks[y,sf] # get time-block catchability index
      srv_q[y,sf] = exp(ln_srv_q[srv_q_blk_idx,sf]) # Input into survey catchability container

      for(s in 1:n_sexes) {
        SrvIAA[y,,s,sf] = NAA[y,,s] * srv_sel[y,,s,sf] # Survey index at age
        SrvIAL[y,,s,sf] = SizeAgeTrans[y,,,s] %*% SrvIAA[y,,s,sf] # Survey index at length
      } # end s loop

      # Get predicted survey index
      if(srv_idx_type[sf] == 0) PredSrvIdx[y,sf] = srv_q[y,sf] * sum(SrvIAA[y,,,sf] * SAA_mid[y,,]) # abundance
      if(srv_idx_type[sf] == 1) PredSrvIdx[y,sf] = srv_q[y,sf] * sum(SrvIAA[y,,,sf] * SAA_mid[y,,] * WAA[y,,]); # biomass

    } # end sf loop
  } # end y loop


# Likelihood Equations -------------------------------------------------------------
  # Fishery Likelihoods -----------------------------------------------------
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {

      ### Fishery Catches ---------------------------------------------------------
      if(!is.na(ObsCatch[y,f])) {
        # ADMB likelihoods
        if(likelihoods == 0) {
          Catch_nLL[y,f] = UseCatch[y,f] * (log(ObsCatch[y,f] + Catch_Constant[f]) - 
                           log(PredCatch[y,f] + Catch_Constant[f]))^2 # SSQ Catch
        } # ADMB likelihoods
        if(likelihoods == 1) {
          Catch_nLL[y,f] = UseCatch[y,f] -1 * dnorm(log(ObsCatch[y,f] + Catch_Constant[f]), 
                           log(PredCatch[y,f] + Catch_Constant[f]), exp(ln_obscatch_sigma[f]), TRUE)
        } # TMB likelihoods
      } # if no NAs for fishery catches

      ### Fishery Indices ---------------------------------------------------------
      if(!is.na(ObsFishIdx[y,f])) {
        if(likelihoods == 0)  {
          FishIdx_nLL[y,f] = UseFishIdx[y,f] * (log(ObsFishIdx[y,f] + 1e-4) - log(PredFishIdx[y,f] + 1e-4))^2 /
                             (2 * (ObsFishIdx_SE[y,f] / ObsFishIdx[y,f])^2) # lognormal fishery index
        } # ADMB likelihoods
        if(likelihoods == 1) {
          FishIdx_nLL[y,f] = UseFishIdx[y,f] -1 * dnorm(log(ObsFishIdx[y,f] + 1e-4), 
                             log(PredFishIdx[y,f] + 1e-4), 0.3, TRUE)
        } # TMB likelihoods
      } # if no NAs for fishery index

      ### Fishery Age Compositions ------------------------------------------------
      if(!is.na(sum(ObsFishAgeComps[y,,1,f]))) {
        tmp_Agg_CAA = CAA[y,,,f] / matrix(data = colSums(CAA[y,,,f]), nrow = nrow(CAA[y,,,f]), ncol = ncol(CAA[y,,,f]), byrow = TRUE) # normalize
        tmp_Agg_CAA = t(as.vector(rowSums(tmp_Agg_CAA) / n_sexes)) %*% AgeingError # aggregate and apply ageing error
        tmp_Agg_CAA = tmp_Agg_CAA / sum(tmp_Agg_CAA) # normalize
        tmp_ObsFishAgeComps = ObsFishAgeComps[y,,1,f] / sum(ObsFishAgeComps[y,,1,f]) # normalize
        ESS_FishAgeComps[y,1,f] = ISS_FishAgeComps[y,1,f] * Wt_FishAgeComps[1,f] # Calculate effective sample size for fishery age compositions
        FishAgeComps_offset_nLL[y,1,f] = -1 * ESS_FishAgeComps[y,1,f] * sum(((ObsFishAgeComps[y,,1,f]  + 0.001)*log(ObsFishAgeComps[y,,1,f]  + 0.001))) # multinomial offset
        FishAgeComps_nLL[y,1,f] = -1 * UseFishAgeComps[y,f] * ESS_FishAgeComps[y,1,f] * sum(((tmp_ObsFishAgeComps + 0.001) * log(tmp_Agg_CAA + 0.001))) # Compute ADMB Multinomial Likelihood
      } # if no NAs for fishery age comps

      ### Fishery Length Compositions ---------------------------------------------
      for(s in 1:n_sexes) {
        if(!is.na(sum(ObsFishLenComps[y,,s,f]))) {
          tmp_CAL = CAL[y,,s,f] / sum(CAL[y,,s,f]) # temporary variable and normalize
          tmp_ObsFishLenComps = ObsFishLenComps[y,,s,f] / sum(ObsFishLenComps[y,,s,f]) # normalize
          ESS_FishLenComps[y,s,f] = ISS_FishLenComps[y,s,f] * Wt_FishLenComps[s,f] # Calculate effective sample size for fishery length compositions
          FishLenComps_offset_nLL[y,s,f] = -1 * ESS_FishLenComps[y,s,f] * sum(((ObsFishLenComps[y,,s,f] + 0.001)*log(ObsFishLenComps[y,,s,f] + 0.001))) # multinomial offset
          FishLenComps_nLL[y,s,f] = -1 * UseFishLenComps[y,f] * ESS_FishLenComps[y,s,f] * sum(((tmp_ObsFishLenComps + 0.001) * log(tmp_CAL + 0.001))) # Compute ADMB Multinomial Likelihood
        } # if no NAs for fishery length comps
      } # end s loop

    } # end f loop
  } # end y loop


  ## Survey Likelihoods ------------------------------------------------------
  for(y in 1:n_yrs) {
    for(sf in 1:n_srv_fleets) {

      ### Survey Indices ---------------------------------------------------------
      if(!is.na(ObsSrvIdx[y,sf])) {
        if(likelihoods == 0) {
          SrvIdx_nLL[y,sf] = UseSrvIdx[y,sf] * (log(ObsSrvIdx[y,sf] + 1e-4) - log(PredSrvIdx[y,sf] + 1e-4))^2 /
                            (2 * (ObsSrvIdx_SE[y,sf] / ObsSrvIdx[y,sf])^2) # lognormal survey index
        } # ADMB likelihoods
        if(likelihoods == 1) {
          SrvIdx_nLL[y,sf] = UseSrvIdx[y,sf] -1 * dnorm(log(ObsSrvIdx[y,sf] + 1e-4), 
                             log(PredSrvIdx[y,sf] + 1e-4), ObsSrvIdx_SE[y,sf] / ObsSrvIdx[y,sf], TRUE)
        } # TMB likelihoods
      } # if no NAs for survey index

      ### Survey Age Compositions ------------------------------------------------
      if(!is.na(sum(ObsSrvAgeComps[y,,1,sf]))) {
        tmp_Agg_SrvIAA = SrvIAA[y,,,sf] / matrix(data = colSums(SrvIAA[y,,,sf]), nrow = nrow(SrvIAA[y,,,sf]), ncol = ncol(SrvIAA[y,,,sf]), byrow = TRUE) # normalize
        tmp_Agg_SrvIAA = t(as.vector(rowSums(tmp_Agg_SrvIAA) / n_sexes)) %*% AgeingError # aggregate and apply ageing error
        tmp_Agg_SrvIAA = tmp_Agg_SrvIAA / sum(tmp_Agg_SrvIAA) # normalize
        tmp_ObsSrvAgeComps = ObsSrvAgeComps[y,,1,sf] / sum(ObsSrvAgeComps[y,,1,sf]) # normalize
        ESS_SrvAgeComps[y,1,sf] = ISS_SrvAgeComps[y,1,sf] * Wt_SrvAgeComps[1,sf] # Calculate effective sample size for survey age compositions
        SrvAgeComps_offset_nLL[y,1,sf] = -1 * ESS_SrvAgeComps[y,1,sf] * sum(((ObsSrvAgeComps[y,,1,sf] + 0.001)*log(ObsSrvAgeComps[y,,1,sf] + 0.001))) # multinomial offset
        SrvAgeComps_nLL[y,1,sf] = -1 * UseSrvAgeComps[y,sf] * ESS_SrvAgeComps[y,1,sf] * sum(((tmp_ObsSrvAgeComps + 0.001) * log(tmp_Agg_SrvIAA + 0.001))) # Compute ADMB Multinomial Likelihood
      } # if no NAs for survey age comps

      # ### Survey Length Compositions ---------------------------------------------
      for(s in 1:n_sexes) {
        if(!is.na(sum(ObsSrvLenComps[y,,s,sf]))) {
          tmp_SrvIAL = SrvIAL[y,,s,sf] / sum(SrvIAL[y,,s,sf]) # temporary variable and normalize
          tmp_ObsSrvLenComps = ObsSrvLenComps[y,,s,sf] / sum(ObsSrvLenComps[y,,s,sf]) # normalize
          ESS_SrvLenComps[y,s,sf] = ISS_SrvLenComps[y,s,sf] * Wt_SrvLenComps[s,sf] # Calculate effective sample size for survey length compositions
          SrvLenComps_offset_nLL[y,s,sf] = -1 * ESS_SrvLenComps[y,s,sf] * sum(((ObsSrvLenComps[y,,s,sf] + 0.001)*log(ObsSrvLenComps[y,,s,sf] + 0.001))) # multinomial offset
          SrvLenComps_nLL[y,s,sf] = -1 * UseSrvLenComps[y,sf] * ESS_SrvLenComps[y,s,sf] * sum(((tmp_ObsSrvLenComps + 0.001) * log(tmp_SrvIAL + 0.001))) # Compute ADMB Multinomial Likelihood
        } # if no NAs for survey length comps
      } # end s loop
    } # end sf loop
  } # end y loop


  ## Priors and Penalties ----------------------------------------------------
    ### Fishing Mortality (Penalty) ---------------------------------------------
    for(f in 1:n_fish_fleets) {
      for(y in 1:n_yrs) {
        if(!is.na(ObsCatch[y,f])) {
          if(likelihoods == 0) Fmort_Pen[f] = Fmort_Pen[f] + ln_F_devs[y, f]^2 # SSQ ADMB
          if(likelihoods == 1) Fmort_Pen[f] = Fmort_Pen[f] -1 * dnorm(ln_F_devs[y, f], 0, exp(ln_sigmaF[f]), TRUE) # TMB likelihoods
        } # end if
      } # y loop
    } # f loop
  
    ### Natural Mortality (Prior) -----------------------------------------------
    if(Use_M_prior == 1) {
      if(likelihoods == 0) M_Pen = (ln_M - log(M_prior[1]))^2 / (2 * (M_prior[2])^2) # ADMB likelihood
      if(likelihoods == 1) M_Pen = -1 * dnorm(ln_M, log(M_prior[1]), M_prior[2], TRUE) # TMB likelihood
    } # end if
  
    ### Recruitment (Penalty) ----------------------------------------------------
    Init_Rec_nLL = (ln_InitDevs / exp(ln_sigmaR_early))^2 # initial age structure penalty
    for(y in 1:(sigmaR_switch-1)) Rec_nLL[y] = (ln_RecDevs[y]/exp(ln_sigmaR_early))^2 + bias_ramp[y]*ln_sigmaR_early # early period
    for(y in (sigmaR_switch:(n_yrs-1))) Rec_nLL[y] = (ln_RecDevs[y]/exp(ln_sigmaR_late))^2 + bias_ramp[y]*ln_sigmaR_late # late period

    # Apply likelihood weights here and compute joint negative log likelihood
    jnLL = (Wt_Catch * sum(Catch_nLL)) + # Catch likelihoods
           (Wt_FishIdx * sum(FishIdx_nLL)) + # Fishery Index likelihood
           sum(FishAgeComps_nLL - FishAgeComps_offset_nLL)  + # Fishery Age likelihood
           sum(FishLenComps_nLL - FishLenComps_offset_nLL)+ # Fishery Length likelihood
          (Wt_SrvIdx * sum(SrvIdx_nLL)) + # Survey Index likelihood
           sum(SrvAgeComps_nLL - SrvAgeComps_offset_nLL)+ # Survey Age likelihood
           sum(SrvLenComps_nLL - SrvLenComps_offset_nLL)+ # Survey Length likelihood
          (Wt_F * sum(Fmort_Pen)) + # Fishery Mortality Penalty
           M_Pen + # Natural Mortality Prior (Penalty)
          (Wt_Rec * 0.5 * sum(Rec_nLL)) + # Recruitment Penalty
          (Wt_Rec * 0.5 * sum(Init_Rec_nLL)); #  Initial Age Penalty


# Report Section ----------------------------------------------------------
  # Biological Processes
  RTMB::REPORT(NAA)
  RTMB::REPORT(ZAA)
  RTMB::REPORT(natmort)
  RTMB::REPORT(bias_ramp)

  # Fishery Processes
  RTMB::REPORT(Fmort)
  RTMB::REPORT(FAA)
  RTMB::REPORT(CAA)
  RTMB::REPORT(CAL)
  RTMB::REPORT(PredCatch)
  RTMB::REPORT(PredFishIdx)
  RTMB::REPORT(fish_sel)
  RTMB::REPORT(fish_q)

  # Survey Processes
  RTMB::REPORT(PredSrvIdx)
  RTMB::REPORT(srv_sel)
  RTMB::REPORT(srv_q)
  RTMB::REPORT(SrvIAA)
  RTMB::REPORT(SrvIAL)

  # Likelihoods
  RTMB::REPORT(Catch_nLL)
  RTMB::REPORT(FishIdx_nLL)
  RTMB::REPORT(SrvIdx_nLL)
  RTMB::REPORT(FishAgeComps_nLL)
  RTMB::REPORT(FishAgeComps_offset_nLL)
  RTMB::REPORT(SrvAgeComps_nLL)
  RTMB::REPORT(SrvAgeComps_offset_nLL)
  RTMB::REPORT(FishLenComps_nLL)
  RTMB::REPORT(FishLenComps_offset_nLL)
  RTMB::REPORT(SrvLenComps_nLL)
  RTMB::REPORT(SrvLenComps_offset_nLL)
  RTMB::REPORT(M_Pen)
  RTMB::REPORT(Fmort_Pen)
  RTMB::REPORT(Rec_nLL)
  RTMB::REPORT(Init_Rec_nLL)
  RTMB::REPORT(Rec_nLL)
  RTMB::REPORT(jnLL)

  # Effective Sample Sizes
  RTMB::REPORT(ESS_FishAgeComps)
  RTMB::REPORT(ESS_SrvAgeComps)
  RTMB::REPORT(ESS_FishLenComps)
  RTMB::REPORT(ESS_SrvLenComps)

  # Report for derived quantities
  RTMB::REPORT(SSB)
  RTMB::REPORT(Rec)
  RTMB::ADREPORT(SSB)
  RTMB::ADREPORT(Rec)
  # can take a while if reported (since its doing delta method on an array of n_yrs, n_ages, n_sexes, n_fish_fleets)
  # RTMB::ADREPORT(fish_sel) 

  return(jnLL)
} # end function