# version 1 - 12/22/24 (M.LH Cheng)
# Bridge model 23.5 from ADMB to RTMB
# Changed code to be more modular, accommodating any number of fishery and survey fleets
# Rectified errors in fitting to length composition data (normalize proportions at length after conversion from age-length matrix)
# Changed survey composition data to be calculated using survival midyear
# Added options for continuous time-varying selectivity
# Added options for TMB / R-like likelihoods (e.g., dnorm) instead of bespoke likelihoods 
# Added options for dirichlet multinomial likelihood

# version 2 - 12/23/24 (M.LH Cheng)
# Incorporated options to fit age and length composition data as sex-aggregated, split by sex (no sex ratio 
# information), and jointly by sex (implicit sex ratio information)
# Added in option for Dirichlet Multinomial likelihood

sabie_RTMB = function(pars) {
  
  require(RTMB); require(here)
  source(here("R", "model", "Get_Selex.R")) # selectivity options
  source(here("R", "model", "Get_Comp_Likelihoods.R")) # selectivity options
  source(here("R", "model", "ddirmult.R")) # dirichlet multinomial
  RTMB::getAll(pars, data) # load in starting values and data

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
  SrvIdx_nLL = matrix(data = 0, nrow = n_yrs, ncol = n_srv_fleets) # Survey Index Likelihoods
  SrvAgeComps_nLL = array(data = 0, dim = c(n_yrs, n_sexes, n_srv_fleets)) # Survey Age Comps Likelihoods
  SrvLenComps_nLL = array(data = 0, dim = c(n_yrs, n_sexes, n_srv_fleets)) # Survey Length Comps Likelihoods
  SrvLenComps_offset_nLL = array(data = 0, dim = c(n_yrs, n_sexes, n_srv_fleets)) # Survey Length Comps Likelihoods multinomial offset
  
  # Penalties and Priors
  Fmort_Pen = rep(0, n_fish_fleets) # Fishing Mortality Deviation penalty
  Rec_nLL = rep(0, n_yrs - 1) # Recruitment penalty
  Init_Rec_nLL = rep(0, n_ages - 2) # Initial Recruitment penalty
  bias_ramp = rep(0, n_yrs) # bias ramp from Methot and Taylor 2011
  M_Pen = 0 # Penalty/Prior for natural mortality
  sel_Pen = 0 # Penalty for selectivity deviations
  jnLL = 0 # Joint negative log likelihood
  
# Model Process Equations -------------------------------------------------
  ## Fishery Selectivity -----------------------------------------------------
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {
      fish_sel_blk_idx = fish_sel_blocks[y,f] # Get fishery selectivity block index
      for(s in 1:n_sexes) {
        
        # if selectivity is either a time block or constant
        if(cont_tv_fish_sel[f] == 0) tmp_fish_sel_vec = ln_fish_fixed_sel_pars[,fish_sel_blk_idx,s,f] # extract temporary selectivity parameters
        # if selectivity is iid or random walk time varying
        if(cont_tv_fish_sel[f] %in% c(1,2)) {
          tmp_fish_sel_vec = c(ln_fish_fixed_sel_pars[1,fish_sel_blk_idx,s,f] + ln_fishsel_dev1[y,s,f],
                               ln_fish_fixed_sel_pars[2,fish_sel_blk_idx,s,f] + ln_fishsel_dev2[y,s,f])
        } # end iid or random walk selectivity

        fish_sel[y,,s,f] = Get_Selex(Selex_Model = fish_sel_model[y,f],
                                       ln_Pars = tmp_fish_sel_vec,
                                       Age = ages) # Calculate selectivity
      } # end s loop
    } # end f loop
  } # end y loop

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
        SrvIAA[y,,s,sf] = NAA[y,,s] * SAA_mid[y,,s] * srv_sel[y,,s,sf] # Survey index at age
        SrvIAL[y,,s,sf] = SizeAgeTrans[y,,,s] %*% SrvIAA[y,,s,sf] # Survey index at length
      } # end s loop

      # Get predicted survey index
      if(srv_idx_type[sf] == 0) PredSrvIdx[y,sf] = srv_q[y,sf] * sum(SrvIAA[y,,,sf]) # abundance
      if(srv_idx_type[sf] == 1) PredSrvIdx[y,sf] = srv_q[y,sf] * sum(SrvIAA[y,,,sf] * WAA[y,,]); # biomass

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
                             log(PredFishIdx[y,f] + 1e-4), 0.35, TRUE)
        } # TMB likelihoods
      } # if no NAs for fishery index

      ### Fishery Age Compositions ------------------------------------------------
      if(UseFishAgeComps[y,f] == 1) {
        FishAgeComps_nLL[y,,f] = Get_Comp_Likelihoods(
          Exp = CAA[y,,,f], Obs = ObsFishAgeComps[y,,,f], # Expected and Observed values
          ISS = ISS_FishAgeComps[y,,f], Wt_Mltnml = Wt_FishAgeComps[,f], # Input sample size and multinomial weight
          Comp_Type = FishAgeComps_Type[f], Likelihood_Type = FishAgeComps_LikeType[f], # Composition and Likelihood Type
          ln_theta = ln_FishAge_DM_theta[f], n_sexes = n_sexes, age_or_len = 0, AgeingError = AgeingError # overdispersion par, Number of sexes, age or length comps, and ageing error
          ) 
      } # if we have fishery age comps

      ### Fishery Length Compositions ---------------------------------------------
      if(UseFishLenComps[y,f] == 1) {
        FishLenComps_nLL[y,,f] = Get_Comp_Likelihoods(
          Exp = CAL[y,,,f], Obs = ObsFishLenComps[y,,,f], # Expected and Observed values
          ISS = ISS_FishLenComps[y,,f], Wt_Mltnml = Wt_FishLenComps[,f], # Input sample size and multinomial weight
          Comp_Type = FishLenComps_Type[f], Likelihood_Type = FishLenComps_LikeType[f], # Composition and Likelihood Type
          ln_theta = ln_FishLen_DM_theta[f], n_sexes = n_sexes, age_or_len = 1, AgeingError = NA # overdispersion, Number of sexes, age or length comps, and ageing error
        ) 
      } # if we have fishery age comps
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
                             log(PredSrvIdx[y,sf] + 1e-4), (ObsSrvIdx_SE[y,sf] / ObsSrvIdx[y,sf]) , TRUE)
        } # TMB likelihoods
      } # if no NAs for survey index

      ### Survey Age Compositions ------------------------------------------------
      if(UseSrvAgeComps[y,sf] == 1) {
        SrvAgeComps_nLL[y,,sf] = Get_Comp_Likelihoods(
          Exp = SrvIAA[y,,,sf], Obs = ObsSrvAgeComps[y,,,sf], # Expected and Observed values
          ISS = ISS_SrvAgeComps[y,,sf], Wt_Mltnml = Wt_SrvAgeComps[,sf], # Input sample size and multinomial weight
          Comp_Type = SrvAgeComps_Type[sf], Likelihood_Type = SrvAgeComps_LikeType[sf], # Composition and Likelihood Type
          ln_theta = ln_SrvAge_DM_theta[sf], n_sexes = n_sexes, age_or_len = 0, AgeingError = AgeingError # overdispersion, Number of sexes, age or length comps, and ageing error
        ) 
      } # if we have survey age comps

      # ### Survey Length Compositions ---------------------------------------------
      if(UseSrvLenComps[y,sf] == 1) {
        SrvLenComps_nLL[y,,sf] = Get_Comp_Likelihoods(
          Exp = SrvIAL[y,,,sf], Obs = ObsSrvLenComps[y,,,sf], # Expected and Observed values
          ISS = ISS_SrvLenComps[y,,sf], Wt_Mltnml = Wt_SrvLenComps[,sf], # Input sample size and multinomial weight
          Comp_Type = SrvLenComps_Type[sf], Likelihood_Type = SrvLenComps_LikeType[sf], # Composition and Likelihood Type
          ln_theta = ln_SrvLen_DM_theta[sf], n_sexes = n_sexes, age_or_len = 1, AgeingError = NA # overdispersion, Number of sexes, age or length comps, and ageing error
        ) 
      } # if we have survey length comps
    } # end sf loop
  } # end y loop


  ## Priors and Penalties ----------------------------------------------------
    ### Fishing Mortality (Penalty) ---------------------------------------------
    for(f in 1:n_fish_fleets) {
      for(y in 1:n_yrs) {
        if(!is.na(ObsCatch[y,f])) {
          if(likelihoods == 0) Fmort_Pen[f] = Fmort_Pen[f] + ln_F_devs[y, f]^2 # SSQ ADMB
          if(likelihoods == 1) Fmort_Pen[f] = Fmort_Pen[f] -1 * dnorm(ln_F_devs[y, f], exp(ln_F_mean[f]), exp(ln_sigmaF[f]), TRUE) # TMB likelihoods
        } # end if
      } # y loop
    } # f loop
  
    ### Natural Mortality (Prior) -----------------------------------------------
    if(Use_M_prior == 1) {
      if(likelihoods == 0) M_Pen = (ln_M - log(M_prior[1]))^2 / (2 * (M_prior[2])^2) # ADMB likelihood
      if(likelihoods == 1) M_Pen = -1 * dnorm(ln_M, log(M_prior[1]), M_prior[2], TRUE) # TMB likelihood
    } # end if
  

    ### Selectivity (Penalty) ---------------------------------------------------
    for(f in 1:n_fish_fleets) {
      # iid selectivity
      if(cont_tv_fish_sel[f] == 1) {
        for(s in 1:n_sexes) {
          for(y in 1:n_yrs) {
            sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev1[y,s,f], 0, exp(ln_fishsel_dev1_sd[s,f]), TRUE)
            sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev2[y,s,f], 0, exp(ln_fishsel_dev2_sd[s,f]), TRUE)
          } # end y loop
        } # end s loop
      } # end if for iid selectivity
      
      # random walk selectivity
      if(cont_tv_fish_sel[f] == 2) {
        for(s in 1:n_sexes) {
          sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev1[1,s,f], 0, 50, TRUE) # initialize first value w/ large prior (prior sd is just arbitrarily large)
          sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev2[1,s,f], 0, 50, TRUE) # initialize first value w/ large prior (prior sd is just arbitrarily large)
          for(y in 2:n_yrs) {sd 
            sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev1[y,s,f], ln_fishsel_dev1[y-1,s,f], exp(ln_fishsel_dev1_sd[s,f]), TRUE)
            sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev2[y,s,f], ln_fishsel_dev2[y-1,s,f], exp(ln_fishsel_dev2_sd[s,f]), TRUE)          
          } # end y loop
        } # end s loop
      } # end if for random walk selectivity
    } # end f loop

    ### Recruitment (Penalty) ----------------------------------------------------
    Init_Rec_nLL = (ln_InitDevs / exp(ln_sigmaR_early))^2 # initial age structure penalty
    for(y in 1:(sigmaR_switch-1)) Rec_nLL[y] = (ln_RecDevs[y]/exp(ln_sigmaR_early))^2 + bias_ramp[y]*ln_sigmaR_early # early period
    for(y in (sigmaR_switch:(n_yrs-1))) Rec_nLL[y] = (ln_RecDevs[y]/exp(ln_sigmaR_late))^2 + bias_ramp[y]*ln_sigmaR_late # late period

    # Apply likelihood weights here and compute joint negative log likelihood
    jnLL = (Wt_Catch * sum(Catch_nLL)) + # Catch likelihoods
           (Wt_FishIdx * sum(FishIdx_nLL)) + # Fishery Index likelihood
           (Wt_SrvIdx * sum(SrvIdx_nLL)) + # Survey Index likelihood
           sum(FishAgeComps_nLL) + # Fishery Age likelihood
           sum(FishLenComps_nLL) + # Fishery Length likelihood
           sum(SrvAgeComps_nLL) + # Survey Age likelihood
           sum(SrvLenComps_nLL) + # Survey Length likelihood
          (Wt_F * sum(Fmort_Pen)) + # Fishery Mortality Penalty
           M_Pen + # Natural Mortality Prior (Penalty)
          (Wt_Rec * 0.5 * sum(Rec_nLL)) + # Recruitment Penalty
          (Wt_Rec * 0.5 * sum(Init_Rec_nLL)) + #  Initial Age Penalty
           sel_Pen; #  selectivity penalty
    
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
  RTMB::REPORT(SrvAgeComps_nLL)
  RTMB::REPORT(FishLenComps_nLL)
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

  return(jnLL)
} # end function