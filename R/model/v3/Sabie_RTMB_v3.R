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

# version 3 - (M.LH Cheng)
# Coded in spatial dimensions
# Parameters (mean recruitment, recruitment devs, initial age devs, 
# selectivity, composition likelihood parameters, catchability, 
# mean fishing mortality, fishing mortlaity deviates) can be estimated spatially
# Incorporated options to allow for estimation of movement parameters across
# years, ages, and sexes


sabie_RTMB = function(pars) {
  
  require(RTMB); require(here)
  source(here("R", "model", "v3", "Get_Selex_v3.R")) # selectivity options
  source(here("R", "model", "v3", "Get_Comp_Likelihoods_v3.R")) # selectivity options
  source(here("R", "model", "ddirmult.R")) # dirichlet multinomial

  RTMB::getAll(pars, data) # load in starting values and data

# Model Set Up (Containers) -----------------------------------------------
  n_ages = length(ages) # number of ages
  n_yrs = length(years) # number of years
  n_lens = length(lens) # number of lengths
  
  # Recruitment stuff
  Rec = array(0, dim = c(n_regions, n_yrs)) # Recruitment
  R0 = rep(0, n_regions) # R0 or mean recruitment
  sigmaR2_early = exp(ln_sigmaR_early)^2 # recruitment variability for early period
  sigmaR2_late = exp(ln_sigmaR_late)^2 # recruitment variability for late period
  
  # Population Dynamics 
  NAA = array(data = 0, dim = c(n_regions, n_yrs + 1, n_ages, n_sexes)) # Numbers at age
  ZAA = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes)) # Total mortality at age
  SAA = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes)) # Survival at age
  SAA_mid = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes)) # Survival at age (midpoint of the year)
  natmort = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes)) # natural mortaltity at age
  Total_Biom = array(0, dim = c(n_regions, n_yrs)) # Total biomass
  SSB = array(0, dim = c(n_regions, n_yrs)) # Spawning stock biomass
  
  # Movement Stuff
  init_iter = n_ages * 10 # Number of times to iterate to equilibrium when movement occurs
  Init_NAA = array(data = 0, dim = c(init_iter, n_regions, n_ages, n_sexes)) # Initial numbers at age
  Movement = array(data = 0, dim = c(n_regions, n_regions, n_yrs, n_ages, n_sexes)) # movement "matrix"
  
  # Fishery Processes
  init_F = init_F_prop * exp(ln_F_mean[1]) # initial F for age structure
  Fmort = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Fishing mortality scalar
  FAA = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets)) # Fishing mortality at age
  CAA = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets)) # Catch at age
  CAL = array(data = 0, dim = c(n_regions, n_yrs, n_lens, n_sexes, n_fish_fleets)) # Catch at length
  PredCatch = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Predicted catch in weight
  PredFishIdx = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Predicted fishery index 
  fish_sel = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_fish_fleets)) # Fishery selectivity
  fish_q = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Fishery catchability
  ESS_FishAgeComps = array(data = 0, dim = c(n_regions, n_yrs, n_sexes, n_fish_fleets)) # Effective sample size for fishery ages (wts * ISS)
  ESS_FishLenComps = array(data = 0, dim = c(n_regions, n_yrs, n_sexes, n_fish_fleets)) # Effective sample size for fishery lengths (wts * ISS)

  # Survey Processes
  SrvIAA = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_srv_fleets)) # Survey index at age
  SrvIAL = array(data = 0, dim = c(n_regions, n_yrs, n_lens, n_sexes, n_srv_fleets)) # Survey index at length
  PredSrvIdx = array(0, dim = c(n_regions, n_yrs, n_srv_fleets)) # Predicted survey index 
  srv_sel = array(data = 0, dim = c(n_regions, n_yrs, n_ages, n_sexes, n_srv_fleets)) # Survey selectivity
  srv_q = array(0, dim = c(n_regions, n_yrs, n_srv_fleets)) # Survey catchability
  ESS_SrvAgeComps = array(data = 0, dim = c(n_regions, n_yrs, n_sexes, n_srv_fleets)) # Effective sample size for survey ages (wts * ISS)
  ESS_SrvLenComps = array(data = 0, dim = c(n_regions, n_yrs, n_sexes, n_srv_fleets)) # Effective sample size for survey lengths (wts * ISS)

  # Likelihoods
  Catch_nLL = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Fishery Catch Likelihoods
  FishIdx_nLL = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Fishery Index Likelihoods
  FishAgeComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_sexes, n_fish_fleets)) # Fishery Age Comps Likelihoods
  FishLenComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_sexes, n_fish_fleets)) # Fishery Length Comps Likelihoods
  SrvIdx_nLL = array(0, dim = c(n_regions, n_yrs, n_srv_fleets)) # Survey Index Likelihoods
  SrvAgeComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_sexes, n_srv_fleets)) # Survey Age Comps Likelihoods
  SrvLenComps_nLL = array(data = 0, dim = c(n_regions, n_yrs, n_sexes, n_srv_fleets)) # Survey Length Comps Likelihoods

  # Penalties and Priors
  Fmort_Pen = array(0, dim = c(n_regions, n_yrs, n_fish_fleets)) # Fishing Mortality Deviation penalty
  Rec_nLL = array(0, dim = c(n_regions, n_yrs)) # Recruitment penalty
  Init_Rec_nLL = array(0, dim = c(n_regions, n_ages - 2)) # Initial Recruitment penalty
  bias_ramp = rep(0, n_yrs) # bias ramp from Methot and Taylor 2011
  M_Pen = 0 # Penalty/Prior for natural mortality
  sel_Pen = 0 # Penalty for selectivity deviations
  jnLL = 0 # Joint negative log likelihood
  
# Model Process Equations -------------------------------------------------
  ## Movement Parameters (Set up) --------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(a in 1:n_ages) {
        for(s in 1:n_sexes) {
          move_tmp = rep(0, n_regions) # temporary movement vector to store values (from - to)
          move_tmp[r] = 0 # Input zero for residency
          move_tmp[-r] = move_pars[r,,y,a,s] # Input parameter values for emigration rates (only if there more than 1 region)
          if(use_fixed_movement == 0) Movement[r,,y,a,s] = exp(move_tmp) / sum(exp(move_tmp)) # multinomial logit transform (basically a softmax) - estimated movement
          if(use_fixed_movement == 1) Movement[r,,y,a,s] = Fixed_Movement[r,,y,a,s] # fixed movement matrix
        } # end s loop
      } # end a loop
    } # end y loop
  } # end r loop

  ## Fishery Selectivity -----------------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {
        fish_sel_blk_idx = fish_sel_blocks[r,y,f] # Get fishery selectivity block index
        for(s in 1:n_sexes) {
          # if selectivity is either a time block or constant
          if(cont_tv_fish_sel[r,f] == 0) tmp_fish_sel_vec = ln_fish_fixed_sel_pars[r,,fish_sel_blk_idx,s,f] # extract temporary selectivity parameters
          # if selectivity is iid or random walk time varying
          if(cont_tv_fish_sel[r,f] %in% c(1,2)) {
            tmp_fish_sel_vec = c(ln_fish_fixed_sel_pars[r,1,fish_sel_blk_idx,s,f] + ln_fishsel_dev1[r,y,s,f],
                                 ln_fish_fixed_sel_pars[r,2,fish_sel_blk_idx,s,f] + ln_fishsel_dev2[r,y,s,f])
          } # end iid or random walk selectivity
          fish_sel[r,y,,s,f] = Get_Selex(Selex_Model = fish_sel_model[r,y,f],
                                         ln_Pars = tmp_fish_sel_vec,
                                         Age = ages) # Calculate selectivity
        } # end s loop
      } # end f loop
    } # end y loop
  } # end r loop

  ## Survey Selectivity ------------------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(sf in 1:n_srv_fleets) {
        srv_sel_blk_idx = srv_sel_blocks[r,y,sf] # Get survey selectivity block index
        for(s in 1:n_sexes) {
          tmp_srv_sel_vec = ln_srv_fixed_sel_pars[r,,srv_sel_blk_idx,s,sf] # extract temporary selectivity parameters
          srv_sel[r,y,,s,sf] = Get_Selex(Selex_Model = srv_sel_model[r,y,sf],
                                         ln_Pars = tmp_srv_sel_vec,
                                         Age = ages) # Calculate selectivity
        } # end s loop
      } # end sf loop
    } # end y loop
  } # end r loop


  ## Mortality ---------------------------------------------------------------
  for(y in 1:n_yrs) {
    for(a in 1:n_ages) {
      for(s in 1:n_sexes) {
        for(r in 1:n_regions) {
          
          # Fishing Mortality at Age calculations
          for(f in 1:n_fish_fleets) {
              if(UseCatch[r,y,f] == 0) {
                Fmort[r,y,f] = 0 # Set F to zero when no catch data
                FAA[r,y,a,s,f] = 0
              } else {
                if(Catch_Type[y,f] == 0 && est_all_regional_F == 0) Fmort[r,y,f] = exp(ln_F_mean_AggCatch[f] + ln_F_devs_AggCatch[y,f]) # If catch is aggregated across regions
                else Fmort[r,y,f] = exp(ln_F_mean[r,f] + ln_F_devs[r,y,f]) # Fully selected F
                # Fmort[r,y,f] = Fmort_dat[r,y,f]
                FAA[r,y,a,s,f] = Fmort[r,y,f] * fish_sel[r,y,a,s,f] # Fishing mortality at age
              }
          } # f loop

          # Population Mortality and Survival
          if(s == 1) natmort[r,y,a,s] = exp(ln_M) # get natural mortality (females or single-sex)
          if(s == 2) natmort[r,y,a,s] = exp(ln_M) + M_offset # natural mortality with offset (males)
          ZAA[r,y,a,s] = sum(FAA[r,y,a,s,]) + natmort[r,y,a,s] # Total Mortality at age
          SAA[r,y,a,s] = exp(-ZAA[r,y,a,s]) # Survival at age
          SAA_mid[r,y,a,s] = exp(-0.5 * ZAA[r,y,a,s]) # Survival at age at midpoint of year
          
        } # r loop
      } # s loop
    } # a loop
  } # y loop


  ## Recruitment Stuff: R0 and Bias Ramp (Methot and Taylor) -------------------------------
  # Set up global virgin recruitment (or mean recruitment)
  if(n_regions == 1) R0 = exp(ln_global_R0) # R0 = global R0 if only 1 region
  if(n_regions > 1) R0 = exp(ln_global_R0) * c(1 - sum(R0_prop), R0_prop) # Multiply a global scaling parameter by estimated proportions if more than 1 region
  
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
  # Iterate equilibrium age structure with mortality, followed by movement
  Init_NAA[1,,1,] = R0 * sexratio # initialize equilibrium population first w/ r0
  for(i in 2:init_iter) {
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        # Recruitment
        Init_NAA[i,r,1,s] = R0[r] * sexratio[s]
        # Iterate from previous "virtual time step"
        Init_NAA[i,r,2:n_ages,s] = Init_NAA[i - 1,r,1:(n_ages-1),s] *
                                    exp(-(natmort[r,1,1:(n_ages-1),s] + init_F * fish_sel[r,1,1:(n_ages-1),s,1])) # mortality
        # Accumulate plus group
        Init_NAA[i,r,n_ages,s] = Init_NAA[i,r,n_ages,s] + # individuals from equilibrium that get aged a year
                                  Init_NAA[i-1,r,n_ages,s] * # decrement and accumulate plus group
                                  exp(-(natmort[r,1,n_ages,s] + init_F * fish_sel[r,1,n_ages,s,1])) # mortality
        } # end s loop
      } # end r loop
    # Apply movement after iteration
    for(a in 1:n_ages) for(s in 1:n_sexes) Init_NAA[i,,a,s] = t(Init_NAA[i,,a,s]) %*% Movement[,,1,a,s]
    
    # Recruits don't move
    if(do_recruits_move == 0) {
      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {
          Init_NAA[i,r,1,s] = R0[r] * sexratio[s]
        } # end s loop
      } # end r loop
    } # end if recruits don't move
  } # end i loop
  
  # Apply initial age deviations
  for(r in 1:n_regions) {
    for(s in 1:n_sexes) {
      Init_NAA[init_iter,r,2:(n_ages-1),s] = Init_NAA[init_iter,r,2:(n_ages-1),s] * exp(ln_InitDevs[r,] - sigmaR2_early/2 * bias_ramp[1]) # add in non-equilibrium age structure
      NAA[r,1,2:n_ages,s] = Init_NAA[init_iter,r,2:n_ages,s] # add in plus group
    } # end s loop
  } # end r loop
  

  ## Annual Recruitment ------------------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(s in 1:n_sexes) {
        if(y < sigmaR_switch) NAA[r,y,1,s] = R0[r] * exp(ln_RecDevs[r,y] - sigmaR2_early/2 * bias_ramp[y]) * sexratio[s] # early period recruitment
        if(y >= sigmaR_switch) NAA[r,y,1,s] = R0[r] * exp(ln_RecDevs[r,y] - sigmaR2_late/2 * bias_ramp[y]) * sexratio[s] # late period recruitment
        Rec[r,y] = sum(NAA[r,y,1,]) # get annual recruitment container here
      } # end s loop
    } # end y loop
  } # end r loop


  ## Population Projection ---------------------------------------------------
  for(y in 1:n_yrs) {
    for(r in 1:n_regions) {
      for(s in 1:n_sexes) {
        for(a in 1:n_ages) {
          # Project ages and years forward and then apply movement
          if(a < n_ages) { 
            # Exponential mortality for individuals not in plus group (recruits experience mortality)
            NAA[r,y+1,a+1,s] = NAA[r,y,a,s] * exp(-ZAA[r,y,a,s])
          } else {
            # Accumulate individuals recently "recruited" into plus group and individuals from previous year
            NAA[r,y+1,n_ages,s] = NAA[r,y+1,n_ages,s] + NAA[r,y,n_ages,s] * exp(-ZAA[r,y,a,s])
          } # end else (calculations for plus group)
        } # end a loop
      } # end s loop
    } # end r loop
    # Recruits don't move
    if(do_recruits_move == 0) {
      # Apply movement after ageing processes - start movement at age 2
      for(a in 2:n_ages) for(s in 1:n_sexes) NAA[,y+1,a,s] = t(NAA[,y+1,a,s]) %*% Movement[,,y,a,s]
      for(r in 1:n_regions) NAA[r,y,1,] = R0[r] * exp(ln_RecDevs[r,y] - sigmaR2_late/2 * bias_ramp[y]) * sexratio
    } # end if recruits don't move
    # Recruits move here
    if(do_recruits_move == 1) for(a in 1:n_ages) for(s in 1:n_sexes) NAA[,y+1,a,s] = t(NAA[,y+1,a,s]) %*% Movement[,,y,a,s]
  } # end y loop

  # Get total biomass and SSB here
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      Total_Biom[r,y] = sum(as.vector(NAA[r,y,,]) * as.vector(WAA[r,y,,])) # Total Biomass
      SSB[r,y] = sum(as.vector(NAA[r,y,,1]) * as.vector(WAA[r,y,,1]) * MatAA[r,y,,1]) # Spawning Stock Biomass 
      if(n_sexes == 1) SSB[r,y] = SSB[r,y] * 0.5 # If single sex model, multiply SSB calculations by 0.5 
    } # end y loop
  } # end r loop
  


# Fishery Observation Model -----------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {

        fish_q_blk_idx = fish_q_blocks[r,y,f] # get time-block catchability index
        fish_q[r,y,f] = exp(ln_fish_q[r,fish_q_blk_idx,f]) # Input into fishery catchability container

        for(s in 1:n_sexes) {
          CAA[r,y,,s,f] = FAA[r,y,,s,f] / ZAA[r,y,,s] * NAA[r,y,,s] * (1 - exp(-ZAA[r,y,,s])) # Catch at age (Baranov's)
          if(fit_lengths == 0) CAL[r,y,,s,f] = SizeAgeTrans[r,y,,,s] %*% CAA[r,y,,s,f] # Catch at length
        } # end s loop
        
        PredCatch[r,y,f] = sum(CAA[r,y,,,f] * WAA[r,y,,]) # get total catch
        
        # Get fishery index
        if(fish_idx_type[r,f] == 0) PredFishIdx[r,y,f] = fish_q[r,y,f] * sum(NAA[r,y,,] * SAA_mid[r,y,,] * fish_sel[r,y,,,f]) # abundance
        if(fish_idx_type[r,f] == 1) {
          # Sablefish specific - for fitting Japanese LL RPW cpue fishery as is implemented in the ADMB assessment
          # (i.e., only using female selex to calculate this)
          if(fish_q_blk_idx == 1 && share_sel == 0) PredFishIdx[r,y,f] = fish_q[r,y,f] * sum(NAA[r,y,,] * SAA_mid[r,y,,] * fish_sel[r,y,,1,f] * WAA[r,y,,]) # first time block
          else PredFishIdx[r,y,f] = fish_q[r,y,f] * sum(NAA[r,y,,] * SAA_mid[r,y,,] * fish_sel[r,y,,,f] * WAA[r,y,,]) # for not first time block
        } # weight

      } # end f loop
    } # end y loop
  } # end r loop


  ## Survey Observation Model ------------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(sf in 1:n_srv_fleets) {

        srv_q_blk_idx = srv_q_blocks[r,y,sf] # get time-block catchability index
        srv_q[r,y,sf] = exp(ln_srv_q[r,srv_q_blk_idx,sf]) # Input into survey catchability container

        for(s in 1:n_sexes) {
          SrvIAA[r,y,,s,sf] = NAA[r,y,,s] * SAA_mid[r,y,,s] * srv_sel[r,y,,s,sf] # Survey index at age
          if(fit_lengths == 0) SrvIAL[r,y,,s,sf] = SizeAgeTrans[r,y,,,s] %*% SrvIAA[r,y,,s,sf] # Survey index at length
        } # end s loop

        # Get predicted survey index
        if(srv_idx_type[r,sf] == 0) PredSrvIdx[r,y,sf] = srv_q[r,y,sf] * sum(SrvIAA[r,y,,,sf]) # abundance
        if(srv_idx_type[r,sf] == 1) PredSrvIdx[r,y,sf] = srv_q[r,y,sf] * sum(SrvIAA[r,y,,,sf] * WAA[r,y,,]) # biomass

      } # end sf loop
    } # end y loop
  } # end r loop


# Likelihood Equations -------------------------------------------------------------
  ## Fishery Likelihoods -----------------------------------------------------
    ### Fishery Catches ---------------------------------------------------------
    for(y in 1:n_yrs) {
      for(f in 1:n_fish_fleets) {
        
        # If we have catch data but it's not resolved on the region scale (sum across regions)
        if(UseCatch[1,y,f] == 1 && Catch_Type[y,f] == 0) {
          # ADMB likelihoods
          if(likelihoods == 0) {
            Catch_nLL[1,y,f] = UseCatch[1,y,f] * (log(ObsCatch[1,y,f] + Catch_Constant[f]) -
                                                    log(sum(PredCatch[,y,f]) + Catch_Constant[f]))^2 # SSQ Catch
          } # ADMB likelihoods
          if(likelihoods == 1) {
            Catch_nLL[1,y,f] = UseCatch[1,y,f] -1 * dnorm(log(ObsCatch[1,y,f] + Catch_Constant[f]),
                                                          log(sum(PredCatch[,y,f]) + Catch_Constant[f]),
                                                          exp(ln_sigmaC[1,f]), TRUE)
          } # TMB likelihoods
        } # if no NAs for fishery catches
        
        for(r in 1:n_regions) {

          # If we have catch data and it's resolved on the region scale
          if(UseCatch[r,y,f] == 1 && Catch_Type[y,f] == 1) {
            # ADMB likelihoods
            if(likelihoods == 0) {
              Catch_nLL[r,y,f] = UseCatch[r,y,f] * (log(ObsCatch[r,y,f] + Catch_Constant[f]) -
                                                    log(PredCatch[r,y,f] + Catch_Constant[f]))^2 # SSQ Catch
            } # ADMB likelihoods
            if(likelihoods == 1) {
              Catch_nLL[r,y,f] = UseCatch[r,y,f] -1 * dnorm(log(ObsCatch[r,y,f] + Catch_Constant[f]),
                                                            log(PredCatch[r,y,f] + Catch_Constant[f]),
                                                            exp(ln_sigmaC[r,f]), TRUE)
            } # TMB likelihoods
          } # if no NAs for fishery catches

        } # end r loop
      } # end f loop
    } # end y loop

  
  ### Fishery Indices ---------------------------------------------------------
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {
      for(r in 1:n_regions) {

        if(UseFishIdx[r,y,f] == 1) {
          if(likelihoods == 0)  {
            FishIdx_nLL[r,y,f] = UseFishIdx[r,y,f] * (log(ObsFishIdx[r,y,f] + 1e-10) - log(PredFishIdx[r,y,f] + 1e-10))^2 /
                                 (2 * (ObsFishIdx_SE[r,y,f] / ObsFishIdx[r,y,f])^2) # lognormal fishery index
          } # ADMB likelihoods
          if(likelihoods == 1) {
            FishIdx_nLL[r,y,f] = UseFishIdx[r,y,f] -1 * dnorm(log(ObsFishIdx[r,y,f] + 1e-10),
                                                              log(PredFishIdx[r,y,f] + 1e-10),
                                                              0.35, TRUE)
          } # TMB likelihoods
        } # if no NAs for fishery index

      } # end r loop
    } # end f loop
  } # end y loop

  
  ### Fishery Compositions ------------------------------------------------
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {
      
      # Fishery Age Compositions
      if(sum(UseFishAgeComps[,y,f]) >= 1) {
        FishAgeComps_nLL[,y,,f] = Get_Comp_Likelihoods(
          Exp = CAA[,y,,,f], Obs = ObsFishAgeComps[,y,,,f], # Expected and Observed values
          ISS = ISS_FishAgeComps[,y,,f], Wt_Mltnml = Wt_FishAgeComps[,,f], # Input sample size and multinomial weight
          Comp_Type = FishAgeComps_Type[y,f], Likelihood_Type = FishAgeComps_LikeType[f], # Composition and Likelihood Type
          ln_theta = ln_FishAge_DM_theta[,f], n_regions =  n_regions, n_sexes = n_sexes, age_or_len = 0, AgeingError = AgeingError, # overdispersion par, Number of sexes, regions, age or length comps, and ageing error
          use = UseFishAgeComps[,y,f], n_bins = n_ages)
      } # if we have fishery age comps

      # Fishery Length Compositions
      if(sum(UseFishLenComps[,y,f]) >= 1 && fit_lengths == 0) {
        FishLenComps_nLL[,y,,f] = Get_Comp_Likelihoods(
          Exp = CAL[,y,,,f], Obs = ObsFishLenComps[,y,,,f], # Expected and Observed values
          ISS = ISS_FishLenComps[,y,,f], Wt_Mltnml = Wt_FishLenComps[r,,f], # Input sample size and multinomial weight
          Comp_Type = FishLenComps_Type[y,f], Likelihood_Type = FishLenComps_LikeType[f], # Composition and Likelihood Type
          ln_theta = ln_FishLen_DM_theta[,f], n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1, AgeingError = NA, # overdispersion, Number of sexes, regions age or length comps, and ageing error
          use = UseFishLenComps[,y,f], n_bins = n_lens) 
      } # if we have fishery length comps

    } # end f loop
  } # end y loop

  
  ## Survey Likelihoods ------------------------------------------------------
    ### Survey Indices ---------------------------------------------------------
    for(y in 1:n_yrs) {
      for(sf in 1:n_srv_fleets) {
        for(r in 1:n_regions) {

          if(UseSrvIdx[r,y,sf] == 1) {
            if(likelihoods == 0) {
              SrvIdx_nLL[r,y,sf] = UseSrvIdx[r,y,sf] * (log(ObsSrvIdx[r,y,sf] + 1e-10) - log(PredSrvIdx[r,y,sf] + 1e-10))^2 /
                                  (2 * (ObsSrvIdx_SE[r,y,sf] / ObsSrvIdx[r,y,sf])^2) # lognormal survey index
            } # ADMB likelihoods
            if(likelihoods == 1) {
              SrvIdx_nLL[r,y,sf] = UseSrvIdx[r,y,sf] -1 * dnorm(log(ObsSrvIdx[r,y,sf] + 1e-10),
                                                                log(PredSrvIdx[r,y,sf] + 1e-10),
                                                                (ObsSrvIdx_SE[r,y,sf] / ObsSrvIdx[r,y,sf]) , TRUE)
            } # TMB likelihoods
          } # if no NAs for survey index

        } # end r loop
      } # end sf loop
    } # end y loop

  
    ### Survey Compositions ---------------------------------------------------------
    for(y in 1:n_yrs) {
      for(sf in 1:n_srv_fleets) {

        # Survey Age Compositions
        if(sum(UseSrvAgeComps[,y,sf]) >= 1) {
          SrvAgeComps_nLL[,y,,sf] = Get_Comp_Likelihoods(
            Exp = SrvIAA[,y,,,sf], Obs = ObsSrvAgeComps[,y,,,sf], # Expected and Observed values
            ISS = ISS_SrvAgeComps[,y,,sf], Wt_Mltnml = Wt_SrvAgeComps[,,sf], # Input sample size and multinomial weight
            Comp_Type = SrvAgeComps_Type[y,sf], Likelihood_Type = SrvAgeComps_LikeType[sf], # Composition and Likelihood Type
            ln_theta = ln_SrvAge_DM_theta[,sf], n_regions = n_regions, n_sexes = n_sexes, age_or_len = 0, AgeingError = AgeingError, # overdispersion, Number of sexes, regions, age or length comps, and ageing error
            use = UseSrvAgeComps[,y,sf], n_bins = n_ages)
        } # if we have survey age comps

        # Survey Length Compositions
        if(sum(UseSrvLenComps[,y,sf]) >= 1 && fit_lengths == 0) {
          SrvLenComps_nLL[,y,,sf] = Get_Comp_Likelihoods(
            Exp = SrvIAL[,y,,,sf], Obs = ObsSrvLenComps[,y,,,sf], # Expected and Observed values
            ISS = ISS_SrvLenComps[,y,,sf], Wt_Mltnml = Wt_SrvLenComps[,,sf], # Input sample size and multinomial weight
            Comp_Type = SrvLenComps_Type[y,sf], Likelihood_Type = SrvLenComps_LikeType[sf], # Composition and Likelihood Type
            ln_theta = ln_SrvLen_DM_theta[,sf], n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1, AgeingError = NA, # overdispersion, Number of sexes, regions, age or length comps, and ageing error
            use = UseSrvLenComps[,y,sf], n_bins = n_lens)
        } # if we have survey length comps

      } # end sf loop
    } # end y loop

  
  ## Priors and Penalties ----------------------------------------------------
    ### Fishing Mortality (Penalty) ---------------------------------------------
    for(f in 1:n_fish_fleets) {
      for(y in 1:n_yrs) {
        for(r in 1:n_regions) {
          
          if(UseCatch[r,y,f] == 1) {
            if(likelihoods == 0) {
              if(Catch_Type[y,f] == 0 && est_all_regional_F == 0) Fmort_Pen[1,y,f] = ln_F_devs_AggCatch[y,f]^2 # Use aggregated catch
              else Fmort_Pen[r,y,f] = ln_F_devs[r,y,f]^2 # SSQ ADMB
            } # ADMB likelihoods
            
            if(likelihoods == 1) {
              if(Catch_Type[y,f] == 0 && est_all_regional_F == 0) Fmort_Pen[1,y,f] = -dnorm(ln_F_devs_AggCatch[y,f], 0, 5, TRUE) # Use aggregated catch
              else Fmort_Pen[r,y,f] = -dnorm(ln_F_devs[r,y,f], 0, 5, TRUE) # TMB likelihood
            } # TMB likelihoods
            
          } # end if have catch
        } # end r loop
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
      if(cont_tv_fish_sel[r,f] == 1) {
        for(s in 1:n_sexes) {
          for(y in 1:n_yrs) {
            for(r in 1:n_regions) {
              sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev1[r,y,s,f], 0, exp(ln_fishsel_dev1_sd[r,s,f]), TRUE)
              sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev2[r,y,s,f], 0, exp(ln_fishsel_dev2_sd[r,s,f]), TRUE)
            } # end r loop
          } # end y loop
        } # end s loop
      } # end if for iid selectivity

      # random walk selectivity
      if(cont_tv_fish_sel[r,f] == 2) {
        for(s in 1:n_sexes) {
          for(r in 1:n_regions) {
            sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev1[r,1,s,f], 0, 50, TRUE) # initialize first value w/ large prior (prior sd is just arbitrarily large)
            sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev2[r,1,s,f], 0, 50, TRUE) # initialize first value w/ large prior (prior sd is just arbitrarily large)
            for(y in 2:n_yrs) {
              sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev1[r,y,s,f], ln_fishsel_dev1[r,y-1,s,f], exp(ln_fishsel_dev1_sd[r,s,f]), TRUE)
              sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev2[r,y,s,f], ln_fishsel_dev2[r,y-1,s,f], exp(ln_fishsel_dev2_sd[r,s,f]), TRUE)
            } # end y loop
          } # end r loop
        } # end s loop
      } # end if for random walk selectivity
    } # end f loop

  
    ### Recruitment (Penalty) ----------------------------------------------------
    if(likelihoods == 0) {
      for(r in 1:n_regions) Init_Rec_nLL[r,] = (ln_InitDevs[r,] / exp(ln_sigmaR_early))^2 # initial age structure penalty
      if(sigmaR_switch > 1) for(r in 1:n_regions) for(y in 1:(sigmaR_switch-1)) Rec_nLL[r,y] = (ln_RecDevs[r,y]/exp(ln_sigmaR_early))^2 + bias_ramp[y]*ln_sigmaR_early # early period
      for(r in 1:n_regions) for(y in (sigmaR_switch:(n_yrs-1))) Rec_nLL[r,y] = (ln_RecDevs[r,y]/exp(ln_sigmaR_late))^2 + bias_ramp[y]*ln_sigmaR_late # late period
      Rec_nLL = 0.5 * sum(Rec_nLL); Init_Rec_nLL = 0.5 * sum(Init_Rec_nLL) # multiply by 0.5 and sum
    }

    if(likelihoods == 1) {
      for(r in 1:n_regions) Init_Rec_nLL[r,] = -dnorm(ln_InitDevs[r,], 0, exp(ln_sigmaR_early), TRUE)
      if(sigmaR_switch > 1) for(r in 1:n_regions) for(y in 1:(sigmaR_switch-1)) Rec_nLL[r,y] = -dnorm(ln_RecDevs[r,y], 0, exp(ln_sigmaR_early), TRUE)
      for(r in 1:n_regions) for(y in sigmaR_switch:(n_yrs)) Rec_nLL[r,y] = -dnorm(ln_RecDevs[r,y], 0, exp(ln_sigmaR_late), TRUE)
    } # Note that this penalizes the terminal year rec devs, which is estiamted in this case

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
          (Wt_Rec * sum(Rec_nLL)) + # Recruitment Penalty
          (Wt_Rec * sum(Init_Rec_nLL)) + #  Initial Age Penalty
           sel_Pen; #  selectivity penalty
    
# Report Section ----------------------------------------------------------
  # Biological Processes
  RTMB::REPORT(R0)
  RTMB::REPORT(Init_NAA)
  RTMB::REPORT(NAA)
  RTMB::REPORT(ZAA)
  RTMB::REPORT(natmort)
  RTMB::REPORT(bias_ramp)
  RTMB::REPORT(Movement)

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
  RTMB::REPORT(Total_Biom)
  RTMB::REPORT(SSB)
  RTMB::REPORT(Rec)
  # Report these in log space because can't be < 0
  RTMB::ADREPORT(log(Total_Biom))
  RTMB::ADREPORT(log(SSB))
  RTMB::ADREPORT(log(Rec))
  RTMB::ADREPORT(Movement)
  
  return(jnLL)
} # end function
