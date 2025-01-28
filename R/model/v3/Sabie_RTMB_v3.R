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
  init_iter = n_ages * 5 # Number of times to iterate to equilibrium when movement occurs
  # Set up initial age structure
  Init_NAA = array(0, dim = c(n_regions, n_ages, n_sexes))
  Init_NAA_next_year = Init_NAA
  Movement = array(data = 0, dim = c(n_regions, n_regions, n_yrs, n_ages, n_sexes)) # movement "matrix"
  
  # Tagging Stuff
  Tags_Avail = array(data = 0, dim = c(max_tag_liberty + 1, n_tag_cohorts, n_regions, n_ages, n_sexes)) # Tags availiable for recapture
  Tag_Reporting = array(data = 0, dim = c(n_regions, n_yrs)) # Tag reporting rate 
  Pred_Tag_Recap = array(data = 0, dim = c(max_tag_liberty, n_tag_cohorts, n_regions, n_ages, n_sexes)) # predicted recaptures

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
  Tag_nLL = array(data = 0, dim = c(max_tag_liberty, n_tag_cohorts, n_regions, n_ages, n_sexes)) # Tagging Likelihoods
  
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


  ## Recruitment: R0 and Bias Ramp (Methot and Taylor) -------------------------------
  if(n_regions > 1) {
    R0_trans = c(0, R0_prop) # set up vector for transformation
    R0_trans = exp(R0_trans) / sum(exp(R0_trans)) # do multinomial logit
    R0 = exp(ln_global_R0) * R0_trans # Multiply a global scaling parameter by estimated proportions
  } else R0 = exp(ln_global_R0)
  
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
  # Set up initial equilibrium age structure, with cumulative sum of selectivity incorporated
  for(r in 1:n_regions) {
    for(s in 1:n_sexes) {
      tmp_cumsum_Z = cumsum(natmort[r,1,1:(n_ages-1),s] + init_F * fish_sel[r,1,1:(n_ages-1),s,1])
      Init_NAA[r,,s] = c(R0[r], R0[r] * exp(-tmp_cumsum_Z)) * sexratio[s]
    } # end s loop
  } # end r loop
  
  # Apply annual cycle and iterate to equilibrium 
  for(i in 1:init_iter) {
    for(s in 1:n_sexes) {
      Init_NAA_next_year[,1,s] = R0 * sexratio[s] # recruitment
      for(a in 1:n_ages) Init_NAA_next_year[,a,s] = t(Init_NAA_next_year[,a,s]) %*% Movement[,,1,a,s] # movement
      # ageing and mortality 
      Init_NAA_next_year[,2:n_ages,s] = Init_NAA[,1:(n_ages-1),s] * 
                                        exp(-(natmort[,1,1:(n_ages-1),s] + 
                                             (init_F * fish_sel[,1,1:(n_ages-1),s,1])))
      # accumulate plus group
      Init_NAA_next_year[,n_ages,s] = (Init_NAA_next_year[,n_ages,s] * exp(-(natmort[,1,n_ages,s] + (init_F * fish_sel[,1,n_ages,s,1])))) +
                                      (Init_NAA[,n_ages,s] * exp(-(natmort[,1,n_ages,s] + (init_F * fish_sel[,1,n_ages,s,1]))))
      Init_NAA = Init_NAA_next_year # iterate to next cycle
    } # end s loop
  } # end i loop
  
  # Apply initial age deviations
  for(r in 1:n_regions) {
    for(s in 1:n_sexes) {
      Init_NAA[r,2:(n_ages-1),s] = Init_NAA[r,2:(n_ages-1),s] * exp(ln_InitDevs[r,] - sigmaR2_early/2 * bias_ramp[1]) # add in non-equilibrium age structure
      NAA[r,1,2:n_ages,s] = Init_NAA[r,2:n_ages,s] # add in plus group
    } # end s loop
  } # end r loop

  # Current Assessment Approach -- FLAG!
  # init_age_idx = 1:(n_ages - 2) # Get initial age indexing
  # for(r in 1:n_regions) {
  #   for(s in 1:n_sexes) {
  #     NAA[r,1,init_age_idx + 1,s] = R0[r] * exp(ln_InitDevs[r,init_age_idx] -
  #                                       (init_age_idx * (natmort[r,1, init_age_idx + 1, s] +
  #                                                          (init_F * fish_sel[r,1, init_age_idx + 1, s, 1])))) * sexratio[s] # not plus group
  #     # Plus group calculations
  #     NAA[r,1,n_ages,s] = R0[r] *  exp( - ((n_ages - 1) * (natmort[r,1, n_ages, s] + (init_F * fish_sel[r,1, n_ages, s, 1]))) ) /
  #       (1 - exp(-(natmort[r,1, n_ages, s] + (init_F * fish_sel[r,1, n_ages, s, 1])))) * sexratio[s]
  # 
  #   } # end s loop
  # }

  ## Annual Recruitment ------------------------------------------------------
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      for(s in 1:n_sexes) {
        if(y < sigmaR_switch) NAA[r,y,1,s] = R0[r] * exp(ln_RecDevs[r,y] - (sigmaR2_early/2 * bias_ramp[y])) * sexratio[s] # early period recruitment
        if(y >= sigmaR_switch && y < n_yrs) NAA[r,y,1,s] = R0[r] * exp(ln_RecDevs[r,y] - (sigmaR2_late/2 * bias_ramp[y])) * sexratio[s] # late period recruitment
        if(y == n_yrs) NAA[r,y,1,s] = R0[r] * sexratio[s] # mean recruitment in terminal year
      } # end s loop
      Rec[r,y] = sum(NAA[r,y,1,]) # get annual recruitment container here
    } # end y loop
  } # end r loop
  

  ## Population Projection ---------------------------------------------------
  for(y in 1:n_yrs) {
    
    # Recruits don't move
    if(do_recruits_move == 0) {
      # Apply movement after ageing processes - start movement at age 2
      for(a in 2:n_ages) for(s in 1:n_sexes) NAA[,y,a,s] = t(NAA[,y,a,s]) %*% Movement[,,y,a,s]
      for(r in 1:n_regions) NAA[r,y,1,] = Rec[r,y] * sexratio
    } # end if recruits don't move
    # Recruits move here
    if(do_recruits_move == 1) for(a in 1:n_ages) for(s in 1:n_sexes) NAA[,y,a,s] = t(NAA[,y,a,s]) %*% Movement[,,y,a,s]
    
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
  } # end y loop

  # Get total biomass and SSB here
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      Total_Biom[r,y] = sum(as.vector(NAA[r,y,,]) * as.vector(WAA[r,y,,])) # Total Biomass
      SSB[r,y] = sum(as.vector(NAA[r,y,,1]) * as.vector(WAA[r,y,,1]) * MatAA[r,y,,1]) # Spawning Stock Biomass 
      if(n_sexes == 1) SSB[r,y] = SSB[r,y] * 0.5 # If single sex model, multiply SSB calculations by 0.5 
    } # end y loop
  } # end r loop
  

  ## Fishery Observation Model -----------------------------------------------
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
  

  ## Tagging Observation Model -----------------------------------------------
  if(UseTagging == 1) {
    for(tc in 1:n_tag_cohorts) {
      tr = tag_release_indicator[tc,1] # extract tag release region
      ty = tag_release_indicator[tc,2] # extract tag release year
      
      for(ry in 1:min(max_tag_liberty, n_yrs - ty + 1)) {
        # Set up variables for tagging dynamics
        y = ty + ry - 1 # Get index for actual year in the model (instead of tag year)
        # get fishing mortality estimates (assumes uniform selex) 
        if(n_fish_fleets == 1) tmp_F = Fmort[,y,] # only a single fishing fleet
        if(n_fish_fleets > 1) tmp_F = rowSums(Fmort[,y,]) # multiple fleets
        
        # Get total mortality and tag shedding (discount if not tagging start of the year for the first recapture year)
        if(ry == 1 && t_tagging != 0) tmp_Z = (natmort[,y,,,drop = FALSE] + tmp_F + exp(ln_Tag_Shed)) * t_tagging # discounting if ry == 1
        else tmp_Z = (natmort[,y,,,drop = FALSE] + tmp_F + exp(ln_Tag_Shed))
        
        # Run tagging dynamics
        if(ry == 1) Tags_Avail[1,tc,tr,,] = Tagged_Fish[tc,,] * exp(-exp(ln_Init_Tag_Mort)) # Tag induced mortality in the first recapture year

        # Move tagged fish around after mortality and ageing (movement only occurs after first release year if tagging occurs mid year - beginning of the year process)
        if(t_tagging != 0) if(ry > 1) for(a in 1:n_ages) for(s in 1:n_sexes) Tags_Avail[ry,tc,,a,s] = t(Tags_Avail[ry,tc,,a,s]) %*% Movement[,,y,a,s]
        else for(a in 1:n_ages) for(s in 1:n_sexes) Tags_Avail[ry,tc,,a,s] = t(Tags_Avail[ry,tc,,a,s]) %*% Movement[,,y,a,s]
        
        # Mortality and ageing of tagged fish
        for(a in 1:n_ages) {
          for(s in 1:n_sexes) {
            if(a < n_ages) Tags_Avail[ry+1,tc,,a+1,s] = Tags_Avail[ry,tc,,a,s] * exp(-tmp_Z[,1,a,s]) # if not plus group
            else Tags_Avail[ry+1,tc,,n_ages,s] = Tags_Avail[ry+1,tc,,n_ages,s] + (Tags_Avail[ry,tc,,n_ages,s] * exp(-tmp_Z[,1,n_ages,s])) # accumulate plus group
          } # end s loop
        } # end a loop
        
        # Get predicted recaptures (Baranov's)
        Tag_Reporting[,y] = plogis(Tag_Reporting_Pars[,y]) # Inverse logit transform tag reporting rate parameters
        Pred_Tag_Recap[ry,tc,,,] = Tag_Reporting[,y] * (tmp_F / tmp_Z[,1,,]) * Tags_Avail[ry,tc,,,] * (1 - exp(-tmp_Z[,1,,])) 
        
      } # end ry loop
    } # end tc loop
  } # end if for using tagging data 

  
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
            FishIdx_nLL[r,y,f] = UseFishIdx[r,y,f] * (log(ObsFishIdx[r,y,f] + 1e-4) - log(PredFishIdx[r,y,f] + 1e-4))^2 /
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
          # Expected and Observed values
          Exp = CAA[,y,,,f], Obs = ObsFishAgeComps[,y,,,f], 
          # Input sample size and multinomial weight
          ISS = ISS_FishAgeComps[,y,,f], Wt_Mltnml = Wt_FishAgeComps[,,f], 
          # Composition and Likelihood Type
          Comp_Type = FishAgeComps_Type[y,f], Likelihood_Type = FishAgeComps_LikeType[f], 
          # overdispersion par, Number of sexes, regions, age or length comps, and ageing error
          ln_theta = ln_FishAge_theta[,f], n_regions =  n_regions, n_sexes = n_sexes, age_or_len = 0, AgeingError = AgeingError, use = UseFishAgeComps[,y,f], n_bins = n_ages)
      } # if we have fishery age comps

      # Fishery Length Compositions
      if(sum(UseFishLenComps[,y,f]) >= 1 && fit_lengths == 0) {
        FishLenComps_nLL[,y,,f] = Get_Comp_Likelihoods(
          # Expected and Observed values
          Exp = CAL[,y,,,f], Obs = ObsFishLenComps[,y,,,f], 
          # Input sample size and multinomial weight
          ISS = ISS_FishLenComps[,y,,f], Wt_Mltnml = Wt_FishLenComps[r,,f], 
          # Composition and Likelihood Type
          Comp_Type = FishLenComps_Type[y,f], Likelihood_Type = FishLenComps_LikeType[f], 
          # overdispersion, Number of sexes, regions age or length comps, and ageing error
          ln_theta = ln_FishLen_theta[,f], n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1, AgeingError = NA, use = UseFishLenComps[,y,f], n_bins = n_lens) 
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
              SrvIdx_nLL[r,y,sf] = UseSrvIdx[r,y,sf] * (log(ObsSrvIdx[r,y,sf] + 1e-4) - log(PredSrvIdx[r,y,sf] + 1e-4))^2 /
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
            # Expected and Observed values
            Exp = SrvIAA[,y,,,sf], Obs = ObsSrvAgeComps[,y,,,sf], 
            # Input sample size and multinomial weight
            ISS = ISS_SrvAgeComps[,y,,sf], Wt_Mltnml = Wt_SrvAgeComps[,,sf],
            # Composition and Likelihood Type
            Comp_Type = SrvAgeComps_Type[y,sf], Likelihood_Type = SrvAgeComps_LikeType[sf], 
            # overdispersion, Number of sexes, regions, age or length comps, and ageing error
            ln_theta = ln_SrvAge_theta[,sf], n_regions = n_regions, n_sexes = n_sexes, age_or_len = 0, AgeingError = AgeingError, use = UseSrvAgeComps[,y,sf], n_bins = n_ages)
        } # if we have survey age comps

        # Survey Length Compositions
        if(sum(UseSrvLenComps[,y,sf]) >= 1 && fit_lengths == 0) {
          SrvLenComps_nLL[,y,,sf] = Get_Comp_Likelihoods(
            # Expected and Observed values
            Exp = SrvIAL[,y,,,sf], Obs = ObsSrvLenComps[,y,,,sf], 
            # Input sample size and multinomial weight
            ISS = ISS_SrvLenComps[,y,,sf], Wt_Mltnml = Wt_SrvLenComps[,,sf], 
            # Composition and Likelihood Type
            Comp_Type = SrvLenComps_Type[y,sf], Likelihood_Type = SrvLenComps_LikeType[sf], 
            # overdispersion, Number of sexes, regions, age or length comps, and ageing error
            ln_theta = ln_SrvLen_theta[,sf], n_regions = n_regions, n_sexes = n_sexes, age_or_len = 1, AgeingError = NA, use = UseSrvLenComps[,y,sf], n_bins = n_lens)
        } # if we have survey length comps

      } # end sf loop
    } # end y loop
  

  ## Tag Likelihoods ---------------------------------------------------------
  if(UseTagging == 1) {
    for(tc in 1:n_tag_cohorts) {
      
      # set up tagging cohort indexing 
      tr = tag_release_indicator[tc,1] # extract tag release region
      ty = tag_release_indicator[tc,2] # extract tag release year
      
      for(ry in mixing_period:min(max_tag_liberty, n_yrs - ty + 1)) { # loop through recapture years
        for(r in 1:n_regions) {
          for(a in 1:n_ages) {
            for(s in 1:n_sexes) {
              
              # Poisson likelihood
              if(Tag_LikeType == 0) Tag_nLL[ry,tc,r,a,s] = -dpois(Obs_Tag_Recap[ry,tc,r,a,s] + 1e-10, 
                                                                  Pred_Tag_Recap[ry,tc,r,a,s] + 1e-10, log = TRUE) 
              
              # Negative binomial likelihood
              if(Tag_LikeType == 1) {
                # set up robust negative binomial
                log_mu = log(Pred_Tag_Recap[ry,tc,r,a,s] + 1e-10) # log mean
                log_var_minus_mu = 2 * log_mu - ln_tag_theta # log_var_minus_mu
                Tag_nLL[ry,tc,r,a,s] = -dnbinom_robust(Obs_Tag_Recap[ry,tc,r,a,s] + 1e-10, log_mu = log_mu, 
                                                       log_var_minus_mu = log_var_minus_mu, log = TRUE) 
              } # end if for negative binomial likelihood
              
            } # end s loop
          } # end a loop
        } # end r loop
        
        # Multinomial likelihood (release conditioned)
        if(Tag_LikeType == 2) {
          # Set up predicted inputs
          tmp_n_tags_released = sum(Tagged_Fish[tc,,] + 1e-10) # number of tags released for a given tag cohort
          tmp_pred_c = (Pred_Tag_Recap[ry,tc,,,] + 1e-10) / tmp_n_tags_released # recaptured
          tmp_pred = c(tmp_pred_c, 1 - sum(tmp_pred_c)) # combine recaptured and non-recaptured
          # Set up observed inputs
          tmp_obs_c = (Obs_Tag_Recap[ry,tc,,,] + 1e-10) / tmp_n_tags_released # recaptured
          tmp_obs = c(tmp_obs_c, 1 - sum(tmp_obs_c)) # combine recaptured and non-recaptured
          Tag_nLL[ry,tc,1,1,1] = -tmp_n_tags_released * sum((tmp_obs) * log(tmp_pred))
        } # end if multinomial release conditioned
        
        # Multinomial likelihood (recapture conditioned)
        if(Tag_LikeType == 3) {
          # set up inputs
          tmp_n_tags_recap = sum(Obs_Tag_Recap[ry,tc,,,] + 1e-10) # number of recaptures
          tmp_obs = (Obs_Tag_Recap[ry,tc,,,] + 1e-10) / tmp_n_tags_recap # get observed probabilities
          tmp_pred = (Pred_Tag_Recap[ry,tc,,,] + 1e-10) / sum(Pred_Tag_Recap[ry,tc,,,] + 1e-10) # get predicted recapture probabilities
          Tag_nLL[ry,tc,1,1,1] = -1 * tmp_n_tags_recap * sum(((tmp_obs) * log(tmp_pred))) # recapture likelihood
        } # end if multinomial recapture conditioned
      } # end ry loop

    } # end tc loop
  } # if we are using tagging data

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
      for(r in 1:n_regions) {
        # iid selectivity
        if(cont_tv_fish_sel[r,f] == 1) {
          for(s in 1:n_sexes) {
            for(y in 1:n_yrs) {
              sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev1[r,y,s,f], 0, exp(ln_fishsel_dev1_sd[r,s,f]), TRUE)
              sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev2[r,y,s,f], 0, exp(ln_fishsel_dev2_sd[r,s,f]), TRUE)
            } # end y loop
          } # end s loop
        } # end if for iid selectivity

        # random walk selectivity
        if(cont_tv_fish_sel[r,f] == 2) {
          for(s in 1:n_sexes) {
            sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev1[r,1,s,f], 0, 50, TRUE) # initialize first value w/ large prior (prior sd is just arbitrarily large)
            sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev2[r,1,s,f], 0, 50, TRUE) # initialize first value w/ large prior (prior sd is just arbitrarily large)
            for(y in 2:n_yrs) {
              sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev1[r,y,s,f], ln_fishsel_dev1[r,y-1,s,f], exp(ln_fishsel_dev1_sd[r,s,f]), TRUE)
              sel_Pen = sel_Pen + -dnorm(ln_fishsel_dev2[r,y,s,f], ln_fishsel_dev2[r,y-1,s,f], exp(ln_fishsel_dev2_sd[r,s,f]), TRUE)
            } # end y loop
          } # end s loop
        } # end if for random walk selectivity

      } # end r loop
    } # end f loop

  
    ### Recruitment (Penalty) ----------------------------------------------------
    if(likelihoods == 0) {
      for(r in 1:n_regions) {
        Init_Rec_nLL[r,] = (ln_InitDevs[r,] / exp(ln_sigmaR_early))^2 # initial age structure penalty
        for(y in 1:(sigmaR_switch-1)) {
          Rec_nLL[r,y] = (ln_RecDevs[r,y]/exp(ln_sigmaR_early))^2 + bias_ramp[y]*ln_sigmaR_early # early period
        } # end first y loop
        for(y in (sigmaR_switch:(n_yrs-1))) {
          Rec_nLL[r,y] = (ln_RecDevs[r,y]/exp(ln_sigmaR_late))^2 + bias_ramp[y]*ln_sigmaR_late # late period
        } # end second y loop
      } # end r loop
      Rec_nLL = 0.5 * sum(Rec_nLL)  # multiply by 0.5 and sum
      Init_Rec_nLL = 0.5 * sum(Init_Rec_nLL) # multiply by 0.5 and sum
    }
  
    if(likelihoods == 1) {
      for(r in 1:n_regions) {
        Init_Rec_nLL[r,] = -dnorm(ln_InitDevs[r,], 0, exp(ln_sigmaR_early), TRUE) # initial age structure penalty
        for(y in 1:(sigmaR_switch-1)) {
          Rec_nLL[r,y] = -dnorm(ln_RecDevs[r,y], 0, exp(ln_sigmaR_early), TRUE)
        } # first y loop
        # Note that this penalizes the terminal year rec devs, which is estimated in this case
        for(y in sigmaR_switch:(n_yrs)) {
          Rec_nLL[r,y] = -dnorm(ln_RecDevs[r,y], 0, exp(ln_sigmaR_late), TRUE)
        } # end second y loop
      }
    } 

    # Apply likelihood weights here and compute joint negative log likelihood
    jnLL = (Wt_Catch * sum(Catch_nLL)) + # Catch likelihoods
           (Wt_FishIdx * sum(FishIdx_nLL)) + # Fishery Index likelihood
           (Wt_SrvIdx * sum(SrvIdx_nLL)) + # Survey Index likelihood
           sum(FishAgeComps_nLL) + # Fishery Age likelihood
           sum(FishLenComps_nLL) + # Fishery Length likelihood
           sum(SrvAgeComps_nLL) + # Survey Age likelihood
           sum(SrvLenComps_nLL) + # Survey Length likelihood
           sum(Tag_nLL) + # Tagging likelihood
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
  RTMB::REPORT(init_F)
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
  
  # Tagging Processes
  RTMB::REPORT(Pred_Tag_Recap)
  RTMB::REPORT(Tags_Avail)
  RTMB::REPORT(Tag_Reporting)
  
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
  RTMB::REPORT(Tag_nLL)
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
