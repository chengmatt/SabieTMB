  # Purpose: To do model development for a spatial model 
  # Creator: Matthew LH. Cheng
  # Date Created: 12/24/24
  
  # Set up ------------------------------------------------------------------
  
  library(here)
  library(R2admb)
  library(tidyverse)
  library(RTMB)
  
  source(here("R", "functions", "Francis_Reweight.R"))
  source(here("R", "model", "v3", "Sabie_RTMB_v3.R"))
  source(here("R", "functions", "Utility_Functions.R"))
  
  sim_out = readRDS(here("sim_out.RDS"))
  
  n_sims <- 100
  ssb_mat <- array(NA, dim = c(length(1:20), 2, n_sims))
  r0_mat <- matrix(NA, nrow = n_sims, ncol = 2)
  tagrep <- vector()
  status <- vector()
  dir <- vector()
  pd = vector()
  grad = vector()
  catch_type <- 1 # aggregated in early period
  est_all_regional_F <- 1 # don't estimate all reigonal F (0) (i.e., some are aggregated) (1) estimate all 
  comp_joint = 1
  # catch type = 0, est_all_reigonal_F = 0, aggregated likelihood initially, followed by aggregated F mean and F devs, and then regional after aggregated period
  # catch type = 0, estimate_all_regional_F = 1, aggregated likelihood initially followed by estimating regional F mean and F devs in all periods
  # catch type = 1, estimate_all_regional_F = 1, regional likelihoods and regional f mean and devs
  move_cv_mat = matrix(NA, nrow = 20, ncol = n_sims)
  
  # Testing with simulated data ---------------------------------------------
  for(sim in 1:n_sims) {
    # Prepare Data ------------------------------------------------------------
    data <- list() # make data list
    
    # Set up dimensions
    data$n_regions <- 2 # number of regions
    data$ages <- 1:15 # ages
    data$lens <- 1 # lengths
    data$years <- as.numeric(1:20) # years
    data$n_sexes <- 1 # number of sexes
    data$n_fish_fleets <- 1 # number of fishery fleets (0 == fixed gear, 1 == trawl gear)
    data$n_srv_fleets <- 1 # number of survey fleets (0 == domestic ll survey, 1 == domestic trawl survey, 2 == coop jp ll survey)
    
    # Recruitment stuff
    data$init_F_prop <- 0 # initial F proportion for initializing population
    data$do_rec_bias_ramp <- 0 # do bias ramp (slot 0 == don't do bias ramp, 1 == do bias ramp)
    data$bias_year <- NA
    data$sigmaR_switch <- 0
    data$sexratio <- as.vector(1) # recruitment sex ratio (assuming 1, since single sex)
    
    # Movement stuff
    data$do_recruits_move = 0 # recruits dont move
    data$use_fixed_movement = 0 # use fixed movement
    data$Fixed_Movement = array(sim_out$movement_matrix[,,,,,sim,drop=FALSE], dim = c(data$n_regions,
                                                                                      data$n_regions,
                                                                                      length(data$years),
                                                                                      length(data$ages),
                                                                                      data$n_sexes))
    
    # data$Fixed_Movement = array(diag(1,2), dim = c(data$n_regions,data$n_regions, length(data$years),length(data$ages),data$n_sexes))
    
    
    data$likelihoods <- 1
    data$Wt_Catch <- 1 # Catch weights
    data$Wt_FishIdx <- 1 # fishery index weights
    data$Wt_SrvIdx <- 1 # survey index weights
    data$Wt_Rec <- 1 # recruitment weights
    data$Wt_F <- 1 # fishing mortality penalty weights
    
    # Biological Processes
    # Natural Mortality
    data$Use_M_prior <- 0 # use natural mortality prior
    data$M_prior <- c(0.3, 0.1) # Mean and CV for M prior
    
    # Weight at age
    data$WAA <- array(aperm(sim_out$WAA[,,,,sim, drop = FALSE], perm = c(2,1,3,4,5)), dim = c(data$n_regions, 
                                                                                              length(data$years),
                                                                                              length(data$ages),
                                                                                              data$n_sexes))
    
    # Maturity at age
    data$MatAA <- array(aperm(sim_out$Maturity_AA[,,,,sim, drop = FALSE], perm = c(2,1,3,4,5)), dim = c(data$n_regions, 
                                                                                                length(data$years),
                                                                                                length(data$ages),
                                                                                                data$n_sexes))
    
    # Ageing error
    data$AgeingError <- diag(1, length(data$ages)) # ageing error matrix
    
    # Size Transition Matrix (Growth)
    data$fit_lengths <- 1
    data$SizeAgeTrans <- NA
    
    ## Observations  -----------------------------------------------------------
    ### Fishery Observations ----------------------------------------------------
    
    # Catches
    if(catch_type == 0) {
      data$ObsCatch <- array(0, dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
      data$ObsCatch[1,1:7,1] <- apply(aperm(sim_out$Obs_Catch[1:7,,,sim, drop = FALSE], perm = c(2,1,3,4)), 2:4, sum)
      data$ObsCatch[,-c(1:7),1] <- aperm(sim_out$Obs_Catch[-c(1:7),,,sim, drop = FALSE], perm = c(2,1,3,4))  
      data$Catch_Type <- array(c(rep(0, 7), rep(1, 13)), dim = c(length(data$years), data$n_fish_fleets))
      data$UseCatch <- array(1, c(data$n_regions, length(data$years), data$n_fish_fleets))
      if(est_all_regional_F == 0) data$est_all_regional_F = 0
      if(est_all_regional_F == 1) data$est_all_regional_F = 1
    } # aggregated in early period
    
    if(catch_type == 1) {
      data$ObsCatch <- array(aperm(sim_out$Obs_Catch[,,,sim, drop = FALSE], perm = c(2,1,3,4)), dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
      data$ObsCatch[data$ObsCatch == 0] <- NA # set 0 catches to NA so we aren't fitting
      data$UseCatch <- array(1, c(data$n_regions, length(data$years), data$n_fish_fleets))
      data$Catch_Type <- array(rep(1,20), dim = c(length(data$years), data$n_fish_fleets))
      data$est_all_regional_F = 1
    }
    
    data$Catch_Constant <- rep(0, data$n_fish_fleets)
    
    # Fishery Indices
    data$ObsFishIdx <- array(NA, c(data$n_regions, length(data$years), data$n_fish_fleets))
    data$ObsFishIdx_SE <- array(NA, c(data$n_regions, length(data$years), data$n_fish_fleets))
    data$UseFishIdx <- array(0, c(data$n_regions, length(data$years), data$n_fish_fleets))
    colnames(data$ObsFishIdx) <- data$years # define row years
    colnames(data$ObsFishIdx_SE) <- data$years # define row years
    
    # fishery CPUE
    data$share_sel <- 0
    data$UseFishIdx[is.na(data$ObsFishIdx)] <- 0 # don't fit if missing data
    
    # Fishery Age Comps
    data$ObsFishAgeComps <- array(aperm(sim_out$Obs_FishAgeComps[,,,,,sim,drop = FALSE], perm = c(2,1,3,4,5,6)), dim = c(data$n_regions, length(data$years), length(data$ages), data$n_sexes, data$n_fish_fleets))

    colnames(data$ObsFishAgeComps) <- data$years # define row years
    data$UseFishAgeComps <- array(1, dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
    colnames(data$UseFishAgeComps) <- data$years # define row years
    # data$UseFishAgeComps[,1:7,] <- 0
    
    # Data weighting for fishery age compositions
    data$ISS_FishAgeComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
    colnames(data$ISS_FishAgeComps) <- data$years # define row years
    if(comp_joint == 1) data$ISS_FishAgeComps[,,1,1] <- apply(data$ObsFishAgeComps, c(2),sum) else {
      for(r in 1:data$n_regions)  data$ISS_FishAgeComps[r,,,1] <- apply(data$ObsFishAgeComps, c(1,2),sum)[r,]
    }
    
    data$Wt_FishAgeComps <- array(NA, dim = c(data$n_regions, data$n_sexes, data$n_fish_fleets)) # weights for fishery age comps
    data$Wt_FishAgeComps[,1,1] <- 1 # Weight for fixed gear age comps
    
    # Fishery Length Comps
    data$ObsFishLenComps <- array(NA, dim = c(data$n_regions, length(data$years), length(data$lens), data$n_sexes, data$n_fish_fleets))
    colnames(data$ObsFishLenComps) <- data$years # define row years
    data$UseFishLenComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
    colnames(data$UseFishLenComps) <- data$years # define row years
    data$ISS_FishLenComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
    colnames(data$ISS_FishLenComps) <- data$years # define row years
    
    # Composition munging stuff
    data$FishAgeComps_LikeType <- array(0, dim = c(data$n_fish_fleets)) # multinomial for both fleets
    data$FishLenComps_LikeType <- array(0, dim = c(data$n_fish_fleets)) # multinomial for both fleets
    if(comp_joint == 1) data$FishAgeComps_Type <- array(3, dim = c(length(data$years), data$n_fish_fleets)) else data$FishAgeComps_Type <- array(2, dim = c(length(data$years), data$n_fish_fleets)) # Joint age comps # FLAG need to figure out what to do when comps are aggregated by region and sex

    if(catch_type == 0) data$FishAgeComps_Type[1:7,] <- 0
    data$FishLenComps_Type <- array(3, dim = c(length(data$years), data$n_fish_fleets)) # Joint length comps
    
    ### Survey Observations -----------------------------------------------------
    # Survey Indices
    data$ObsSrvIdx <- array(aperm(sim_out$Obs_SrvIdx[,,,sim,drop = FALSE], perm = c(2,1,3,4)), dim = c(data$n_regions, length(data$years), data$n_srv_fleets))
    data$ObsSrvIdx_SE <- (data$ObsSrvIdx * 0.3)
    data$UseSrvIdx <- array(1, dim = c(data$n_regions, length(data$years), data$n_srv_fleets))
    colnames(data$ObsSrvIdx) <- data$years # define row years
    colnames(data$ObsSrvIdx_SE) <- data$years # define row years
    
    # Survey Age Comps
    data$ObsSrvAgeComps <- array(aperm(sim_out$Obs_SrvAgeComps[,,,,,sim,drop=FALSE], perm = c(2,1,3,4,5,6)), dim = c(data$n_regions, length(data$years), length(data$ages), data$n_sexes, data$n_srv_fleets))
    colnames(data$ObsSrvAgeComps) <- data$years # define row years
    data$UseSrvAgeComps <- array(1, dim = c(data$n_regions, length(data$years), data$n_srv_fleets))
    
    # Data weighting for survey age compositions
    data$ISS_SrvAgeComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_srv_fleets))
    colnames(data$ISS_SrvAgeComps) <- data$years # define row years
    
    if(comp_joint == 1) data$ISS_SrvAgeComps[,,1,1] <- apply(data$ObsSrvAgeComps, c(2),sum) else {
      for(r in 1:data$n_regions)  data$ISS_SrvAgeComps[r,,,1] <- apply(data$ObsSrvAgeComps, c(1,2),sum)[r,]
    }
    
    data$Wt_SrvAgeComps <- array(NA, dim = c(data$n_regions, data$n_sexes, data$n_srv_fleets)) # weights for survey age comps
    data$Wt_SrvAgeComps[,1,1] <- 1 # Weight for domestic survey ll gear age comps
    
    # Survey Length Comps
    data$ObsSrvLenComps <- array(NA, dim = c(data$n_regions,length(data$years), length(data$lens), data$n_sexes, data$n_srv_fleets))
    colnames(data$ObsSrvLenComps) <- data$years # define row years
    data$UseSrvLenComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_srv_fleets))
    colnames(data$UseSrvLenComps) <- data$years # define row years
    
    # Data weighting for survey length compositions
    data$ISS_SrvLenComps <- array(0, dim = c(data$n_regions,length(data$years), data$n_sexes, data$n_srv_fleets))
    data$Wt_SrvLenComps <- array(0, dim = c(data$n_regions, data$n_sexes, data$n_srv_fleets)) # weights for survey age comps
    
    # Composition munging stuff
    data$SrvAgeComps_LikeType <- array(0, dim = c(data$n_srv_fleets)) # multinomial for both survey fleet
    data$SrvLenComps_LikeType <- array(0, dim = c(data$n_srv_fleets)) # multinomial for both survey fleet
    data$SrvLenComps_Type <- array(3, dim = c(length(data$years), data$n_srv_fleets)) # split for both survey fleet length
    if(comp_joint == 1) data$SrvAgeComps_Type <- array(3, dim = c(length(data$years), data$n_srv_fleets)) else data$SrvAgeComps_Type <- array(2, dim = c(length(data$years), data$n_srv_fleets)) # split for both survey age length
    
    ### Fishery Stuff -----------------------------------------------------
    # Selectivity
    data$cont_tv_fish_sel <- array(0, dim = c(data$n_regions, data$n_fish_fleets)) # no timevarying selex continously
    # Time Block Specification
    data$fish_sel_blocks <- array(NA, dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
    data$fish_sel_blocks[] <- 0 # block one fishery ll selex
    # Selectivity Model
    data$fish_sel_model <- array(NA, dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
    data$fish_sel_model[] <- 0 # Logistic selectivity for fixed gear
    # Catchability
    # Time Block Specification
    data$fish_q_blocks <- data$fish_sel_blocks # catchability blocks are same as selectivity blocks
    data$fish_idx_type <- array(1, dim = c(data$n_regions, data$n_fish_fleets)) # fishery index type (0 == abundance, 1 == biomass)
    
    ### Survey Selectivity ------------------------------------------------------
    # Time Block Specification
    data$srv_sel_blocks <- array(NA, dim = c(data$n_regions, length(data$years), data$n_srv_fleets))
    data$srv_sel_blocks[] <- 0 # block one survey ll selex
    # Selectivity Model
    data$srv_sel_model <- array(NA, dim = c(data$n_regions, length(data$years), data$n_srv_fleets))
    data$srv_sel_model[] <- 0 # Logistic selectivity for longline survey
    
    # Time Block Specification
    data$srv_q_blocks <- array(0, dim = c(data$n_regions, length(data$years), data$n_srv_fleets)) # catchability blocks are same as selectivity blocks
    data$srv_idx_type <- array(0, dim = c(data$n_regions, data$n_srv_fleets))
    

    # Tagging Stuff -----------------------------------------------------------
    data$UseTagging = 1
    data$tag_release_indicator <- sim_out$Tag_Release_Ind # tag release indicator, with rows = cohorts, first col = release regino, second col = release year
    data$n_tag_cohorts <- nrow(data$tag_release_indicator) # number of tag cohorts
    data$max_tag_liberty <- 30 # maximum liberty to track cohorts
    data$Tagged_Fish <- array(sim_out$Tag_Fish[,,,sim], dim = c(data$n_tag_cohorts, length(data$ages), data$n_sexes)) # tagged fish
    data$Obs_Tag_Recap <- array(sim_out$Obs_Tag_Recap[,,,,,sim], dim = c(data$max_tag_liberty, data$n_tag_cohorts, data$n_regions, length(data$ages), data$n_sexes))
    data$Tag_LikeType <- 1 # poisson likelihood
    data$mixing_period <- 2 # when to start mixing period after first release year
    data$t_tagging <- 0.5 # discounting for tagging 
    
    # Prepare Parameters ------------------------------------------------------
    parameters <- list()
    parameters$dummy <- 1
    
    ### Fishery Stuff ---------------------------------------------------------
    
    # Fishing Mortality
    parameters$ln_sigmaC <- array(log(1e-3), dim = c(data$n_regions, data$n_fish_fleets))
    parameters$ln_F_mean <- array(0, dim = c(data$n_regions, data$n_fish_fleets)) # mean fishing mortality
    parameters$ln_F_devs <- array(0, dim = c(data$n_regions, length(data$years), data$n_fish_fleets)) # fishing mortality deviations from mean
    
    parameters$ln_F_devs_AggCatch = array(0, dim = c(sum(data$Catch_Type == 0), data$n_fish_fleets))
    parameters$ln_F_mean_AggCatch = log(0.1)
    
    # Set up continuous fishery selectivity stuff
    parameters$ln_fishsel_dev1 <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
    parameters$ln_fishsel_dev2 <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
    parameters$ln_fishsel_dev1_sd <- array(0.1, dim = c(data$n_regions, data$n_sexes, data$n_fish_fleets))
    parameters$ln_fishsel_dev2_sd <- array(0.1, dim = c(data$n_regions, data$n_sexes, data$n_fish_fleets))
    
    # Fixed Gear Fishery three time blocks
    max_fish_blks <- 1 # maximum number of fishery blocks for any fleet
    max_fish_pars <- 2 # maximum number of fishery fixed parameters for any fleet
    parameters$ln_fish_fixed_sel_pars <- array(0, dim = c(data$n_regions, max_fish_pars, max_fish_blks, data$n_sexes, data$n_fish_fleets))
    parameters$ln_fish_fixed_sel_pars[,1,,,] = log(3.333333)
    parameters$ln_fish_fixed_sel_pars[,2,,,] = log(3)
    
    # Fishery Catchability
    # Fixed Gear Fishery Catchability
    parameters$ln_fish_q <- array(0, dim = c(data$n_regions, max_fish_blks, data$n_fish_fleets))
    
    ### Survey Stuff -----------------------------------------------------
    # Survey Selectivity
    max_srv_blks <- 1 # maximum number of survey blocks for any fleet
    max_srv_pars <- 2 # maximum number of survey fixed parameters for any fleet
    parameters$ln_srv_fixed_sel_pars <- array(0, dim = c(data$n_regions, max_srv_pars, max_srv_blks, data$n_sexes, data$n_srv_fleets))
    parameters$ln_srv_fixed_sel_pars[,1,,,] = log(2)
    parameters$ln_srv_fixed_sel_pars[,2,,,] = log(1)
    
    # Survey Catchability
    max_q_srv_blks <- 1 # maximum catchability blocks
    # Longline Survey Catchability
    parameters$ln_srv_q <- array(log(1), dim = c(data$n_regions,max_q_srv_blks, data$n_srv_fleets))
    
    ### Natural Mortality -------------------------------------------------------
    parameters$ln_M <- log(0.15) # Female M (base)
    parameters$M_offset <- 0 # Male M offset (accidently jittered from OG assessment)
    
    # Recruitment -------------------------------------------------------------
    parameters$ln_global_R0 <- log(1e7 + 5e6) # mean recruitment
    parameters$R0_prop <- array(c(0.3, 0.5, 0.3), dim = c(data$n_regions - 1))
    parameters$ln_InitDevs <- array(0, dim = c(data$n_regions, length(data$ages) - 2))
    parameters$ln_RecDevs <- array(0, dim = c(data$n_regions, length(data$years)))
    parameters$ln_sigmaR_early <- log(0.5) # early sigma R
    parameters$ln_sigmaR_late <- log(0.5)  # late sigma R
    
    # Comp Likelihood Stuff ---------------------------------------------------
    parameters$ln_FishAge_theta <- array(log(0.1), dim = c(data$n_regions, data$n_fish_fleets))
    parameters$ln_FishLen_theta <- array(0, dim = c(data$n_regions, data$n_fish_fleets))
    parameters$ln_SrvAge_theta <- array(0, dim = c(data$n_regions, data$n_srv_fleets))
    parameters$ln_SrvLen_theta <- array(0, dim = c(data$n_regions, data$n_srv_fleets))
    
    # Movement Stuff ---------------------------------------------------
    parameters$move_pars <- array(0, dim = c(data$n_regions, data$n_regions - 1, length(data$years), length(data$ages), data$n_sexes))

    # Tagging Stuff -----------------------------------------------------------
    parameters$ln_Init_Tag_Mort <- log(1e-10) # initial tag induced mortality
    parameters$ln_Tag_Shed <- log(1e-10) # annual tag sheeding
    parameters$Tag_Reporting_Pars <- array(log(0.2 / (1 - 0.2)), dim = c(data$n_regions, length(data$years)))
    parameters$ln_tag_theta <- log(0.1)
    
    # Mapping -----------------------------------------------------------------
    mapping <- list()
    mapping$dummy <- factor(NA)
    mapping$ln_sigmaR_late <- factor(NA)
    mapping$ln_sigmaR_early <- factor(NA) # fix early sigma R
    mapping$M_offset <- factor(NA) # fix natural mortality offset
    mapping$ln_fish_q <- factor(rep(NA, length(parameters$ln_fish_q)))
    
    # Fixing sigmas for fishery catch and Fdevs here
    mapping$ln_sigmaC <- factor(rep(NA, length(parameters$ln_sigmaC)))
    
    # Fixing sel pars
    # mapping$ln_fish_fixed_sel_pars = factor(c(1,1,1,1,2,2,2,2))
    # mapping$ln_srv_fixed_sel_pars = factor(c(1,1,1,1,2,2,2,2))
    mapping$ln_fish_fixed_sel_pars = factor(c(1,1,2,2))
    mapping$ln_srv_fixed_sel_pars = factor(c(1,1,2,2))
    
    # Fixing continuous time-varying selecitvity stuff
    mapping$ln_fishsel_dev1 <- factor(rep(NA, length(parameters$ln_fishsel_dev1)))
    mapping$ln_fishsel_dev2 <- factor(rep(NA, length(parameters$ln_fishsel_dev2)))
    mapping$ln_fishsel_dev1_sd <- factor(rep(NA, length(parameters$ln_fishsel_dev1_sd)))
    mapping$ln_fishsel_dev2_sd <- factor(rep(NA, length(parameters$ln_fishsel_dev2_sd)))
    
    # Fixing dirichlet mutlinomial stuff
    mapping$ln_FishAge_theta <- factor(rep(NA, length(parameters$ln_FishAge_theta)))
    mapping$ln_FishLen_theta <- factor(rep(NA, length(parameters$ln_FishLen_theta)))
    mapping$ln_SrvAge_theta <- factor(rep(NA, length(parameters$ln_SrvAge_theta)))
    mapping$ln_SrvLen_theta <- factor(rep(NA, length(parameters$ln_SrvLen_theta)))

    # mapping$ln_FishAge_theta <- factor(c(1,NA)) # joint
    # mapping$ln_FishAge_theta <- factor(c(1,1)) # split 
    # mapping$ln_SrvAge_theta <- factor(c(1,2))

    # Fixing movement stuff
    # parameters$move_pars[,,1,1,1] <- 0.1
    # parameters$move_pars[,,1,2,1] <- 0.2
    # parameters$move_pars[,,1,3,1] <- 0.3
    # 
    #     mapping$move_pars = factor(c(rep(1:3, length.out = data$n_regions * length(data$years) * length(1:3)),
    #                                  rep(4:6, length.out = data$n_regions * length(data$years) * length(4:6)),
    #                                  rep(7:10, length.out = data$n_regions * length(data$years) * length(7:10))))

    mapping$move_pars = factor(rep(1:2, length.out = prod(dim(parameters$move_pars))))
    # mapping$move_pars = factor(rep(1:6, length.out = prod(dim(parameters$move_pars))))
    
    # parameters$move_pars[,,1,1,1] <- 0.1
    # parameters$move_pars[2,,1,1,1] <- 0.2
    # parameters$move_pars[3,,1,1,1] <- 0.3
    # parameters$move_pars[4,,1,1,1] <- 0.4
    
    # mapping$move_pars = factor(rep(NA, length.out = prod(dim(parameters$move_pars))))

    # Fixing survey catchability
    mapping$ln_srv_q <- factor(rep(1, length(parameters$ln_srv_q)))
    # mapping$ln_srv_q <- factor(rep(NA, length(parameters$ln_srv_q)))
    
    # mapping$ln_global_R0 <- factor(NA)
    # mapping$R0_prop <- factor(NA)
    
    # Fixing M
    # mapping$ln_M <- factor(NA)

    # Fixing fishing mortlaity stuff
    # mapping$ln_F_devs = factor(rep(NA, length(parameters$ln_F_devs)))
    # mapping$ln_F_mean = factor(rep(NA, length(parameters$ln_F_mean)))
    # data$Fmort_dat = array(aperm(sim_out$Fmort[,,,1], perm = c(2,1)), dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
    # mapping$ln_fish_fixed_sel_pars = factor(rep(NA, length(parameters$ln_fish_fixed_sel_pars)))
    
    # Fix survey selex
    # mapping$ln_srv_fixed_sel_pars = factor(rep(NA, length(parameters$ln_srv_fixed_sel_pars)))
    
    # fixing tagging stuff
    mapping$ln_Init_Tag_Mort <- factor(NA)
    mapping$ln_Tag_Shed <- factor(NA)
    mapping$Tag_Reporting_Pars <- factor(rep(1, length.out = length(parameters$Tag_Reporting_Pars)))
    # mapping$ln_tag_theta <- factor(NA)
    
    # global density dependence
    map_recdevs = parameters$ln_RecDevs
    map_recdevs[1,] = 1:length(map_recdevs[1,])
    map_recdevs[2,] = map_recdevs[1,]
    # map_recdevs[3,] = map_recdevs[1,]
    # map_recdevs[4,] = map_recdevs[1,]
    mapping$ln_RecDevs = factor(map_recdevs)
    
    map_initdevs = parameters$ln_InitDevs
    map_initdevs[1,] = 1:length(map_initdevs[1,])
    map_initdevs[2,] = map_initdevs[1,]
    # map_initdevs[3,] = map_initdevs[1,]
    # map_initdevs[4,] = map_initdevs[1,]
    mapping$ln_InitDevs = factor(map_initdevs)

    # mapping$ln_RecDevs = factor(rep(NA, length(parameters$ln_RecDevs)))
    # mapping$ln_InitDevs = factor(rep(NA, length(parameters$ln_InitDevs)))
    
    data$srv_q_blocks = data$srv_q_blocks + 1
    data$fish_q_blocks = data$fish_q_blocks + 1
    data$fish_sel_blocks = data$fish_sel_blocks + 1
    data$srv_sel_blocks = data$srv_sel_blocks + 1
    data$bias_year = data$bias_year + 1
    data$sigmaR_switch = data$sigmaR_switch + 1
    
    if(catch_type == 0 && est_all_regional_F == 0) mapping$ln_F_devs = factor(c(rep(NA, 14), 1:26))
    if((catch_type == 1 || catch_type == 0) && est_all_regional_F == 1) {
      mapping$ln_F_devs_AggCatch <- factor(rep(NA, data$n_fish_fleets * sum(data$Catch_Type == 0)))
      mapping$ln_F_mean_AggCatch <- factor(rep(NA, data$n_fish_fleets))
    }
    
    # make AD model function
    sabie_rtmb_model <- RTMB::MakeADFun(sabie_RTMB, parameters = parameters, map = mapping)
    
    # Now, optimize the function
    sabie_optim <- stats::nlminb(sabie_rtmb_model$par, sabie_rtmb_model$fn, sabie_rtmb_model$gr,
                                 control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
    # newton steps
    try_improve <- tryCatch(expr =
                              for(i in 1:3) {
                                g = as.numeric(sabie_rtmb_model$gr(sabie_optim$par))
                                h = optimHess(sabie_optim$par, fn = sabie_rtmb_model$fn, gr = sabie_rtmb_model$gr)
                                sabie_optim$par = sabie_optim$par - solve(h,g)
                                sabie_optim$objective = sabie_rtmb_model$fn(sabie_optim$par)
                              }
                            , error = function(e){e}, warning = function(w){w})
    
    sabie_rtmb_model$optim <- sabie_optim # Save optimized model results
    sabie_rtmb_model$sd_rep <- RTMB::sdreport(sabie_rtmb_model) # Get sd report
    sabie_rtmb_model$rep <- sabie_rtmb_model$report(sabie_rtmb_model$env$last.par.best) # Get report
    
    
    sabie_rtmb_model$rep$Pred_Tag_Recap[1,1,,,]
    
    r0_mat[sim,] <- sabie_rtmb_model$rep$R0
    tagrep[sim] <- mean(sabie_rtmb_model$rep$Tag_Reporting)
    PD = sabie_rtmb_model$sd_rep$pdHess
    gradient = max(sabie_rtmb_model$sd_rep$gradient.fixed)
    pd[sim] = PD
    grad[sim] = gradient
    if(PD == TRUE && gradient < 0.005) status[sim] = TRUE else status[sim] = FALSE
    ssb_mat[,1,sim] <- (sabie_rtmb_model$rep$SSB[1,] - sim_out$SSB[,1,sim]) / sim_out$SSB[,1,sim]
    ssb_mat[,2,sim] <- (sabie_rtmb_model$rep$SSB[2,] - sim_out$SSB[,2,sim]) / sim_out$SSB[,2,sim]
    # dir[sim] <- exp(sabie_rtmb_model$sd_rep$par.fixed[names(sabie_rtmb_model$sd_rep$par.fixed) == 'ln_FishAge_theta'])
    
    # o = cbind(sim_out$Obs_FishAgeComps[,1,,,,sim], sim_out$Obs_FishAgeComps[,2,,,,sim])
    # o = (o) / rowSums(o)
    
    # for(i in 1:20) {
    #   tmpy = sim_out$Obs_FishAgeComps[i,,,,,sim,drop = FALSE]
    #   if(any(tmpy == 0)) tmpy = (tmpy + min(tmpy[tmpy != 0]) * 0.01) / sum(tmpy + min(tmpy[tmpy != 0]) * 0.01)
    #   else tmpy = (tmpy) / sum(tmpy)
    #   o1 = array(tmpy, dim = c(data$n_regions,length(data$ages)))
    #   data$ObsFishAgeComps[,i,,,] = o1
    # }
    
    # o <- cbind(data$ObsFishAgeComps[1,,,,], data$ObsFishAgeComps[2,,,,])
    # e = cbind(sabie_rtmb_model$rep$CAA[1,,,,], sabie_rtmb_model$rep$CAA[2,,,,])
    # e = (e) / rowSums(e)
    # tmp_e = log(e[,-length(e)])
    # tmp_e = tmp_e - log(e[,nrow(e)])
    # res = compResidual::reslogistN(t(o), t(tmp_e), diag(dir[sim]^2, 19), do_mult = F)
    # plot(res)
    # res = compResidual::resMulti(t(o), t(e))
    
    # o = cbind(sim_out$Obs_FishAgeComps[,1,,,,sim], sim_out$Obs_FishAgeComps[,2,,,,sim]) * 5
    # o = (o) / rowSums(o)
    # e = cbind(sabie_rtmb_model$rep$CAA[1,,,,], sabie_rtmb_model$rep$CAA[2,,,,])
    # e = (e) / rowSums(e)
    # res = compResidual::resDirM(t(o), t(e) * dir[sim] * 500)
    # 
    # plot(res)
    # sd(res, na.rm = T)

    move_cv = unique(sabie_rtmb_model$sd_rep$sd[names(sabie_rtmb_model$sd_rep$value) == "Movement"] / sabie_rtmb_model$sd_rep$value[names(sabie_rtmb_model$sd_rep$value) == "Movement"])
    move_cv_mat[1:length(move_cv),sim] <- move_cv
    
    par(mfrow = c(2,4))
    #
    plot(sim_out$CAA[20,2,,1,1,sim], type = 'l', col = 'red')
    lines(sabie_rtmb_model$rep$CAA[2,20,,1,1], type = 'l')
    #
    # plot(sim_out$Init_NAA[100,1,,1,sim], type = 'l', col = 'red')
    # lines(sabie_rtmb_model$rep$Init_NAA[100,1,,1], type = 'l')
    #
    # plot(sim_out$NAA[20,2,,1,sim], type = 'l', col = 'red')
    # lines(sabie_rtmb_model$rep$NAA[2,20,,1], type = 'l')
    #
    
    plot(sim_out$Total_Biom[,1,sim], col = 'red', type = 'l', ylab = 'Biomass region 1 (red = simulation, black = est)')
    lines(sabie_rtmb_model$rep$Total_Biom[1,], type = 'l')
    
    plot(sim_out$Total_Biom[,2,sim], col = 'red', type = 'l', ylab = 'Biomass region 2 (red = simulation, black = est)')
    lines(sabie_rtmb_model$rep$Total_Biom[2,], type = 'l')
    #
    # # plot(sabie_rtmb_model$rep$PredCatch[1,,1], type = 'l', ylab = 'Catch region 1 (red = simulation, black = est)')
    # # lines(sim_out$Obss_Catch[,1,,sim], type = 'l', col = 'red')
    #
    # plot(sabie_rtmb_model$rep$PredCatch[2,,1], type = 'l', ylab = 'Catch region 2 (red = simulation, black = est)')
    # lines(sim_out$Obs_Catch[,2,,sim], type = 'l', col = 'red')
    #
    # plot(colSums(sabie_rtmb_model$rep$PredCatch[,,1]), type = 'l', ylab = 'Aggregated catch (red = simulation, black = est)')
    # lines(rowSums(sim_out$Obs_Catch[,,,sim]), type = 'l', col = 'red')
    #
    plot(sabie_rtmb_model$rep$Fmort[1,,1], type = 'l', ylab = 'Fmort Region 1 (red = simulation, black = est)', xlab = 'Year')
    lines(sim_out$Fmort[,1,1,sim], col = 'red')
    #
    # # plot(sabie_rtmb_model$rep$Fmort[2,,1], type = 'l', ylab = 'Fmort Region 2 (red = simulation, black = est)', xlab = 'Year')
    # # lines(sim_out$Fmort[,2,1,sim], col = 'red')
    #
    plot(sim_out$fish_sel[1,1,,1,1,sim], col = 'red')
    lines(sabie_rtmb_model$rep$fish_sel[1,1,,1,1])
    
    # par(mfrow = c(1,3))
    hist((r0_mat[,1] - 5e6) / 5e6, main = round(median((r0_mat[,1] - 5e6) / 5e6, na.rm = T),2), xlab = 'R0 bias region 1')
    hist((r0_mat[,2] - 1e7) / 1e7, main = round(median((r0_mat[,2] - 1e7) / 1e7, na.rm = T),2), xlab = 'R0 bias region 2')
    hist((tagrep - 0.2) / 0.2, main = round(median((tagrep - 0.2) / 0.2, na.rm = T),2), xlab = 'R0 bias region 2')
    
    # hist((r0_mat[,3] - 100) / 100, main = round(median((r0_mat[,3] - 100) / 100, na.rm = T),2), xlab = 'R0 bias region 3')
    # hist((r0_mat[,4] - 800) / 800, main = round(median((r0_mat[,4] - 800) / 800, na.rm = T),2), xlab = 'R0 bias region 4')
  }

  median(move_cv_mat, na.rm = TRUE)
  
  dev.off()
  max(sabie_rtmb_model$sd_rep$gradient.fixed)
  sabie_rtmb_model$sd_rep$par.fixed[which.max(sabie_rtmb_model$sd_rep$gradient.fixed)]

  sabie_rtmb_model$rep$jnLL

  sabie_rtmb_model$rep$R0
  exp(sabie_rtmb_model$sd_rep$par.fixed[names(sabie_rtmb_model$sd_rep$par.fixed) == 'ln_srv_q'])
  exp(sabie_rtmb_model$sd_rep$par.fixed[names(sabie_rtmb_model$sd_rep$par.fixed) == 'ln_M'])

  sabie_rtmb_model$rep$Movement[,,1,1,1]
  movement_matrix[,,1,2,1,1]
