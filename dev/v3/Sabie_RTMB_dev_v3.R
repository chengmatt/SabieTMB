# Purpose: To bridge conduct testing in developing Sablefish RTMB model (v3)
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

# Read in data
tem_dat <- dget(here('2. Base (23.5)_final model', 'tem.rdat'))
ageing_dat <- dget(here('2. Base (23.5)_final model', 'test.rdat'))
tem_admb_dat <- readLines(here("2. Base (23.5)_final model", "tem_2023_na_wh.dat"))
tem_rep <- readLines(here('2. Base (23.5)_final model', 'sable.rep'))
tem_par <- read_pars(here('2. Base (23.5)_final model', 'tem'))
# Read in disaggregated sex composition data which goes until 2021
compdata_2021 <- readRDS(here('2. Base (23.5)_final model', "SabieRTMB_2021compdata.RDS"))

# Testing with simulated data ---------------------------------------------
# Prepare Data ------------------------------------------------------------
data <- list() # make data list

# Set up dimensions
data$n_regions <- 2 # number of regions
data$ages <- 1:15 # ages
data$lens <- 1 # lengths
data$years <- as.numeric(1:20) # years
data$n_sexes <- 2 # number of sexes
data$n_fish_fleets <- 1 # number of fishery fleets (0 == fixed gear, 1 == trawl gear)
data$n_srv_fleets <- 1 # number of survey fleets (0 == domestic ll survey, 1 == domestic trawl survey, 2 == coop jp ll survey)

# Recruitment stuff
data$init_F_prop <- 0 # initial F proportion for initializing population
data$do_rec_bias_ramp <- 0 # do bias ramp (slot 0 == don't do bias ramp, 1 == do bias ramp)
data$bias_year <- NA
data$sigmaR_switch <- 0
data$sexratio <- as.vector(c(0.5, 0.5)) # recruitment sex ratio (assuming 50,50)

# Movement stuff
data$do_recruits_move = 0 # recruits dont move
data$use_fixed_movement = 0 # use fixed movement
data$Fixed_Movement = sim_out$movement_matrix[,,,,,1]

data$likelihoods <- 1
data$Wt_Catch <- 1 # Catch weights
data$Wt_FishIdx <- 1 # fishery index weights
data$Wt_SrvIdx <- 1 # survey index weights
data$Wt_Rec <- 1 # recruitment weights
data$Wt_F <- 1 # fishing mortality penalty weights

# Biological Processes
# Natural Mortality
data$Use_M_prior <- 0 # use natural mortality prior
data$M_prior <- c(0.1, 0.1) # Mean and CV for M prior

# Weight at age
data$WAA <- aperm(sim_out$WAA[,,,,1], perm = c(2,1,3,4))

# Maturity at age
data$MatAA <- aperm(sim_out$Maturity_AA[,,,,1], perm = c(2,1,3,4))

# Ageing error
data$AgeingError <- diag(1, length(data$ages)) # ageing error matrix

# Size Transition Matrix (Growth)
data$fit_lengths <- 1
data$SizeAgeTrans <- NA

## Observations  -----------------------------------------------------------
### Fishery Observations ----------------------------------------------------

# Catches
data$ObsCatch <- array(aperm(sim_out$Obs_Catch[,,,1], perm = c(2,1)), dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
data$ObsCatch[data$ObsCatch == 0] <- NA # set 0 catches to NA so we aren't fitting
data$UseCatch <- array(1, c(data$n_regions, length(data$years), data$n_fish_fleets))
data$Catch_Constant <- c(0, 0) 

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
data$ObsFishAgeComps <- array(aperm(sim_out$Obs_FishAgeComps[,,,,,1], perm = c(2,1,3,4)), dim = c(data$n_regions, length(data$years), length(data$ages), data$n_sexes, data$n_fish_fleets))
colnames(data$ObsFishAgeComps) <- data$years # define row years
data$UseFishAgeComps <- array(1, dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
colnames(data$UseFishAgeComps) <- data$years # define row years

# Data weighting for fishery age compositions
data$ISS_FishAgeComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
colnames(data$ISS_FishAgeComps) <- data$years # define row years
data$ISS_FishAgeComps[,,1,1] <- 1e2
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
data$FishAgeComps_Type <- array(3, dim = c(data$n_fish_fleets)) # Joint age comps
data$FishLenComps_Type <- array(3, dim = c(data$n_fish_fleets)) # Joint length comps

### Survey Observations -----------------------------------------------------
# Survey Indices
data$ObsSrvIdx <- array(aperm(sim_out$Obs_SrvIdx[,,,1], perm = c(2,1)), dim = c(data$n_regions, length(data$years), data$n_srv_fleets))
data$ObsSrvIdx_SE <- (data$ObsSrvIdx * 0.1) 
data$UseSrvIdx <- array(1, dim = c(data$n_regions, length(data$years), data$n_srv_fleets))
colnames(data$ObsSrvIdx) <- data$years # define row years
colnames(data$ObsSrvIdx_SE) <- data$years # define row years

# Survey Age Comps
data$ObsSrvAgeComps <- array(aperm(sim_out$Obs_SrvAgeComps[,,,,,1], perm = c(2,1,3,4)), dim = c(data$n_regions, length(data$years), length(data$ages), data$n_sexes, data$n_srv_fleets))
colnames(data$ObsSrvAgeComps) <- data$years # define row years
data$UseSrvAgeComps <- array(1, dim = c(data$n_regions, length(data$years), data$n_srv_fleets))

# Data weighting for survey age compositions
data$ISS_SrvAgeComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_srv_fleets))
colnames(data$ISS_SrvAgeComps) <- data$years # define row years
data$ISS_SrvAgeComps[,,1,1] <- 1e2
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
data$SrvLenComps_Type <- array(3, dim = c(data$n_srv_fleets)) # split for both survey fleet length
data$SrvAgeComps_Type <- array(3, dim = c(data$n_srv_fleets)) # split for both survey age length

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

# Prepare Parameters ------------------------------------------------------
parameters <- list()
parameters$dummy <- 1

### Fishery Stuff ---------------------------------------------------------

# Fishing Mortality
parameters$ln_sigmaC <- array(log(0.001), dim = c(data$n_regions, data$n_fish_fleets))
parameters$ln_F_mean <- array(0, dim = c(data$n_regions, data$n_fish_fleets)) # mean fishing mortality
parameters$ln_F_devs <- array(0, dim = c(data$n_regions, length(data$years), data$n_fish_fleets)) # fishing mortality deviations from mean

# Set up continuous fishery selectivity stuff
parameters$ln_fishsel_dev1 <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
parameters$ln_fishsel_dev2 <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
parameters$ln_fishsel_dev1_sd <- array(0.1, dim = c(data$n_regions, data$n_sexes, data$n_fish_fleets))
parameters$ln_fishsel_dev2_sd <- array(0.1, dim = c(data$n_regions, data$n_sexes, data$n_fish_fleets))

# Fixed Gear Fishery three time blocks
max_fish_blks <- 1 # maximum number of fishery blocks for any fleet
max_fish_pars <- 2 # maximum number of fishery fixed parameters for any fleet
parameters$ln_fish_fixed_sel_pars <- array(0, dim = c(data$n_regions, max_fish_pars, max_fish_blks, data$n_sexes, data$n_fish_fleets))
parameters$ln_fish_fixed_sel_pars[,1,,,] = log(5)
parameters$ln_fish_fixed_sel_pars[,2,,,] = log(1)

# Fishery Catchability
# Fixed Gear Fishery Catchability
parameters$ln_fish_q <- array(0, dim = c(data$n_regions, max_fish_blks, data$n_fish_fleets))

### Survey Stuff -----------------------------------------------------
# Survey Selectivity
max_srv_blks <- 1 # maximum number of survey blocks for any fleet
max_srv_pars <- 2 # maximum number of survey fixed parameters for any fleet
parameters$ln_srv_fixed_sel_pars <- array(0, dim = c(data$n_regions, max_srv_pars, max_srv_blks, data$n_sexes, data$n_srv_fleets))
parameters$ln_srv_fixed_sel_pars[,1,,,] = log(4)
parameters$ln_srv_fixed_sel_pars[,2,,,] = log(1)

# Survey Catchability
max_q_srv_blks <- 1 # maximum catchability blocks
# Longline Survey Catchability
parameters$ln_srv_q <- array(log(1), dim = c(data$n_regions,max_q_srv_blks, data$n_srv_fleets))

### Natural Mortality -------------------------------------------------------
parameters$ln_M <- log(0.1) # Female M (base)
parameters$M_offset <- 0 # Male M offset (accidently jittered from OG assessment)

# Recruitment -------------------------------------------------------------
parameters$ln_global_R0 <- log(100) # mean recruitment
parameters$R0_prop <- array(0.5, dim = c(data$n_regions - 1))
parameters$ln_InitDevs <- array(0, dim = c(data$n_regions, length(data$ages) - 2))
parameters$ln_RecDevs <- array(0, dim = c(data$n_regions, length(data$years) - 1))
parameters$ln_sigmaR_early <- log(0.75) # early sigma R
parameters$ln_sigmaR_late <- log(0.75)  # late sigma R

# Comp Likelihood Stuff ---------------------------------------------------
parameters$ln_FishAge_DM_theta <- array(0, dim = c(data$n_regions, data$n_fish_fleets))
parameters$ln_FishLen_DM_theta <- array(0, dim = c(data$n_regions, data$n_fish_fleets))
parameters$ln_SrvAge_DM_theta <- array(0, dim = c(data$n_regions, data$n_srv_fleets))
parameters$ln_SrvLen_DM_theta <- array(0, dim = c(data$n_regions, data$n_srv_fleets))

# Movement Stuff ---------------------------------------------------
parameters$move_pars <- array(0, dim = c(data$n_regions, data$n_regions - 1, length(data$years), length(data$ages), data$n_sexes))
parameters$move_pars[1,,1,1,1] = 0
parameters$move_pars[2,,1,1,1] = 0

# Mapping -----------------------------------------------------------------
mapping <- list()
mapping$dummy <- factor(NA)
mapping$ln_sigmaR_late <- factor(NA)
mapping$ln_sigmaR_early <- factor(NA) # fix early sigma R
mapping$M_offset <- factor(NA) # fix natural mortality offset
mapping$ln_fish_q <- factor(rep(NA, length(parameters$ln_fish_q))) 

# Fixing sigmas for fishery catch and Fdevs here
mapping$ln_sigmaC <- factor(rep(NA, length(parameters$ln_sigmaC)))

# Fixing continuous time-varying selecitvity stuff
mapping$ln_fishsel_dev1 <- factor(rep(NA, length(parameters$ln_fishsel_dev1)))
mapping$ln_fishsel_dev2 <- factor(rep(NA, length(parameters$ln_fishsel_dev2)))
mapping$ln_fishsel_dev1_sd <- factor(rep(NA, length(parameters$ln_fishsel_dev1_sd)))
mapping$ln_fishsel_dev2_sd <- factor(rep(NA, length(parameters$ln_fishsel_dev2_sd)))

# Fixing dirichlet mutlinomial stuff
mapping$ln_FishAge_DM_theta <- factor(rep(NA, length(parameters$ln_FishAge_DM_theta)))
mapping$ln_FishLen_DM_theta <- factor(rep(NA, length(parameters$ln_FishLen_DM_theta)))
mapping$ln_SrvAge_DM_theta <- factor(rep(NA, length(parameters$ln_SrvAge_DM_theta)))
mapping$ln_SrvLen_DM_theta <- factor(rep(NA, length(parameters$ln_SrvLen_DM_theta)))

# Fixing movement stuff
mapping$move_pars = factor(rep(1:2, length.out = prod(dim(parameters$move_pars))))
# mapping$move_pars = factor(rep(NA, length.out = prod(dim(parameters$move_pars))))

# Fixing survey catchability
mapping$ln_srv_q <- factor(rep(1, length(parameters$ln_srv_q)))

# Fixing M
# mapping$ln_M <- factor(NA)

# Fixing fishing mortlaity stuff
# mapping$ln_F_devs = factor(rep(NA, length(parameters$ln_F_devs)))
# mapping$ln_F_mean = factor(rep(NA, length(parameters$ln_F_mean)))
# data$Fmort_dat = array(aperm(sim_out$Fmort[,,,1], perm = c(2,1)), dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
# mapping$ln_fish_fixed_sel_pars = factor(rep(NA, length(parameters$ln_fish_fixed_sel_pars)))

# Fix survey selex
# mapping$ln_srv_fixed_sel_pars = factor(rep(NA, length(parameters$ln_srv_fixed_sel_pars)))

# global density dependence
map_recdevs = parameters$ln_RecDevs
map_recdevs[1,] = 1:length(map_recdevs[1,])
map_recdevs[2,] = map_recdevs[1,]
mapping$ln_RecDevs = factor(map_recdevs)
map_initdevs = parameters$ln_InitDevs
map_initdevs[1,] = 1:length(map_initdevs[1,])
map_initdevs[2,] = map_initdevs[1,]
mapping$ln_InitDevs = factor(map_initdevs)

# mapping$ln_RecDevs = factor(rep(NA, length(parameters$ln_RecDevs)))
# mapping$ln_InitDevs = factor(rep(NA, length(parameters$ln_InitDevs)))

data$srv_q_blocks = data$srv_q_blocks + 1
data$fish_q_blocks = data$fish_q_blocks + 1
data$fish_sel_blocks = data$fish_sel_blocks + 1
data$srv_sel_blocks = data$srv_sel_blocks + 1
data$bias_year = data$bias_year + 1
data$sigmaR_switch = data$sigmaR_switch + 1

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

max(sabie_rtmb_model$sd_rep$gradient.fixed)
sabie_rtmb_model$sd_rep$par.fixed[which.max(sabie_rtmb_model$sd_rep$gradient.fixed)]

sabie_rtmb_model$rep$jnLL 

sabie_rtmb_model$rep$R0
exp(sabie_rtmb_model$sd_rep$par.fixed[names(sabie_rtmb_model$sd_rep$par.fixed) == 'ln_srv_q'])
exp(sabie_rtmb_model$sd_rep$par.fixed[names(sabie_rtmb_model$sd_rep$par.fixed) == 'ln_M'])

par(mfrow = c(2,3))
plot(sabie_rtmb_model$rep$Rec[1,], type = 'l')
lines(rowSums(sim_out$NAA[-21,1,1,,1]), col = 'red')

plot(sabie_rtmb_model$rep$Rec[2,], type = 'l')
lines(rowSums(sim_out$NAA[-21,2,1,,1]), col = 'red')

plot(sabie_rtmb_model$rep$Total_Biom[1,], type = 'l')
lines(sim_out$Total_Biom[,1,1], col = 'red')

plot(sabie_rtmb_model$rep$Total_Biom[2,], type = 'l')
lines(sim_out$Total_Biom[,2,1], col = 'red')

plot(sabie_rtmb_model$rep$SSB[1,], type = 'l')
lines(sim_out$SSB[,1,1], col = 'red')

plot(sabie_rtmb_model$rep$SSB[2,], type = 'l')
lines(sim_out$SSB[,2,1], col = 'red')

sabie_rtmb_model$rep$Movement[,,1,1,1]
movement_matrix[,,1,1,1,1]
#  
# plot(sabie_rtmb_model$rep$srv_sel[1,1,,1,1])
# lines(sim_out$srv_sel[1,1,,1,1,1])
# 
# plot(sabie_rtmb_model$rep$fish_sel[1,1,,1,1])
# lines(sim_out$fish_sel[1,1,,1,1,1])


# Get population dynamics -------------------------------------------------
# ages <- 1:30 # ages
# years <- 1960:2021 # years
# n_regions <- 1 # number of regions
# waa_f <- tem_dat$growthmat$wt.f.block1 # female WAA
# waa_m <- tem_dat$growthmat$wt.m.block1 # male WAA
# mat <- tem_dat$growthmat$mage.block2 # maturity
# m_f <- exp(tem_par$coefficients[names(tem_par$coefficients) == "logm"]) # log natural mortality
# m_m <- m_f-0.00819813327862 # male M (adding offset from accidental jitter)
# 
# ### Recruitment -------------------------------------------------------------
# rec <- tem_dat$t.series$Recr # total recruitment
# init_age <- matrix(cbind(tem_dat$natage.female[1,-1] , tem_dat$natage.male[1,-1] ), ncol = 2) # initial ages
# mean_rec <- tem_par$coefficients[names(tem_par$coefficients) == "log_mean_rec"] # mean recruitment
# init_rec_devs <- rev(tem_par$coefficients[str_detect(names(tem_par$coefficients), "rec_dev") ][1:28]) # initial devs
# rec_devs <- tem_par$coefficients[str_detect(names(tem_par$coefficients), "rec_dev")][c(29:91)] # rec devs
# 
# ### Fishing Mortality -------------------------------------------------------
# total_f <- tem_dat$t.series$fmort # total F
# mean_ll_fish <- tem_par$coefficients[names(tem_par$coefficients) == "log_avg_F_fish1"] # Average F LL
# devs_ll_fish <- tem_par$coefficients[str_detect(names(tem_par$coefficients), "log_F_devs_fish1")] # Devs LL
# mean_ll_trwl <- tem_par$coefficients[names(tem_par$coefficients) == "log_avg_F_fish3"] # Average F Trawl
# devs_ll_trwl <- tem_par$coefficients[str_detect(names(tem_par$coefficients), "log_F_devs_fish3")] # Devs Trawl
# 
# # Prepare Data ------------------------------------------------------------
# data <- list() # make data list
# 
# # Set up dimensions
# data$n_regions <- n_regions # number of regions
# data$ages <- ages # ages
# data$lens <- seq(41, 99, 2) # lengths
# data$years <- as.numeric(years) # years
# data$n_sexes <- 2 # number of sexes
# data$n_fish_fleets <- 2 # number of fishery fleets (0 == fixed gear, 1 == trawl gear)
# data$n_srv_fleets <- 3 # number of survey fleets (0 == domestic ll survey, 1 == domestic trawl survey, 2 == coop jp ll survey)
# data$fit_lengths <- 0
# 
# # Recruitment stuff
# data$init_iter <- 100 # Number of iterations to get equilibrium age structure w/ movement
# data$init_F_prop <- 0.1 # initial F proportion for initializing population
# data$do_rec_bias_ramp <- 1 # do bias ramp (slot 0 == don't do bias ramp, 1 == do bias ramp)
# # breakpoints for bias ramp (0 == no bias ramp - 1960 - 1980, 1 == ascending limb of bias ramp - 1980 - 1990,
# # 2 == full bias correction - 1990 - 2022, == 3 no bias correction - terminal year of recruitment estimate)
# data$bias_year <- c(length(1960:1979), length(1960:1989), (length(1960:2022) - 5), length(1960:2021) - 2)
# data$sigmaR_switch <- as.integer(length(1960:1975)) - 1 # when to switch sigmaR from 0.4 to a larger value
# data$sexratio <- as.vector(c(0.5, 0.5)) # recruitment sex ratio (assuming 50,50)
# 
# # Movement stuff
# data$do_recruits_move = 0 # recruits dont move
# 
# data$likelihoods <- 0
# # Weights for likelihoods
# if(data$likelihoods == 0) {
#   data$Wt_Catch <- 50 # Catch weights
#   data$Wt_FishIdx <- 0.448 # fishery index weights
#   data$Wt_SrvIdx <- 0.448 # survey index weights
#   data$Wt_Rec <- 1.5 # recruitment weights
#   data$Wt_F <- 0.1 # fishing mortality penalty weights
# }
# 
# if(data$likelihoods == 1) {
#   data$Wt_Catch <- 1 # Catch weights
#   data$Wt_FishIdx <- 1 # fishery index weights
#   data$Wt_SrvIdx <- 1 # survey index weights
#   data$Wt_Rec <- 1 # recruitment weights
#   data$Wt_F <- 1 # fishing mortality penalty weights
# }
# 
# # Biological Processes
# # Natural Mortality
# data$Use_M_prior <- 1 # use natural mortality prior
# data$M_prior <- c(0.1, 0.1) # Mean and CV for M prior
# 
# # Weight at age
# data$WAA <- array(NA, dim = c(data$n_regions, length(data$years), length(data$ages), data$n_sexes)) # weight-at-age array
# data$WAA[1,,,1] <- matrix(waa_f, nrow = length(1:62), ncol = length(1:30), byrow = TRUE) # female weight-at-age
# data$WAA[1,,,2] <- matrix(waa_m, nrow = length(1:62), ncol = length(1:30), byrow = TRUE) # male weight-at-age
# 
# # Maturity at age
# data$MatAA <- array(0,  dim = c(data$n_regions, length(data$years), length(data$ages), data$n_sexes)) # maturity at age
# data$MatAA[1,1:62,,1] <- aperm(replicate(length(1:62), mat), c(2,1))
# data$MatAA[1,1:62,,2] <- aperm(replicate(length(1:62), mat), c(2,1))
# 
# # Ageing error
# data$AgeingError <- as.matrix(ageing_dat$age_error) # ageing error matrix
# 
# # Size Transition Matrix (Growth)
# data$SizeAgeTrans <- array(0, dim = c(data$n_regions, length(data$years), length(data$lens), length(data$ages), data$n_sexes)) # size age transition matrix
# 
# # Size age transition, first time block
# data$SizeAgeTrans[1,1:35,,,1] <- aperm(replicate(length(1:35), tem_dat$sizeage.f.block1), perm = c(3,2,1)) # female
# data$SizeAgeTrans[1,1:35,,,2] <- aperm(replicate(length(1:35), tem_dat$sizeage.m.block1), perm = c(3,2,1)) # male
# 
# # Size age transition, second time block
# data$SizeAgeTrans[1,36:62,,,1] <- aperm(replicate(length(36:62), tem_dat$sizeage.f.block2), perm = c(3,2,1)) # female
# data$SizeAgeTrans[1,36:62,,,2] <- aperm(replicate(length(36:62), tem_dat$sizeage.m.block2), perm = c(3,2,1)) # male
# 
# # movement stuff
# data$use_fixed_movement = 1
# data$Fixed_Movement = array(1, dim = c(data$n_regions, data$n_regions, length(data$years), length(data$ages), data$n_sexes))
# 
# ## Observations  -----------------------------------------------------------
# ### Fishery Observations ----------------------------------------------------
# 
# # Catches
# data$ObsCatch <- array(NA, c(data$n_regions, length(1:62), data$n_fish_fleets))
# data$ObsCatch[1,,1] <- as.numeric(strsplit(tem_admb_dat[46], split = " ")[[1]])[1:62] # fixed gear catches
# data$ObsCatch[1,,2] <- as.numeric(strsplit(tem_admb_dat[48], split = " ")[[1]])[1:62] # trawl gear catches
# data$ObsCatch[data$ObsCatch == 0] <- NA # set 0 catches to NA so we aren't fitting
# data$UseCatch <- array(1, c(data$n_regions, length(1:62), data$n_fish_fleets))
# data$UseCatch[1,1:3,2] <- 0 # don't use first 3 years of trawl catches (0)
# data$Catch_Constant <- c(0.01, 0.8) # catch constants to add to likeilhood (== 0 fixed gear adds 0.01 and == 1, trawl gear adds 0.8)
# 
# # Fishery Indices
# data$ObsFishIdx <- array(NA, c(data$n_regions, length(1:62), data$n_fish_fleets))
# data$ObsFishIdx_SE <- array(NA, c(data$n_regions, length(1:62), data$n_fish_fleets))
# data$UseFishIdx <- array(1, c(data$n_regions, length(1:62), data$n_fish_fleets))
# colnames(data$ObsFishIdx) <- data$years # define row years
# colnames(data$ObsFishIdx_SE) <- data$years # define row years
# 
# # Fixed gear fishery CPUE
# # Domestic LL + Japanese Fishery (pre 1995)
# data$share_sel <- 0 # share fishery index selectivity for the first time block and fleet across sexes (specific to sablefish)
# # Basically, the pre 1995 selectivity only utilizes the female selectivity from the first time block - this should be changed to
# # use both sexes however... (oh well!)
# data$ObsFishIdx[1,colnames(data$ObsFishIdx) %in% rownames(tem_dat$obssrv6),1] <- as.numeric(strsplit(tem_admb_dat[144], split = " ")[[1]])
# data$ObsFishIdx_SE[1,colnames(data$ObsFishIdx_SE) %in% rownames(tem_dat$obssrv6),1] <- as.numeric(strsplit(tem_admb_dat[146], split = " ")[[1]])
# 
# # Domestic LL Fishery (after 1995)
# data$ObsFishIdx[1,colnames(data$ObsFishIdx) %in% rownames(tem_dat$obssrv5)[1:27],1] <- as.numeric(strsplit(tem_admb_dat[127], split = " ")[[1]])[1:27]
# data$ObsFishIdx_SE[1,colnames(data$ObsFishIdx_SE) %in% rownames(tem_dat$obssrv5)[1:27],1] <- as.numeric(strsplit(tem_admb_dat[129], split = " ")[[1]])[1:27]
# data$UseFishIdx[is.na(data$ObsFishIdx)] <- 0 # don't fit if missing data
# 
# # Fishery Age Comps
# # Note that NA is in trawl fishery slot so it doesn't fit
# data$ObsFishAgeComps <- array(NA, dim = c(data$n_regions, length(data$years), length(data$ages), data$n_sexes, data$n_fish_fleets))
# colnames(data$ObsFishAgeComps) <- data$years # define row years
# 
# # Fishery Age Comps (Sex Specific)
# data$ObsFishAgeComps[1,colnames(data$ObsFishAgeComps) %in% rownames(compdata_2021$FishAge),,,] <- compdata_2021$FishAge # Observed fixed gear fishery age comps (sex specific)
# data$UseFishAgeComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
# colnames(data$UseFishAgeComps) <- data$years # define row years
# data$UseFishAgeComps[1,which(rowSums(compdata_2021$FishAge) != 0),1] <- 1 # only fit if have age comp data
# 
# # Data weighting for fishery age compositions
# data$ISS_FishAgeComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
# colnames(data$ISS_FishAgeComps) <- data$years # define row years
# data$ISS_FishAgeComps[1,which(rowSums(compdata_2021$FishAge) != 0),1,1] <- 20 # Assuming constant ISS of 20 for fixed gear fishery ages (sex-specific - still indexing 1 because fitting single vector right now)
# data$Wt_FishAgeComps <- array(NA, dim = c(data$n_regions, data$n_sexes, data$n_fish_fleets)) # weights for fishery age comps
# data$Wt_FishAgeComps[1,1,1] <- 1.195228 # Weight for fixed gear age comps
# 
# # Fishery Length Comps
# data$ObsFishLenComps <- array(NA, dim = c(data$n_regions, length(data$years), length(data$lens), data$n_sexes, data$n_fish_fleets))
# colnames(data$ObsFishLenComps) <- data$years # define row years
# 
# # observed fixed and trawl gear fishery length comps (joint)
# data$ObsFishLenComps[1,colnames(data$ObsFishLenComps) %in% rownames(compdata_2021$FishLen),,,] <- compdata_2021$FishLen
# 
# # Use fishery length comp controls
# data$UseFishLenComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
# colnames(data$UseFishLenComps) <- data$years # define row years
# data$UseFishLenComps[1,which(rowSums(compdata_2021$FishLen[,,,1]) != 0),1] <- 1 # only fit if have length comp data (fixed gear)
# data$UseFishLenComps[1,which(rowSums(compdata_2021$FishLen[,,,2]) != 0),2] <- 1 # only fit if have length comp data (trawl gear)
# 
# # Data weighting for fishery length compositions
# data$ISS_FishLenComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
# colnames(data$ISS_FishLenComps) <- data$years # define row years
# 
# # joint iss
# data$ISS_FishLenComps[1,which(rowSums(compdata_2021$FishLen[,,,1]) != 0),1,1] <- 20 # Assuming constant ISS of 20 for fixed gear fishery
# data$ISS_FishLenComps[1,which(rowSums(compdata_2021$FishLen[,,,2]) != 0),1,2] <- 20 # Assuming constant ISS of 20 for trawl gear fishery
# data$Wt_FishLenComps <- array(0, dim = c(data$n_regions, data$n_sexes, data$n_fish_fleets)) # weights for fishery age comps
# data$Wt_FishLenComps[1,1,1] <- 2.099563 # Weight for fixed gear len comps females
# data$Wt_FishLenComps[1,1,2] <- 0.6589962 # Weight for trawl gear len comps females
# 
# # Composition munging stuff
# data$FishAgeComps_LikeType <- array(0, dim = c(data$n_fish_fleets)) # multinomial for both fleets
# data$FishLenComps_LikeType <- array(0, dim = c(data$n_fish_fleets)) # multinomial for both fleets
# data$FishAgeComps_Type <- array(2, dim = c(data$n_fish_fleets)) # Joint age comps
# data$FishLenComps_Type <- array(2, dim = c(data$n_fish_fleets)) # Joint length comps
# 
# ### Survey Observations -----------------------------------------------------
# # Survey Indices
# data$ObsSrvIdx <- array(NA, dim = c(data$n_regions, length(1:62), data$n_srv_fleets))
# data$ObsSrvIdx_SE <- array(NA, dim = c(data$n_regions, length(1:62), data$n_srv_fleets))
# data$UseSrvIdx <- array(1, dim = c(data$n_regions, length(1:62), data$n_srv_fleets))
# colnames(data$ObsSrvIdx) <- data$years # define row years
# colnames(data$ObsSrvIdx_SE) <- data$years # define row years
# 
# # Domestic LL Survey post 1995
# data$ObsSrvIdx[1,colnames(data$ObsSrvIdx) %in% rownames(tem_dat$obssrv3)[1:32],1] <- as.numeric(strsplit(tem_admb_dat[93], split = " ")[[1]])[1:32]
# data$ObsSrvIdx_SE[1,colnames(data$ObsSrvIdx_SE) %in% rownames(tem_dat$obssrv3)[1:32],1] <- as.numeric(strsplit(tem_admb_dat[95], split = " ")[[1]])[1:32]
# 
# # Domestic Trawl Survey
# data$ObsSrvIdx[1,colnames(data$ObsSrvIdx) %in% rownames(tem_dat$obssrv7)[1:14],2] <- as.numeric(strsplit(tem_admb_dat[161], split = " ")[[1]])[1:14]
# data$ObsSrvIdx_SE[1,colnames(data$ObsFishIdx_SE) %in% rownames(tem_dat$obssrv7)[1:14],2] <- as.numeric(strsplit(tem_admb_dat[163], split = " ")[[1]])[1:14]
# 
# # Coop LL Survey pre 1995
# data$ObsSrvIdx[1,colnames(data$ObsSrvIdx) %in% rownames(tem_dat$obssrv4),3] <- as.numeric(strsplit(tem_admb_dat[110], split = " ")[[1]])
# data$ObsSrvIdx_SE[1,colnames(data$ObsSrvIdx_SE) %in% rownames(tem_dat$obssrv4),3] <- as.numeric(strsplit(tem_admb_dat[112], split = " ")[[1]])
# data$UseSrvIdx[is.na(data$ObsSrvIdx)] <- 0 # don't fit if missing data
# 
# # Survey Age Comps
# data$ObsSrvAgeComps <- array(NA, dim = c(data$n_regions, length(data$years), length(data$ages), data$n_sexes, data$n_srv_fleets))
# colnames(data$ObsSrvAgeComps) <- data$years # define row years
# 
# # joint ages
# data$ObsSrvAgeComps[1,colnames(data$ObsSrvAgeComps) %in% rownames(compdata_2021$SrvAge),,,] <- compdata_2021$SrvAge
# 
# data$UseSrvAgeComps <- array(0, dim = c(data$n_regions, length(1:62), data$n_srv_fleets))
# colnames(data$UseSrvAgeComps) <- data$years # define row years
# data$UseSrvAgeComps[1,which(rowSums(compdata_2021$SrvAge[,,,1]) != 0),1] <- 1
# data$UseSrvAgeComps[1,which(rowSums(compdata_2021$SrvAge[,,,3]) != 0),3] <- 1
# 
# # Data weighting for fishery age compositions
# data$ISS_SrvAgeComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_srv_fleets))
# colnames(data$ISS_SrvAgeComps) <- data$years # define row years
# 
# data$ISS_SrvAgeComps[1,which(rowSums(compdata_2021$SrvAge[,,,1]) != 0),1,1] <- 20 # Assuming constant ISS of 20 for domestic ll survey ages (aggregated)
# data$ISS_SrvAgeComps[1,which(rowSums(compdata_2021$SrvAge[,,,3]) != 0),1,3] <- 20 # Assuming constant ISS of 20 for coop jp survey ages (aggregated)
# 
# data$Wt_SrvAgeComps <- array(NA, dim = c(data$n_regions, data$n_sexes, data$n_srv_fleets)) # weights for survey age comps
# data$Wt_SrvAgeComps[1,1,1] <- 5.399427 # Weight for domestic survey ll gear age comps
# data$Wt_SrvAgeComps[1,1,3] <- 3.240157 # Weight for coop jp survey ll gear age comps
# 
# # Survey Length Comps
# data$ObsSrvLenComps <- array(NA, dim = c(data$n_regions,length(data$years), length(data$lens), data$n_sexes, data$n_srv_fleets))
# colnames(data$ObsSrvLenComps) <- data$years # define row years
# 
# # observed domestic survey ll length comps
# data$ObsSrvLenComps[1,colnames(data$ObsSrvLenComps) %in% rownames(tem_dat$olc.srv1.f)[1:32],,1,1] <- tem_dat$olc.srv1.f[1:32,] # females
# data$ObsSrvLenComps[1,colnames(data$ObsSrvLenComps) %in% rownames(tem_dat$olc.srv1.m)[1:32],,2,1] <- tem_dat$olc.srv1.m[1:32,] # males
# 
# # observed domestic trawl survey length comps
# data$ObsSrvLenComps[1,colnames(data$ObsSrvLenComps) %in% rownames(tem_dat$olc.srv7.f)[1:14],,1,2] <- tem_dat$olc.srv7.f[1:14,] # females
# data$ObsSrvLenComps[1,colnames(data$ObsSrvLenComps) %in% rownames(tem_dat$olc.srv7.m)[1:14],,2,2] <- tem_dat$olc.srv7.m[1:14,] # males
# 
# # observed coop jp ll survey length comps
# srv_trawl_iss <- as.numeric(strsplit(tem_admb_dat[580], split = " ")[[1]]) # get ISS from trawl survey
# data$ObsSrvLenComps[1,rownames(data$ObsSrvLenComps) %in% rownames(tem_dat$olc.srv2.f),,1,3] <- tem_dat$olc.srv2.f # females
# data$ObsSrvLenComps[1,rownames(data$ObsSrvLenComps) %in% rownames(tem_dat$olc.srv2.m),,2,3] <- tem_dat$olc.srv2.m # males
# 
# # Use survey length comp controls
# data$UseSrvLenComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_srv_fleets))
# colnames(data$UseSrvLenComps) <- data$years # define row years
# data$UseSrvLenComps[1,colnames(data$UseSrvLenComps) %in% rownames(tem_dat$olc.srv1.f)[1:32],1] <- 0 # only fit if have len comp data
# data$UseSrvLenComps[1,colnames(data$UseSrvLenComps) %in% rownames(tem_dat$olc.srv7.f)[1:14],2] <- 1 # only fit if have len comp data
# data$UseSrvLenComps[1,colnames(data$UseSrvLenComps) %in% rownames(tem_dat$olc.srv2.f),3] <- 0 # only fit if have len comp data
# 
# # Data weighting for survey length compositions
# data$ISS_SrvLenComps <- array(0, dim = c(data$n_regions,length(data$years), data$n_sexes, data$n_srv_fleets))
# colnames(data$ISS_SrvLenComps) <- data$years # define row years
# data$ISS_SrvLenComps[1,colnames(data$ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv1.f)[1:32],1,1] <- 20 # Assuming constant ISS of 20 for domestic ll survey females
# data$ISS_SrvLenComps[1,colnames(data$ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv1.m)[1:32],2,1] <- 20 # Assuming constant ISS of 20 for domestic ll survey males
# data$ISS_SrvLenComps[1,colnames(data$ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv7.f)[1:14],1,2] <- srv_trawl_iss[1:14] # Assuming constant ISS of 20 for domestic trawl survey females
# data$ISS_SrvLenComps[1,colnames(data$ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv7.m)[1:14],2,2] <- srv_trawl_iss[1:14] # Assuming constant ISS of 20 for domestic trawl survey males
# data$ISS_SrvLenComps[1,colnames(data$ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv2.f),1,3] <- 20 # Assuming constant ISS of 20 for coop jp ll survey females
# data$ISS_SrvLenComps[1,colnames(data$ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv2.m),2,3] <- 20 # Assuming constant ISS of 20 for coop jp ll survey males
# 
# # data weights
# data$Wt_SrvLenComps <- array(0, dim = c(data$n_regions, data$n_sexes, data$n_srv_fleets)) # weights for survey age comps
# data$Wt_SrvLenComps[1,1,1] <- 1 # Weight for domestic ll survey len comps females
# data$Wt_SrvLenComps[1,2,1] <- 1 # Weight for domestic ll survey len comps males
# data$Wt_SrvLenComps[1,1,2] <- 0.5916596  # Weight for domestic trawl survey len comps females
# data$Wt_SrvLenComps[1,2,2] <- 0.3986298  # Weight for domestic trawl survey len comps males
# data$Wt_SrvLenComps[1,1,3] <- 1 # Weight for coop jp ll survey len comps females
# data$Wt_SrvLenComps[1,2,3] <- 1  # Weight for coop jp ll survey len comps males
# 
# # Composition munging stuff
# data$SrvAgeComps_LikeType <- array(0, dim = c(data$n_srv_fleets)) # multinomial for both survey fleet
# data$SrvLenComps_LikeType <- array(0, dim = c(data$n_srv_fleets)) # multinomial for both survey fleet
# data$SrvLenComps_Type <- array(1, dim = c(data$n_srv_fleets)) # split for both survey fleet length
# data$SrvAgeComps_Type <- array(2, dim = c(data$n_srv_fleets)) # split for both survey age length
# 
# ### Fishery Stuff -----------------------------------------------------
# # Selectivity
# data$cont_tv_fish_sel <- array(0, dim = c(data$n_regions, data$n_fish_fleets)) # no timevarying selex continously
# # Time Block Specification
# data$fish_sel_blocks <- array(NA, dim = c(data$n_regions, length(1:62), data$n_fish_fleets))
# data$fish_sel_blocks[1,1:35,1] <- 0 # block one fishery ll selex
# data$fish_sel_blocks[1,36:56,1] <- 1 # block two fishery ll selex
# data$fish_sel_blocks[1,57:62,1] <- 2 # block three fishery ll selex
# data$fish_sel_blocks[1,,2] <- 0 # block two fishery ll selex
# 
# # Selectivity Model
# data$fish_sel_model <- array(NA, dim = c(data$n_regions, length(1:62), data$n_fish_fleets))
# data$fish_sel_model[1,,1] <- 0 # Logistic selectivity for fixed gear
# data$fish_sel_model[1,,2] <- 1 # Gamma selectivity for trawl fishery
# 
# # Catchability
# # Time Block Specification
# data$fish_q_blocks <- data$fish_sel_blocks # catchability blocks are same as selectivity blocks
# data$fish_idx_type <- array(1, dim = c(data$n_regions, data$n_fish_fleets)) # fishery index type (0 == abundance, 1 == biomass)
# 
# ### Survey Selectivity ------------------------------------------------------
# # Time Block Specification
# data$srv_sel_blocks <- array(NA, dim = c(data$n_regions, length(1:62), data$n_srv_fleets))
# data$srv_sel_blocks[1,1:56,1] <- 0 # block one survey ll selex
# data$srv_sel_blocks[1,57:62,1] <- 1 # block two survey ll selex
# data$srv_sel_blocks[1,,2] <- 0 # block one trawl survey selex
# data$srv_sel_blocks[1,,3] <- 0 # block one coop jp ll survey selex
# 
# # Selectivity Model
# data$srv_sel_model <- array(NA, dim = c(data$n_regions, length(1:62), data$n_srv_fleets))
# data$srv_sel_model[1,,1] <- 0 # Logistic selectivity for longline survey
# data$srv_sel_model[1,,2] <- 2 # power function selectivity for trawl survey
# data$srv_sel_model[1,,3] <- 0 # logistic selectivity for jp ll survey
# 
# ### Survey Catchability ----------------------------------------------------
# # Time Block Specification
# data$srv_q_blocks <- array(0, dim = c(data$n_regions, length(1:62), data$n_srv_fleets)) # catchability blocks are same as selectivity blocks
# data$srv_idx_type <- array(NA, dim = c(data$n_regions, data$n_srv_fleets))
# data$srv_idx_type[] <- c(0,1,0) # survey index type (0 == abundance, 1 == biomass)
# 
# # Prepare Parameters ------------------------------------------------------
# parameters <- list()
# parameters$dummy <- 1
# 
# ### Fishery Stuff ---------------------------------------------------------
# 
# # Fishing Mortality
# parameters$ln_sigmaC <- array(0.1, dim = c(data$n_regions, data$n_fish_fleets))
# parameters$ln_F_mean <- array(c(mean_ll_fish, mean_ll_trwl), dim = c(data$n_regions, data$n_fish_fleets)) # mean fishing mortality
# parameters$ln_F_devs <- array(0, dim = c(data$n_regions, length(data$years), data$n_fish_fleets)) # fishing mortality deviations from mean
# parameters$ln_F_devs[1,,1] <- devs_ll_fish[1:62] # longline fishery Fs
# parameters$ln_F_devs[1,,2] <- c(rep(0, 3), devs_ll_trwl[1:59] ) # trwl fishery Fs
# 
# # Set up continuous fishery selectivity stuff
# parameters$ln_fishsel_dev1 <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
# parameters$ln_fishsel_dev2 <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
# parameters$ln_fishsel_dev1_sd <- array(0.1, dim = c(data$n_regions, data$n_sexes, data$n_fish_fleets))
# parameters$ln_fishsel_dev2_sd <- array(0.1, dim = c(data$n_regions, data$n_sexes, data$n_fish_fleets))
# 
# # Fixed Gear Fishery three time blocks
# max_fish_blks <- 3 # maximum number of fishery blocks for any fleet
# max_fish_pars <- 2 # maximum number of fishery fixed parameters for any fleet
# 
# parameters$ln_fish_fixed_sel_pars <- array(0, dim = c(data$n_regions, max_fish_pars, max_fish_blks, data$n_sexes, data$n_fish_fleets))
# parameters$ln_fish_fixed_sel_pars[1,,1,1,1] <- c(1.4226e+00, -6.8808e-01) # Female longline fishery first time block parameters (a50 and then delta)
# parameters$ln_fish_fixed_sel_pars[1,,2,1,1] <- c(1.2219e+00, 5.3916e-01) # Female longline fishery second time block parameters (a50 and then delta)
# parameters$ln_fish_fixed_sel_pars[1,,3,1,1] <- c(1.2219e+00, 5.3916e-01) # Female longline fishery third time block parameters (a50 and then delta)
# 
# parameters$ln_fish_fixed_sel_pars[1,,1,2,1] <- c(1.9339e+00, -6.8808e-01) # Male longline fishery first time block parameters (a50 and then delta)
# parameters$ln_fish_fixed_sel_pars[1,,2,2,1] <- c(1.4907e+00, -9.5989e-02) # Male longline fishery second time block parameters (a50 and then delta)
# parameters$ln_fish_fixed_sel_pars[1,,3,2,1] <- c(1.0702e+00, -3.8718e-01) # Male longline fishery third time block parameters (a50 and then delta)
# 
# # Trawl Fishery kept the same across blocks
# parameters$ln_fish_fixed_sel_pars[1,,1,1,2] <- c(1.7735e+00, 3) # Female trawl fishery first time block parameters (amax and then power)
# parameters$ln_fish_fixed_sel_pars[1,,2,1,2] <- c(1.7735e+00, 3) # Female trawl fishery first time block parameters (amax and then power)
# parameters$ln_fish_fixed_sel_pars[1,,3,1,2] <- c(1.7735e+00, 3) # Female trawl fishery first time block parameters (amax and then power)
# 
# # Trawl Fishery kept the same across blocks
# parameters$ln_fish_fixed_sel_pars[1,,1,2,2] <- c(2.1060e+00, 3) # Male trawl fishery first time block parameters (amax and then power)
# parameters$ln_fish_fixed_sel_pars[1,,2,2,2] <- c(2.1060e+00, 3) # Male trawl fishery first time block parameters (amax and then power)
# parameters$ln_fish_fixed_sel_pars[1,,3,2,2] <- c(2.1060e+00, 3) # Male trawl fishery first time block parameters (amax and then power)
# 
# # Fishery Catchability
# # Fixed Gear Fishery Catchability
# parameters$ln_fish_q <- array(0, dim = c(data$n_regions, max_fish_blks, data$n_fish_fleets))
# parameters$ln_fish_q[1,1,1] <- tem_par$coefficients[names(tem_par$coefficients) == "log_q_srv6"] # first fixed gear fishery time block
# parameters$ln_fish_q[1,2,1] <- tem_par$coefficients[names(tem_par$coefficients) == "log_q_srv8"] # second fixed gear fishery time block
# parameters$ln_fish_q[1,3,1] <- tem_par$coefficients[227] # third fixed gear fishery time block
# 
# # Trawl Gear Fishery Catchability
# parameters$ln_fish_q[1,,2] <- 0 # trawl gear fishery time block (none and no index used)
# 
# ### Survey Stuff -----------------------------------------------------
# # Survey Selectivity
# max_srv_blks <- 2 # maximum number of survey blocks for any fleet
# max_srv_pars <- 2 # maximum number of survey fixed parameters for any fleet
# parameters$ln_srv_fixed_sel_pars <- array(0, dim = c(data$n_regions, max_srv_pars, max_srv_blks, data$n_sexes, data$n_srv_fleets))
# 
# # Longline Survey - 2 time blocks
# parameters$ln_srv_fixed_sel_pars[1,,1,1,1] <- c(tem_par$coefficients[names(tem_par$coefficients) == "log_a50_srv1_f"], tem_par$coefficients[names(tem_par$coefficients) == "log_delta_srv1_f"]) # Female longline survey first time block parameters (a50 and then delta)
# parameters$ln_srv_fixed_sel_pars[1,,2,1,1] <- c(tem_par$coefficients[names(tem_par$coefficients) == "log_a50_srv10_f"], tem_par$coefficients[names(tem_par$coefficients) == "log_delta_srv1_f"]) # Female longline survey second time block parameters (a50 and then delta)
# parameters$ln_srv_fixed_sel_pars[1,,1,2,1] <- c(tem_par$coefficients[names(tem_par$coefficients) == "log_a50_srv1_m"], tem_par$coefficients[names(tem_par$coefficients) == "log_delta_srv1_m"]) # Male longline survey first time block parameters (a50 and then delta)
# parameters$ln_srv_fixed_sel_pars[1,,2,2,1] <- c(tem_par$coefficients[names(tem_par$coefficients) == "log_a50_srv10_m"], tem_par$coefficients[names(tem_par$coefficients) == "log_delta_srv1_m"]) # Male longline survey second time block parameters (a50 and then delta)
# 
# # Trawl Survey - power function, single time block
# parameters$ln_srv_fixed_sel_pars[1,,,1,2] <- tem_par$coefficients[names(tem_par$coefficients) == "log_a50_srv7_f"] # Female trawl survey first time block parameters (power)
# parameters$ln_srv_fixed_sel_pars[1,,,2,2] <- tem_par$coefficients[names(tem_par$coefficients) == "log_a50_srv7_m"] # Male trawl survey first time block parameters (power)
# 
# # Coop JP Survey (Logistic) Single time block
# parameters$ln_srv_fixed_sel_pars[1,,,1,3] <- c(0.953479618491, tem_par$coefficients[names(tem_par$coefficients) == "log_delta_srv1_f"]) # Female trawl survey first time block parameters (a50, delta)
# parameters$ln_srv_fixed_sel_pars[1,,,2,3] <- c(1.23386822176, tem_par$coefficients[names(tem_par$coefficients) == "log_delta_srv1_m"]) # Male trawl survey first time block parameters (a50, delta)
# 
# # Survey Catchability
# max_q_srv_blks <- 1 # maximum catchability blocks
# 
# # Longline Survey Catchability
# parameters$ln_srv_q <- array(0, dim = c(data$n_regions,max_q_srv_blks, data$n_srv_fleets))
# parameters$ln_srv_q[1,,1] <- tem_par$coefficients[names(tem_par$coefficients) == "log_q_srv1"] # first longline survey (domestic) time block
# 
# # Trawl Survey Catchability
# parameters$ln_srv_q[1,,2] <- tem_par$coefficients[names(tem_par$coefficients) == "log_q_srv7"] # trawl survey time block
# 
# # Coop JP Survey Catchability
# parameters$ln_srv_q[1,,3] <- tem_par$coefficients[names(tem_par$coefficients) == "log_q_srv2"]
# 
# ### Natural Mortality -------------------------------------------------------
# parameters$ln_M <- tem_par$coefficients[names(tem_par$coefficients) == "logm"] # Female M (base)
# parameters$M_offset <- -0.00819813327862 # Male M offset (accidently jittered from OG assessment)
# 
# # Recruitment -------------------------------------------------------------
# parameters$ln_global_R0 <- log(100) # mean recruitment
# parameters$R0_prop <- array(1, dim = c(data$n_regions - 1))
# parameters$ln_InitDevs <- array(rev(tem_par$coefficients[str_detect(names(tem_par$coefficients), "rec_dev")][c(1:28)]),
#                                 dim = c(data$n_regions, length(data$ages) - 2))
# parameters$ln_RecDevs <- array(tem_par$coefficients[str_detect(names(tem_par$coefficients), "rec_dev")][c(29:89)],
#                                dim = c(data$n_regions, length(data$years) - 1))
# parameters$ln_sigmaR_early <- log(0.4) # early sigma R
# parameters$ln_sigmaR_late <- tem_par$coefficients[names(tem_par$coefficients) == "log_sigr"]  # late sigma R
# 
# # Comp Likelihood Stuff ---------------------------------------------------
# parameters$ln_FishAge_DM_theta <- array(0, dim = c(data$n_regions, data$n_fish_fleets))
# parameters$ln_FishLen_DM_theta <- array(0, dim = c(data$n_regions, data$n_fish_fleets))
# parameters$ln_SrvAge_DM_theta <- array(0, dim = c(data$n_regions, data$n_srv_fleets))
# parameters$ln_SrvLen_DM_theta <- array(0, dim = c(data$n_regions, data$n_srv_fleets))
# 
# # Movement Stuff ---------------------------------------------------
# parameters$move_pars <- array(-999, dim = c(data$n_regions, data$n_regions - 1, length(data$years), length(data$ages), data$n_sexes))
# 
# # RTMB Port ---------------------------------------------------------------
# mapping <- list()
# mapping$ln_sigmaR_late <- factor(NA)
# mapping$ln_sigmaR_early <- factor(NA) # fix early sigma R
# mapping$R0_prop = factor(rep(NA, length(parameters$R0_prop)))
# mapping$M_offset <- factor(NA) # fix natural mortality offset
# mapping$dummy <- factor(NA)
# mapping$ln_F_devs <- factor(c(1:length(data$years), rep(NA, 3), 65:(126-3))) # fix all fishery mort parameters for now
# mapping$ln_fish_q <- factor(c(1,2,3,rep(NA,3))) # estimate catchabilities only for first fishery with index
# # sharing delta across sexes from early domestic fishery (first time block)
# # also fixing parameters so that no time block for trawl fishery
# mapping$ln_fish_fixed_sel_pars <- factor(c(1:7, 2, 8:10, 6, rep(12:13,3), rep(c(14,13),3)))
# # mapping$ln_fish_fixed_sel_pars <- factor(c(1:12, rep(13:14,3), rep(c(15,16),3)))
# 
# # ll survey, share delta female (index 2) across time blocks and to the coop jp ll survey delta
# # ll survey, share delta male (index 5) across time blocks and to the coop jp ll survey delta
# # coop jp survey does not estimate parameters and shares deltas with longline survey
# # single time block with trawl survey and only one parameter hence, only one parameter estimated across blocks (indices 7 and 8)
# mapping$ln_srv_fixed_sel_pars <- factor(c(1:3, 2, 4:6, 5,
#                                           rep(7,4), rep(8, 4),
#                                           rep(c(NA,2), 2), rep(c(NA, 5), 2)))
# 
# # Fixing sigmas for fishery catch and Fdevs here
# mapping$ln_sigmaC <- factor(rep(NA, length(parameters$ln_sigmaC)))
# 
# # Fixing continuous time-varying selecitvity stuff
# mapping$ln_fishsel_dev1 <- factor(rep(NA, length(parameters$ln_fishsel_dev1)))
# mapping$ln_fishsel_dev2 <- factor(rep(NA, length(parameters$ln_fishsel_dev2)))
# mapping$ln_fishsel_dev1_sd <- factor(rep(NA, length(parameters$ln_fishsel_dev1_sd)))
# mapping$ln_fishsel_dev2_sd <- factor(rep(NA, length(parameters$ln_fishsel_dev2_sd)))
# 
# # Fixing dirichlet mutlinomial stuff
# mapping$ln_FishAge_DM_theta <- factor(rep(NA, length(parameters$ln_FishAge_DM_theta)))
# mapping$ln_FishLen_DM_theta <- factor(rep(NA, length(parameters$ln_FishLen_DM_theta)))
# mapping$ln_SrvAge_DM_theta <- factor(rep(NA, length(parameters$ln_SrvAge_DM_theta)))
# mapping$ln_SrvLen_DM_theta <- factor(rep(NA, length(parameters$ln_SrvLen_DM_theta)))
# 
# # Movement
# mapping$move_pars <- factor(rep(NA, length(parameters$move_pars)))
# 
# data$srv_q_blocks = data$srv_q_blocks + 1
# data$fish_q_blocks = data$fish_q_blocks + 1
# data$fish_sel_blocks = data$fish_sel_blocks + 1
# data$srv_sel_blocks = data$srv_sel_blocks + 1
# data$bias_year = data$bias_year + 1
# data$sigmaR_switch = data$sigmaR_switch + 1
# 
# # make AD model function
# sabie_rtmb_model <- RTMB::MakeADFun(sabie_RTMB, parameters = parameters, map = mapping)
# 
# # Now, optimize the function
# sabie_optim <- stats::nlminb(sabie_rtmb_model$par, sabie_rtmb_model$fn, sabie_rtmb_model$gr,
#                              control = list(iter.max = 1e5, eval.max = 1e5, rel.tol = 1e-15))
# # newton steps
# try_improve <- tryCatch(expr =
#                           for(i in 1:4) {
#                             g = as.numeric(sabie_rtmb_model$gr(sabie_optim$par))
#                             h = optimHess(sabie_optim$par, fn = sabie_rtmb_model$fn, gr = sabie_rtmb_model$gr)
#                             sabie_optim$par = sabie_optim$par - solve(h,g)
#                             sabie_optim$objective = sabie_rtmb_model$fn(sabie_optim$par)
#                           }
#                         , error = function(e){e}, warning = function(w){w})
# 
# sabie_rtmb_model$optim <- sabie_optim # Save optimized model results
# sabie_rtmb_model$rep <- sabie_rtmb_model$report(sabie_rtmb_model$env$last.par.best) # Get report
# sabie_rtmb_model$sd_rep <- RTMB::sdreport(sabie_rtmb_model) # Get sd report
# 
# # Check consistency -------------------------------------------------------
# 
# # Get time series
# rec_series <- data.frame(Par = "Recruitment",
#                          Year = 1960:2021,
#                          TMB = sabie_rtmb_model$sd_rep$value[names(sabie_rtmb_model$sd_rep$value) == "Rec"],
#                          ADMB = rec[1:62])
# 
# f_series <- data.frame(Par = "Total F",
#                        Year = 1960:2021,
#                        TMB = apply(sabie_rtmb_model$rep$Fmort,2,sum),
#                        ADMB = tem_dat$t.series$fmort[1:62])
# 
# females_series <- data.frame(Par = "Total Females",
#                              Year = 1960:2021,
#                              TMB = rowSums(sabie_rtmb_model$rep$NAA[1,-63,,1]),
#                              ADMB = tem_dat$t.series$numbers.f[1:62])
# 
# males_series <- data.frame(Par = "Total Males",
#                            Year = 1960:2021,
#                            TMB = rowSums(sabie_rtmb_model$rep$NAA[1,-63,,2]),
#                            ADMB = tem_dat$t.series$numbers.m[1:62])
# 
# ssb_se_series <- data.frame(Par = "SSB (SE)",
#                             Year = 1960:2021,
#                             TMB = sabie_rtmb_model$sd_rep$sd[names(sabie_rtmb_model$sd_rep$value) == "SSB"],
#                             ADMB = tem_par$se[str_detect(names(tem_par$se), "ssb")][1:62])
# 
# ts_df <- rbind(ssb_se_series, rec_series, rec_se_series, f_series, females_series, males_series)
# 
# # Get selectivities
# dom_ll_fish_f1 <- data.frame(Age = 1:30,
#                              TMB = sabie_rtmb_model$rep$fish_sel[1,1,,1,1],
#                              ADMB = tem_dat$agesel$fish1sel.f,
#                              Type = "Domestic LL Fishery Female Block 1")
# 
# dom_ll_fish_m1 <- data.frame(Age = 1:30,
#                              TMB = sabie_rtmb_model$rep$fish_sel[1,1,,2,1],
#                              ADMB = tem_dat$agesel$fish1sel.m,
#                              Type = "Domestic LL Fishery Male Block 1")
# 
# dom_ll_fish_f2 <- data.frame(Age = 1:30,
#                              TMB = sabie_rtmb_model$rep$fish_sel[1,40,,1,1],
#                              ADMB = tem_dat$agesel$fish4sel.f,
#                              Type = "Domestic LL Fishery Female Block 2")
# 
# dom_ll_fish_m2 <- data.frame(Age = 1:30,
#                              TMB = sabie_rtmb_model$rep$fish_sel[1,40,,2,1],
#                              ADMB = tem_dat$agesel$fish4sel.m,
#                              Type = "Domestic LL Fishery Male Block 2")
# 
# dom_ll_fish_f3 <- data.frame(Age = 1:30,
#                              TMB = sabie_rtmb_model$rep$fish_sel[1,60,,1,1],
#                              ADMB = tem_dat$agesel$fish5sel.f,
#                              Type = "Domestic LL Fishery Female Block 3")
# 
# dom_ll_fish_m3 <- data.frame(Age = 1:30,
#                              TMB = sabie_rtmb_model$rep$fish_sel[1,60,,2,1],
#                              ADMB = tem_dat$agesel$fish5sel.m,
#                              Type = "Domestic LL Fishery Male Block 3")
# 
# dom_trwl_fish_f <- data.frame(Age = 1:30,
#                               TMB = sabie_rtmb_model$rep$fish_sel[1,1,,1,2],
#                               ADMB = tem_dat$agesel$fish3sel.f,
#                               Type = "Domestic Trawl Female")
# 
# dom_trwl_fish_m <- data.frame(Age = 1:30,
#                               TMB = sabie_rtmb_model$rep$fish_sel[1,1,,2,2],
#                               ADMB = tem_dat$agesel$fish3sel.m,
#                               Type = "Domestic Trawl Male")
# 
# dom_ll_srv_f1 <- data.frame(Age = 1:30,
#                             TMB = sabie_rtmb_model$rep$srv_sel[1,1,,1,1],
#                             ADMB = tem_dat$agesel$srv1sel.f,
#                             Type = "Domestic LL Survey Female Block 1")
# 
# dom_ll_srv_m1 <- data.frame(Age = 1:30,
#                             TMB = sabie_rtmb_model$rep$srv_sel[1,1,,2,1],
#                             ADMB = tem_dat$agesel$srv1sel.m,
#                             Type = "Domestic LL Survey Male Block 1")
# 
# dom_ll_srv_f2 <- data.frame(Age = 1:30,
#                             TMB = sabie_rtmb_model$rep$srv_sel[1,60,,1,1],
#                             ADMB = tem_dat$agesel$srv10sel.f,
#                             Type = "Domestic LL Survey Female Block 2")
# 
# dom_ll_srv_m2 <- data.frame(Age = 1:30,
#                             TMB = sabie_rtmb_model$rep$srv_sel[1,60,,2,1],
#                             ADMB = tem_dat$agesel$srv10sel.m,
#                             Type = "Domestic LL Survey Male Block 2")
# 
# dom_trwl_srv_f2 <- data.frame(Age = 1:30,
#                               TMB = sabie_rtmb_model$rep$srv_sel[1,60,,1,2],
#                               ADMB = tem_dat$agesel$srv7sel.f,
#                               Type = "Domestic Trawl Survey Female")
# 
# dom_trwl_srv_m2 <- data.frame(Age = 1:30,
#                               TMB = sabie_rtmb_model$rep$srv_sel[1,60,,2,2],
#                               ADMB = tem_dat$agesel$srv7sel.m,
#                               Type = "Domestic Trawl Survey Male")
# 
# coop_ll_srv_f2 <- data.frame(Age = 1:30,
#                              TMB = sabie_rtmb_model$rep$srv_sel[1,60,,1,3],
#                              ADMB = tem_dat$agesel$srv2sel.f,
#                              Type = "Coop LL Survey Female")
# 
# coop_ll_srv_m2 <- data.frame(Age = 1:30,
#                              TMB = sabie_rtmb_model$rep$srv_sel[1,60,,2,3],
#                              ADMB = tem_dat$agesel$srv2sel.m,
#                              Type = "Coop LL Survey Male")
# 
# combined_sel <- rbind(
#   dom_ll_fish_m1,
#   dom_ll_fish_f2,
#   dom_ll_fish_m2,
#   dom_ll_fish_f3,
#   dom_ll_fish_m3,
#   dom_trwl_fish_f,
#   dom_trwl_fish_m,
#   dom_ll_srv_f1,
#   dom_ll_srv_m1,
#   dom_ll_srv_f2,
#   dom_ll_srv_m2,
#   dom_trwl_srv_f2,
#   dom_trwl_srv_m2,
#   coop_ll_srv_f2,
#   coop_ll_srv_m2
# )
# 
# # Plots -------------------------------------------------------------------
# 
# # Time Series Estimated
# ggplot() +
#   geom_line(ts_df, mapping  = aes(x = Year, y = TMB, color = "RTMB"), size = 1.3, lty = 1) +
#   geom_line(ts_df, mapping  = aes(x = Year, y = ADMB, color = "ADMB"), size = 1.3, lty = 2) +
#   facet_wrap(~Par, scales = "free") +
#   labs(x = "Year", color = 'Model', y = "Value") +
#   ggthemes::scale_color_hc() +
#   theme_sablefish()
# 
# 
# # Selectivity curves
# ggplot() +
#   geom_line(combined_sel, mapping  = aes(x = Age, y = TMB, color = "RTMB"), size = 1.3, lty = 1) +
#   geom_line(combined_sel, mapping  = aes(x = Age, y = ADMB, color = "ADMB"), size = 1.3, lty = 2) +
#   ggthemes::scale_color_hc() +
#   facet_wrap(~Type) +
#   labs(y = "Selex", color = "Model") +
#   theme_sablefish()

