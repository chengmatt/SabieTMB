# Purpose: To bridge to the spatial assessment in Cheng and Marsh et al. 2025 
# Creator: Matthew LH. Cheng
# Date Created: 2/5/25

# Set up ------------------------------------------------------------------

library(here)
library(tidyverse)
library(RTMB)

source(here("R", "functions", "Francis_Reweight.R"))
source(here("R", "model", "v3", "Sabie_RTMB_v3.R"))
source(here("R", "functions", "Utility_Functions.R"))

# Read in spatial datas
spatial_data <- readRDS(here("dev", "Spatial Sablefish Model", "data.RDS"))
spatial_pars <- readRDS(here("dev", "Spatial Sablefish Model", "parameters.RDS")) 
spatial_rep <- readRDS(here("dev", "Spatial Sablefish Model", "mle_report.RDS")) 
spatial_map <- readRDS(here("dev", "Spatial Sablefish Model", "map_fixed_pars.RDS")) 

ageing_dat <- dget(here("dev",'2023 Base (23.5)_final model', 'test.rdat')) # for getting ageing error
tag_rel <- readRDS(here("dev", "Spatial Sablefish Model", "Tag_release_summarised.RDS")) 
tag_rec <- readRDS(here("dev", "Spatial Sablefish Model", "Tag_recovery_summarised.RDS")) 


# Prepare Data ------------------------------------------------------------
data <- list() # make data list

# Set up dimensions
data$n_regions <- spatial_data$n_regions # number of regions
data$ages <- spatial_data$ages # ages
data$lens <- spatial_data$length_bins # lengths
data$years <- spatial_data$years # years
data$n_sexes <- 2 # number of sexes
data$n_fish_fleets <- 2 # number of fishery fleets (0 == fixed gear, 1 == trawl gear)
data$n_srv_fleets <- 2 # number of survey fleets

# Recruitment and population initializaiton stuff
data$init_F_prop <- 0 # initial F proportion for initializing population
data$do_rec_bias_ramp <- 0 # do bias ramp (slot 0 == don't do bias ramp, 1 == do bias ramp)
data$bias_year <- NA
data$sigmaR_switch <- 1 # when to switch to late sigmaR
data$sexratio <- c(0.5, 0.5) # recruitment sex ratio (assuming 1, since single sex)
data$init_age_strc <- 0 # iterative approach to calculate initial age structure

# Movement stuff
data$do_recruits_move <- 0 # recruits dont move
data$use_fixed_movement <- 1 # use fixed movement
data$Fixed_Movement <- array(0, dim = c(data$n_regions, data$n_regions, length(data$years), length(data$ages), data$n_sexes))
data$Fixed_Movement[,,,1:6,] <- spatial_rep$movement_matrix[,,,1] # age block 1
data$Fixed_Movement[,,,7:15,] <- spatial_rep$movement_matrix[,,,2] # age block 2
data$Fixed_Movement[,,,16:30,] <- spatial_rep$movement_matrix[,,,3] # age block 3

# Set up movement blocks
data$move_age_blocks <- list(c(1:6), c(7:15), c(16:30))
# data$move_age_blocks <- list(1:30)
data$move_sex_blocks <- list(c(1:2))

# Likelihood and data weighting stuff
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
data$WAA <- array(0, dim = c(data$n_regions,length(data$years),length(data$ages), data$n_sexes))

# Maturity at age
data$MatAA <- array(0, dim = c(data$n_regions,length(data$years),length(data$ages), data$n_sexes))

for(r in 1:data$n_regions) {
  data$WAA[r,,,1] <- t(spatial_data$female_mean_weight_by_age) # female weight at age
  data$WAA[r,,,2] <- t(spatial_data$male_mean_weight_by_age) # male weight at age
  # maturity
  data$MatAA[r,,,1] <- t(spatial_data$maturity)
  data$MatAA[r,,,2] <- t(spatial_data$maturity)
}

# Ageing Error
data$AgeingError <- as.matrix(ageing_dat$age_error) # ageing error matrix

# Size Transition Matrix
data$fit_lengths <- 0
data$SizeAgeTrans <- NA
data$SizeAgeTrans <- array(0, dim = c(data$n_regions, length(data$years), length(data$lens), length(data$ages), data$n_sexes)) # size age transition matrix
for(y in 1:length(data$years)) {
  for(r in 1:data$n_regions) {
    data$SizeAgeTrans[r,y,,,1] <- t(spatial_data$female_age_length_transition[,,y])
    data$SizeAgeTrans[r,y,,,1] <- apply(data$SizeAgeTrans[r,y,,,1], 2, function(x) x / sum(x))
    data$SizeAgeTrans[r,y,,,2] <- t(spatial_data$male_age_length_transition[,,y]) / colSums(t(spatial_data$male_age_length_transition[,,y]))
    data$SizeAgeTrans[r,y,,,2] <- apply(data$SizeAgeTrans[r,y,,,2], 2, function(x) x / sum(x))
  } # end r loop
} # end y loop

sum(t(spatial_data$female_age_length_transition[,,y])[,3] / colSums(t(spatial_data$female_age_length_transition[,,y]))[3])
sum(data$SizeAgeTrans[r,y,,3,2])

# Fishery Observations
data$Use_F_pen <- 1
data$ObsCatch <- array(0, dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
data$ObsCatch[,,1] <- spatial_data$fixed_fishery_catch # fixed gear catch
data$ObsCatch[,,2] <- spatial_data$trwl_fishery_catch # trawl gear catch
data$ObsCatch[,1:3,2] <- 0 # input 0 for first three years
data$UseCatch <- array(1, c(data$n_regions, length(data$years), data$n_fish_fleets)) # fit catch data everywhere
data$UseCatch[,1:3,2] <- 0
data$Catch_Type <- array(1, dim = c(length(data$years), data$n_fish_fleets)) # regional catch is availiable
data$est_all_regional_F <- 1 # estimate all regional Fs
data$Catch_Constant <- rep(0, data$n_fish_fleets) # don't add any constants

# Fishery Indices (Not fit to)
data$ObsFishIdx <- array(NA, c(data$n_regions, length(data$years), data$n_fish_fleets))
data$ObsFishIdx_SE <- array(NA, c(data$n_regions, length(data$years), data$n_fish_fleets))
data$UseFishIdx <- array(0, c(data$n_regions, length(data$years), data$n_fish_fleets))
colnames(data$ObsFishIdx) <- data$years # define row years
colnames(data$ObsFishIdx_SE) <- data$years # define row years

# Not using single area ADMB aggregation methods and calculations
data$sablefish_ADMB = 0
data$FishAge_comp_agg_type = NA
data$FishLen_comp_agg_type = NA
data$SrvAge_comp_agg_type = NA
data$SrvLen_comp_agg_type = NA

# Fishery Age Compositions (Joint by sex, split by region)
data$ObsFishAgeComps <- array(0, dim = c(data$n_regions, length(data$years), length(data$ages), data$n_sexes, data$n_fish_fleets))
for(y in 1:length(data$years)) {
  for(r in 1:data$n_regions) {
    if(spatial_data$fixed_catchatage_indicator[r,y] == 1) {
      data$ObsFishAgeComps[r,y,,1,1] <- spatial_data$obs_fixed_catchatage[-(1:length(data$ages)),r,y] / sum(spatial_data$obs_fixed_catchatage[,r,y]) # females
      data$ObsFishAgeComps[r,y,,2,1] <- spatial_data$obs_fixed_catchatage[1:length(data$ages),r,y] / sum(spatial_data$obs_fixed_catchatage[,r,y]) # males
    } # end if
  } # end r loop
} # end y loop

colnames(data$ObsFishAgeComps) <- data$years # define row years
data$UseFishAgeComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
data$UseFishAgeComps[,,1] <- spatial_data$fixed_catchatage_indicator # fixed gear use indicator
colnames(data$UseFishAgeComps) <- data$years # define row years

data$ISS_FishAgeComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
data$ISS_FishAgeComps[,,1,1] <- apply(spatial_data$obs_fixed_catchatage, c(2,3), sum) # ISS for fixed gear fishery
colnames(data$ISS_FishAgeComps) <- data$years # define row years

data$Wt_FishAgeComps <- array(NA, dim = c(data$n_regions, data$n_sexes, data$n_fish_fleets)) # weights for fishery age comps
data$Wt_FishAgeComps[,1,1] <- 1 # Weight for fixed gear age comps

# Fishery Length Comps (Joint by sex, split by region)
data$ObsFishLenComps <- array(0, dim = c(data$n_regions, length(data$years), length(data$lens), data$n_sexes, data$n_fish_fleets))
for(y in 1:length(data$years)) {
  for(r in 1:data$n_regions) {
    
    # fixed gear
    if(spatial_data$fixed_catchatlgth_indicator[r,y] == 1) {
      data$ObsFishLenComps[r,y,,1,1] <- spatial_data$obs_fixed_catchatlgth[-(1:length(data$ages)),r,y] / sum(spatial_data$obs_fixed_catchatlgth[,r,y]) # females
      data$ObsFishLenComps[r,y,,2,1] <- spatial_data$obs_fixed_catchatlgth[1:length(data$ages),r,y] / sum(spatial_data$obs_fixed_catchatlgth[,r,y]) # males
    } # end if
    
    # trawl gear
    if(spatial_data$trwl_catchatlgth_indicator[r,y] == 1) {
      data$ObsFishLenComps[r,y,,1,2] <- spatial_data$obs_trwl_catchatlgth[-(1:length(data$ages)),r,y] / sum(spatial_data$obs_trwl_catchatlgth[,r,y]) # females
      data$ObsFishLenComps[r,y,,2,2] <- spatial_data$obs_trwl_catchatlgth[1:length(data$ages),r,y] / sum(spatial_data$obs_trwl_catchatlgth[,r,y]) # males
    } # end if
    
  } # end r loop
} # end y loop

colnames(data$ObsFishLenComps) <- data$years # define row years
data$UseFishLenComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
data$UseFishLenComps[,,1] <- spatial_data$fixed_catchatlgth_indicator # fixed gear use indicator
data$UseFishLenComps[,,2] <- spatial_data$trwl_catchatlgth_indicator # trawl gear use indicator
colnames(data$UseFishLenComps) <- data$years # define row years

data$ISS_FishLenComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
data$ISS_FishLenComps[,,1,1] <- apply(spatial_data$obs_fixed_catchatlgth, c(2,3), sum) # ISS for fixed gear fishery
data$ISS_FishLenComps[,,1,2] <- apply(spatial_data$obs_trwl_catchatlgth, c(2,3), sum) # ISS for trawl gear fishery
colnames(data$ISS_FishLenComps) <- data$years # define row years

data$Wt_FishLenComps <- array(NA, dim = c(data$n_regions, data$n_sexes, data$n_fish_fleets)) # weights for fish len comps
data$Wt_FishLenComps[,1,] <- 1 # Weight for fish len comps

# Composition munging stuff
data$FishAgeComps_LikeType <- array(0, dim = c(data$n_fish_fleets)) # multinomial for both fleets
data$FishLenComps_LikeType <- array(0, dim = c(data$n_fish_fleets)) # multinomial for both fleets
data$FishAgeComps_Type <- array(2, dim = c(length(data$years), data$n_fish_fleets)) # joint by sex, split by region
data$FishLenComps_Type <- array(2, dim = c(length(data$years), data$n_fish_fleets)) # joint by sex, split by region

# Survey Observations (== 1, coop, == 2, domestic)
# Survey Indices
data$ObsSrvIdx <- array(0, dim = c(data$n_regions, length(data$years), data$n_srv_fleets))
data$ObsSrvIdx <- spatial_data$obs_srv_bio
data$ObsSrvIdx_SE <- spatial_data$obs_srv_se

data$UseSrvIdx <- array(0, dim = c(data$n_regions, length(data$years), data$n_srv_fleets))
data$UseSrvIdx <- spatial_data$srv_bio_indicator # survey indicator
colnames(data$ObsSrvIdx) <- data$years # define row years
colnames(data$ObsSrvIdx_SE) <- data$years # define row years

# Survey Age Comps
data$ObsSrvAgeComps <- array(0, dim = c(data$n_regions, length(data$years), length(data$ages), data$n_sexes, data$n_srv_fleets))

for(y in 1:length(data$years)) {
  for(r in 1:data$n_regions) {
    
    if(spatial_data$srv_catchatage_indicator[r,y,1] == 1) { # coop survey
      data$ObsSrvAgeComps[r,y,,1,1] <- spatial_data$obs_srv_catchatage[-(1:length(data$ages)),r,y,1] / sum(spatial_data$obs_srv_catchatage[,r,y,1]) # females
      data$ObsSrvAgeComps[r,y,,2,1] <- spatial_data$obs_srv_catchatage[1:length(data$ages),r,y,1] / sum(spatial_data$obs_srv_catchatage[,r,y,1]) # males
    } # end if
    
    if(spatial_data$srv_catchatage_indicator[r,y,2] == 1) { # domestic survey
      data$ObsSrvAgeComps[r,y,,1,2] <- spatial_data$obs_srv_catchatage[-(1:length(data$ages)),r,y,2] / sum(spatial_data$obs_srv_catchatage[,r,y,2]) # females
      data$ObsSrvAgeComps[r,y,,2,2] <- spatial_data$obs_srv_catchatage[1:length(data$ages),r,y,2] / sum(spatial_data$obs_srv_catchatage[,r,y,2]) # males
    } # end if
    
  } # end r loop
} # end y loop

colnames(data$ObsSrvAgeComps) <- data$years # define row years
data$UseSrvAgeComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_srv_fleets))
data$UseSrvAgeComps <- spatial_data$srv_catchatage_indicator

# Data weighting for survey age compositions
data$ISS_SrvAgeComps <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_srv_fleets))
data$ISS_SrvAgeComps[,,1,1] <- apply(spatial_data$obs_srv_catchatage[,,,1], c(2,3), sum) # ISS for coop survey
data$ISS_SrvAgeComps[,,1,2] <- apply(spatial_data$obs_srv_catchatage[,,,2], c(2,3), sum) # ISS for domestic survey
colnames(data$ISS_SrvAgeComps) <- data$years # define row years
data$Wt_SrvAgeComps <- array(NA, dim = c(data$n_regions, data$n_sexes, data$n_srv_fleets)) # weights for survey age comps
data$Wt_SrvAgeComps[,1,] <- 1 # Weight for surveys

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
data$SrvAgeComps_Type <- array(2, dim = c(length(data$years), data$n_srv_fleets)) # joint by sex, split by region
data$SrvLenComps_Type <- array(2, dim = c(length(data$years), data$n_srv_fleets)) # joint by sex, split by region

# Fishery Selectivity
data$cont_tv_fish_sel <- array(0, dim = c(data$n_regions, data$n_fish_fleets)) # no timevarying selex continously
# Time Block Specification
data$fish_sel_blocks <- array(NA, dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
for(r in 1:data$n_regions) data$fish_sel_blocks[r,,1] <- spatial_data$fixed_sel_by_year_indicator + 1 # fishery blocks for fixed gear
data$fish_sel_blocks[,,2] <- 1 # fishery blocks for trawl gear

# Selectivity Model
data$fish_sel_model <- array(NA, dim = c(data$n_regions, length(data$years), data$n_fish_fleets))
data$fish_sel_model[,,1] <- 0 # Logistic selectivity for fixed gear
data$fish_sel_model[,,2] <- 1 # Gamma selectivity for trawl gear

# Catchability (not used for fishery Q)
# Time Block Specification
data$fish_q_blocks <- data$fish_sel_blocks # catchability blocks are same as selectivity blocks
data$fish_idx_type <- array(1, dim = c(data$n_regions, data$n_fish_fleets)) # fishery index type (0 == abundance, 1 == biomass)

# Survey Selectivity
# Time Block Specification
data$srv_sel_blocks <- array(1, dim = c(data$n_regions, length(data$years), data$n_srv_fleets)) # no survey blocks
# Selectivity Model
data$srv_sel_model <- array(0, dim = c(data$n_regions, length(data$years), data$n_srv_fleets)) # logistic selectivity

# Time Block Specification
data$srv_q_blocks <- array(1, dim = c(data$n_regions, length(data$years), data$n_srv_fleets)) # catchability blocks are same as selectivity blocks
data$srv_idx_type <- array(0, dim = c(data$n_regions, data$n_srv_fleets)) # abundance based for both surveys

# Tagging and Movement Stuff
data$UseTagging <- 1 # use tagging data

# Figure out tag release indicator
# name regions correctly
tag_rel <- tag_rel %>% mutate(region_num = case_when(
  region_release == "BS" ~ 1,
  region_release == "AI" ~ 2,
  region_release == "WGOA" ~ 3,
  region_release == "CGOA" ~ 4,
  region_release == "EGOA" ~ 5
)) %>% filter(release_year %in% c(data$years))

# name regions correctly
tag_rec <- tag_rec %>% mutate(region_num = case_when(
  recovery_region == "BS" ~ 1,
  recovery_region == "AI" ~ 2,
  recovery_region == "WGOA" ~ 3,
  recovery_region == "CGOA" ~ 4,
  recovery_region == "EGOA" ~ 5
)) %>% filter(recovery_year %in% c(data$years)) %>%
  mutate(tag_lib = (recovery_year - release_year) + 1)

# loop through to get releases
rel_cohorts_df <- tag_rel %>% distinct(release_event_id, release_year, region_release, region_num)
tag_release_ind <- data.frame()
for(i in 1:nrow(rel_cohorts_df)) {
  tmp <- data.frame(regions = rel_cohorts_df[i,4], tag_yrs = rel_cohorts_df[i,2]-1959)
  tag_release_ind <- rbind(tag_release_ind, tmp)
} # end i

data$tag_release_indicator <- as.matrix(tag_release_ind) # input releases of cohorts in
data$n_tag_cohorts <- nrow(data$tag_release_indicator) # number of tag cohorts
data$max_tag_liberty <- 15 # maximum liberty to track cohorts

# Set up tagged fish
data$Tagged_Fish <- array(0, dim = c(data$n_tag_cohorts, length(data$ages), data$n_sexes)) # tagged fish
for(i in 1:nrow(data$tag_release_indicator )) {
  tmp_rel_f <- tag_rel %>% filter(release_event_id == i, sex == "F") # filter to a given cohort females
  tmp_rel_m <- tag_rel %>% filter(release_event_id == i, sex == "M") # filter to a given cohort males
  data$Tagged_Fish[i,,1] <- tmp_rel_f$Nage_at_release # females
  data$Tagged_Fish[i,,2] <- tmp_rel_m$Nage_at_release # males
} # end i for tag cohorts

data$Obs_Tag_Recap <- array(0, dim = c(data$max_tag_liberty, data$n_tag_cohorts, data$n_regions, length(data$ages), data$n_sexes))

# pre filter data frame before looping
tag_rec_f <- tag_rec %>% filter(sex == "F")
tag_rec_m <- tag_rec %>% filter(sex == "M")

for(ry in 1:data$max_tag_liberty) {
  for(i in 1:data$n_tag_cohorts) {
    for(r in 1:data$n_regions) {
      tmp_rec_f <- tag_rec_f %>% filter(release_event_id == i, region_num == r, tag_lib == ry)
      tmp_rec_m <- tag_rec_m %>% filter(release_event_id == i, region_num == r, tag_lib == ry)
      if(nrow(tmp_rec_f) > 0) data$Obs_Tag_Recap[ry,i,r,,1] <- tmp_rec_f$Nage_at_recovery
      if(nrow(tmp_rec_m) > 0) data$Obs_Tag_Recap[ry,i,r,,2] <- tmp_rec_m$Nage_at_recovery
    }
  }
}

data$Tag_LikeType <- 1 # poisson likelihood
data$mixing_period <- 3 # when to start mixing period after first release year
data$t_tagging <- 0.5 # discounting for tagging

# tag reporting rate priors
data$Use_TagRep_Prior = 1 # use tag reporting rate prior
data$TagRep_PriorType = 0 # symmetric beta prior with upper and lower bounds at 1 and 0 (away from edges)
data$TagRep_mu = NA # penalize mu
data$TagRep_sd = 10 # sd of tag reporting rate prio

# movement rate priors
data$Use_Movement_Prior = 1 # use mvoement rate prior
data$Movement_prior = array(1.2, dim = c(data$n_regions, data$n_regions, length(data$years), length(data$ages), data$n_sexes)) # away form edges

# Recruitment
data$rec_model <- 0 # mean recruitment
data$rec_lag <- 1 # recruitment ssb lag
data$Use_h_prior <- 0 # don't use steepness prior
data$h_mu <- NA
data$h_sd <- NA


# Parameters --------------------------------------------------------------
parameters <- list()
parameters$dummy <- 1

# Fishing Mortality
parameters$ln_sigmaC <- array(log(0.02), dim = c(data$n_regions, data$n_fish_fleets))
parameters$ln_F_mean <- array(0, dim = c(data$n_regions, data$n_fish_fleets)) # mean fishing mortality
parameters$ln_F_mean[,1] <- log(rowMeans(spatial_rep$annual_F_fixed)) # fixed gear
parameters$ln_F_mean[,2] <- log(rowMeans(spatial_rep$annual_F_trwl)) # trawl gear

parameters$ln_F_devs <- array(0, dim = c(data$n_regions, length(data$years), data$n_fish_fleets)) # fishing mortality deviations from mean
parameters$ln_F_devs[,,1] <- log(spatial_rep$annual_F_fixed) - parameters$ln_F_mean[,1] # fixed gear
parameters$ln_F_devs[,,2] <- log(spatial_rep$annual_F_trwl) - parameters$ln_F_mean[,2] # trawl gear

# Aggregated fishing mortality (not used)
parameters$ln_F_devs_AggCatch = array(0, dim = c(sum(data$Catch_Type == 0), data$n_fish_fleets))
parameters$ln_F_mean_AggCatch = 0

# Set up continuous fishery selectivity stuff (not used)
parameters$ln_fishsel_dev1 <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
parameters$ln_fishsel_dev2 <- array(0, dim = c(data$n_regions, length(data$years), data$n_sexes, data$n_fish_fleets))
parameters$ln_fishsel_dev1_sd <- array(0, dim = c(data$n_regions, data$n_sexes, data$n_fish_fleets))
parameters$ln_fishsel_dev2_sd <- array(0, dim = c(data$n_regions, data$n_sexes, data$n_fish_fleets))

# Fixed Gear Fishery three time blocks
max_fish_blks <- 2 # maximum number of fishery blocks for any fleet
max_fish_pars <- 2 # maximum number of fishery fixed parameters for any fleet
parameters$ln_fish_fixed_sel_pars <- array(0, dim = c(data$n_regions, max_fish_pars, max_fish_blks, data$n_sexes, data$n_fish_fleets))

for(r in 1:data$n_regions) {
  parameters$ln_fish_fixed_sel_pars[r,,,,1] <- log(spatial_rep$fixed_sel_pars) # fixed gear pars
  parameters$ln_fish_fixed_sel_pars[r,,1:2,,2] <- log(spatial_rep$trwl_sel_pars) # trawl gear pars
}

# Fishery Catchability
# Fixed Gear Fishery Catchability
parameters$ln_fish_q <- array(0, dim = c(data$n_regions, max_fish_blks, data$n_fish_fleets)) # not used

# Survey Stuff
max_srv_blks <- 1 # maximum number of survey blocks for any fleet
max_srv_pars <- 2 # maximum number of survey fixed parameters for any fleet
parameters$ln_srv_fixed_sel_pars <- array(log(3), dim = c(data$n_regions, max_srv_pars, max_srv_blks, data$n_sexes, data$n_srv_fleets))

# Survey Catchability
max_q_srv_blks <- 1 # maximum catchability blocks
# Longline Survey Catchability
parameters$ln_srv_q <- array(0, dim = c(data$n_regions, max_q_srv_blks, data$n_srv_fleets))
parameters$ln_srv_q <- log(spatial_rep$srv_q)

# Natural Mortality
parameters$ln_M <- log(0.104884) # Female M (base)
parameters$M_offset <- 0 # M offset not used

# Recruitment 
parameters$ln_global_R0 <- log(sum(spatial_rep$mean_rec)) # mean recruitment
parameters$R0_prop <- array(log(spatial_rep$mean_rec[-1]) - log(spatial_rep$mean_rec[1]) , dim = c(data$n_regions - 1))

parameters$h <- c(10, 10) # not used
parameters$ln_InitDevs <- array(0, dim = c(data$n_regions, length(data$ages) - 2))
for(i in 1:data$n_regions) parameters$ln_InitDevs[i,] <- spatial_rep$init_rec_dev
# don't estimate last 2 year rec devs
parameters$ln_RecDevs <- array(0, dim = c(data$n_regions, length(data$years) - 1))
parameters$ln_RecDevs <- spatial_rep$recruitment_devs[,-c(62)]
parameters$ln_sigmaR_early <- log(0.4) # early sigma R
parameters$ln_sigmaR_late <- log(1.2)  # late sigma R

# Comp likelihood stuff (not used)
parameters$ln_FishAge_theta <- array(0, dim = c(data$n_regions, data$n_fish_fleets))
parameters$ln_FishLen_theta <- array(0, dim = c(data$n_regions, data$n_fish_fleets))
parameters$ln_SrvAge_theta <- array(0, dim = c(data$n_regions, data$n_srv_fleets))
parameters$ln_SrvLen_theta <- array(0, dim = c(data$n_regions, data$n_srv_fleets))

# Movement stuff
parameters$move_pars <- array(0, dim = c(data$n_regions, data$n_regions - 1, length(data$years), length(data$ages), data$n_sexes))

for(i in 1:data$n_regions) {
  for(t in 1:length(data$years)) {
    for(a in 1:length(data$ages)) {
      for(s in 1:data$n_sexes) {
        parameters$move_pars[,,t,a,s] = SpatialSablefishAssessment::simplex(data$Fixed_Movement[i,,t,a,s])
      }
    } # end a block
  } # end t block
} # end i region

# Tagging Stuff
parameters$ln_Init_Tag_Mort <- log(unique(spatial_data$initial_tag_induced_mortality)) # initial tag induced mortality (not estimated)
parameters$ln_Tag_Shed <- log(spatial_data$annual_tag_shedding_rate) # annual tag shedding (not estimated)
parameters$Tag_Reporting_Pars <- array(log(0.2 / (1 - 0.2)), dim = c(data$n_regions, length(data$years)))
parameters$ln_tag_theta <- log(0.1)

# Mapping -----------------------------------------------------------------
mapping <- list()
mapping$dummy <- factor(NA)
mapping$ln_sigmaR_late <- factor(NA)
mapping$ln_sigmaR_early <- factor(NA) # fix early sigma R
mapping$ln_M <- factor(NA) # fix natural mortality 
mapping$M_offset <- factor(NA) # fix natural mortality offset
mapping$ln_fish_q <- factor(rep(NA, length(parameters$ln_fish_q))) # no fishery q estimated

# Fixing sigmas for fishery catch 
mapping$ln_sigmaC <- factor(rep(NA, length(parameters$ln_sigmaC)))
mapping$ln_F_devs_AggCatch <- factor(rep(NA, data$n_fish_fleets * sum(data$Catch_Type == 0)))
mapping$ln_F_mean_AggCatch <- factor(rep(NA, length(data$n_fish_fleets)))

map_F_devs <- parameters$ln_F_devs
map_F_devs[] <- 1:length(map_F_devs)
map_F_devs[,1:3,2] <- NA # input NA for trawl devs
map_F_devs[,,2] <- map_F_devs[,,2] - length(map_F_devs[,1:3,2]) # subtract the number of NAs out
mapping$ln_F_devs <- factor(map_F_devs)

# Fishery selectivity stuff
map_ln_fish_fixed_sel_pars <- parameters$ln_fish_fixed_sel_pars # mapping fishery selectivity

# Fixed gear fleet, unique parameters for each sex (time block 1)
map_ln_fish_fixed_sel_pars[,1,1,1,1] <- 1 # a50, female, time block 1, fixed, gear
map_ln_fish_fixed_sel_pars[,2,1,1,1] <- 2 # delta, female, time block 1, fixed, gear
map_ln_fish_fixed_sel_pars[,1,1,2,1] <- 3 # a50, male, time block 1, fixed, gear
map_ln_fish_fixed_sel_pars[,2,1,2,1] <- 2 # delta, male, time block 1, fixed, gear

# time block 2, fixed gear fishery
map_ln_fish_fixed_sel_pars[,1,2,1,1] <- 4 # a50, female, time block 2, fixed, gear
map_ln_fish_fixed_sel_pars[,2,2,1,1] <- 2 # delta, female, time block 2, fixed, gear
map_ln_fish_fixed_sel_pars[,1,2,2,1] <- 5 # a50, male, time block 2, fixed, gear
map_ln_fish_fixed_sel_pars[,2,2,2,1] <- 2 # delta, male, time block 2, fixed, gear

# time block 1 and 2, trawl gear fishery
map_ln_fish_fixed_sel_pars[,1,,1,2] <- 6 # amax, female, time block 1, trawl, gear
map_ln_fish_fixed_sel_pars[,2,,1,2] <- 7 # delta, female, time block 1, trawl, gear
map_ln_fish_fixed_sel_pars[,1,,2,2] <- 8 # amax, male, time block 1, trawl, gear
map_ln_fish_fixed_sel_pars[,2,,2,2] <- 7 # delta, male, time block 1, trawl, gear

# fix fishery selectivity
mapping$ln_fish_fixed_sel_pars <- factor(as.vector(map_ln_fish_fixed_sel_pars)) 

# Fixing continuous time-varying selecitvity stuff (not used)
mapping$ln_fishsel_dev1 <- factor(rep(NA, length(parameters$ln_fishsel_dev1)))
mapping$ln_fishsel_dev2 <- factor(rep(NA, length(parameters$ln_fishsel_dev2)))
mapping$ln_fishsel_dev1_sd <- factor(rep(NA, length(parameters$ln_fishsel_dev1_sd)))
mapping$ln_fishsel_dev2_sd <- factor(rep(NA, length(parameters$ln_fishsel_dev2_sd)))

# Fix survey selectivity stuff
map_ln_srv_fixed_sel_pars <- parameters$ln_srv_fixed_sel_pars # set up mapping factor stuff

# Coop survey
map_ln_srv_fixed_sel_pars[,1,1,1,1] <- 1 # a50, coop survey, time block 1, female
map_ln_srv_fixed_sel_pars[,2,1,1,1] <- 2 # delta, coop survey, time block 1, female
map_ln_srv_fixed_sel_pars[,1,1,2,1] <- 3 # a50, coop survey, time block 1, male
map_ln_srv_fixed_sel_pars[,2,1,2,1] <- 4 # delta, coop survey, time block 1, male

# domestic survey
map_ln_srv_fixed_sel_pars[,1,1,1,2] <- 5 # a50, domestic survey, time block 1, female
map_ln_srv_fixed_sel_pars[,2,1,1,2] <- 2 # delta, domestic survey, time block 1, female
map_ln_srv_fixed_sel_pars[,1,1,2,2] <- 6 # a50, domestic survey, time block 1, male
map_ln_srv_fixed_sel_pars[,2,1,2,2] <- 4 # delta, domestic survey, time block 1, male

mapping$ln_srv_fixed_sel_pars <- factor(as.vector(map_ln_srv_fixed_sel_pars)) # fix survey selectivity

# Fixing dirichlet mutlinomial stuff (not used)
mapping$ln_FishAge_theta <- factor(rep(NA, length(parameters$ln_FishAge_theta)))
mapping$ln_FishLen_theta <- factor(rep(NA, length(parameters$ln_FishLen_theta)))
mapping$ln_SrvAge_theta <- factor(rep(NA, length(parameters$ln_SrvAge_theta)))
mapping$ln_SrvLen_theta <- factor(rep(NA, length(parameters$ln_SrvLen_theta)))

# Fixing survey catchability stuff
map_ln_srv_q <- parameters$ln_srv_q
map_ln_srv_q[,,1] <- 1 # coop survey
map_ln_srv_q[,,2] <- 2 # domestic survey
mapping$ln_srv_q <- factor(as.vector(map_ln_srv_q))

# Fixing Tagging Stuff
mapping$ln_Init_Tag_Mort <- factor(NA)
mapping$ln_Tag_Shed <- factor(NA)

# Tag reporting rates
map_Tag_Reporting_Pars <- parameters$Tag_Reporting_Pars
map_Tag_Reporting_Pars[,(1960:1994)-1959] <- 1 # pre-IFQ
map_Tag_Reporting_Pars[,(1995:2016)-1959] <- 2 # post-IFQ
map_Tag_Reporting_Pars[,(2017:2021)-1959] <- 2 # pot fleet transition (based on earlier analyses, probably no need to acocunt for that)

# map these off
mapping$Tag_Reporting_Pars <- factor(as.vector(map_Tag_Reporting_Pars))
data$map_Tag_Reporting_Pars <- array(as.numeric(mapping$Tag_Reporting_Pars), dim = dim(parameters$Tag_Reporting_Pars))
# mapping$Tag_Reporting_Pars <- factor(rep(NA, length(parameters$Tag_Reporting_Pars)))

# map theta off
# mapping$ln_tag_theta <- factor(NA)

# Movement
map_Movement_Pars <- parameters$move_pars 
map_Movement_Pars[,,1:length(data$years),1:6,] <- 1:20 # age block 1-6 for all years and sexes
map_Movement_Pars[,,1:length(data$years),7:15,] <- 21:40 # age block 7-15 for all years and sexes
map_Movement_Pars[,,1:length(data$years),16:30,] <- 41:60 # age block 16-30 for all years and sexes
mapping$move_pars <- factor(map_Movement_Pars)
# Currently not estimating movement at all
# mapping$move_pars <- factor(rep(NA, length(parameters$move_pars)))
data$map_Movement_Pars <- array(as.numeric(mapping$move_pars), dim = dim(parameters$move_pars))

# Steepness (not used)
mapping$h <- factor(rep(NA, length(parameters$h)))
data$map_h_Pars <- as.numeric(mapping$h)

# Recruitment stuff
map_init_devs <- parameters$ln_InitDevs
for(r in 1:data$n_regions) map_init_devs[r,] <- 1:length(map_init_devs[r,])
mapping$ln_InitDevs <- factor(as.vector(map_init_devs))

# Change Log -------------------------------------------------------------------
# This section is when no tagging data are used and movement is estimated with a movement prior
# to push estimation away from edges
# Getting trawl selectivity to estimate properly ....
# Increasing ISS and inputing a recruitment bias ramp helps lock down those
# early recruiments in the BS and helps with estimation ... 

# Exploration 1
# Changing ISS to generally reflect ssampling imprecision
# Lowest for Trawl fishery lengths, followed by fishery ages, and Survey ages 
data$ISS_FishLenComps[,,,1] <- 40 # fixed gear lens
data$ISS_FishLenComps[,,,2] <- 40 # trawl gear lens
data$ISS_FishAgeComps[,,,1] <- 50 # fishery iss ages
data$ISS_SrvAgeComps[,,,] <- 60 # survey iss ages

# do bias ramp, helps stabalize early recruitment and helps w/ selex
data$do_rec_bias_ramp <- 1 # do bias ramp (slot 0 == don't do bias ramp, 1 == do bias ramp)
data$bias_year <- c(length(1960:1979), length(1960:1989), (length(1960:2019) - 5), length(1960:2021) - 2)

# Exploration 2
# Next, try to see what happens when we incorporate tagging data
# Here, we are also using a vague prior to push movement away from the edges as well
# as a penalty to discourage extreme values for tag reporting rates
# Movement is time and age-invariant, and tag reporting rates are estimated as time-invariant for now
# Once we incorporated tagging data, it did pretty well, converged with a gradient of ~1e-9 (pretty sweet!)

# Exploration 3
# Next, we can start increasing complexity.
# For the next set of investigations, we will allow for age-varying movement.
# Allowing for age-varying movement changes the recruitment dynamics a bunch.
# They generally make more sense ... movement in the intermediate ages are a bit wonky ... 
# but could make sense?
# Trawl selectivity is more reasonable now - i.e., doesn't go to a bound.
# It converges, with relatively low gradients which is great!

# Exploration 4
# We will continue with our exploration and allow for time-varying tag reporting rates here
# Basically allow for blocks to accommodate fishery changes. It converges, but movement in
# the intermediate ages are still a bit off. Could be due to difference in tagging data 
# used (Craig did a lot of filtering). Again, converges with pretty good gradients 
# Weird movement rates ... probably due to forgetting to accumulate when pooling tagging data when fitting. 

# Exploration 5
# Continuing our exploration and now allowing for over dispersion in tagging data
# via the use of a negative binomial likelihood - this could help alleviate
# some of the weirdness with the movement probabilities 
# Converged, but gradients are a tad bit on the higher side ~ 0.09

# Exploration 6 
# Following the negative binomial, we will then free up some of the selectivities, and first
# allow for sex-specific deltas in the survey fleet. 
# Generally seems to work well, with survey fleet selectivity estimated relatively well.
# Conerged, with good gradients ~ 1e-5

# Exploration 7 
# Following exploration 6, we will now allow for survey and sex-specific deltas for the survey fleet
# Probably converged, but coop survey results in knife edged selex

# Exploration 8
# Following exploration 7, but now allowing for sex-specific deltas for the fishery fleet
# Probably converged, but coop survey and male fixed-gear fishery results in knife edged selex

# Exploration 9 
# Following exploration 8, but now keeping coop survey delta shared between sexes, but domestic
# survey selex sex-specific
# Estimated reasonable survey selectivities
# Converges, with good gradients

# Exploration 10
# Following exploration 9, but now keeping coop survey delta shared between sexes, but domestic
# survey selex sex-specific. Also fixed-gear fishery share deltas w/ only 1 time block for now
# Also estimated decently reasonalbe fishery selectivities 
# Converges, but gradients are high ~ 0.1

# Exploration 11
# Following exploration 9, but allowing for multinomial likelihood instead of negative binomial
# Movement rates here seem to make a bit more sense compared to previous models
# Converged, with good gradients, although survey coop selex is knife edged ... 

# Exploration 12
# Following exploration 8, but fixing likelihood so it's accumulating instead of just once ... 
# Movement rates are now much more reasonable and in line with previous estimates
# However, survey coop selex is a bit knife edged, which we will need to resolve
# Converges, but gradients are high ~ 0.82

# Exploration 13 
# Following exploration 12, but allowing for fixed-gear time-varying a50s
# Movement rates are now much more reasonable and in line with previous estimates
# However, survey coop selex is a bit knife edged, which we will need to resolve
# Fishery selex estimates seem to generally make sense
# Converges, and gradients are good!

# Exploration 14 
# Following exploration 13, but sharing survey selectivity parameters of deltas across
# sexes. i.e., reverting back to the version of exploration 5-6 where survey delta are
# not survey-specific but are sex-specific
# This helps resolve the knife edged selectivity issue!
# Converges, with good gradients

# Exploration 15 
# Exploration 14, but Don't estimate rec devs only in the last year. 
# Converges, with good gradients

# Exploration 16 
# Exploration 15, but change mixing period to 3 and max tag liberty to 15
# Converges with good gradients

# Exploration 17 
# Exploration 16, but lock down those early recruitments by using sigmaRswitch to 1960:1977
data$sigmaR_switch <- length(1960:1975) # when to switch to late sigmaR 
# Helps lock down those early recruitments
# Converges with good gradients

# Exploration 18 
# Exploration 17, but change recruitment to be multinomial == 0
# Works, this should also be the more "correct" way to use a multinomial
# Converges with good gradients

# Exploration 19 (current parameterization)
# Exploration 18, but changing tag reporting rate to only do the 2 fishery transitions (pre and post IFQ)
# Converged, with good gradients ~ 1e-10

# Exploration 20 
# Exploration 19, but allowing for age and time-varying movement
# map_Movement_Pars[,,1:(length(data$years) / 2),1:6,] <- 1:20 # age block 1-6 for years 1 - 31 and sexes
# map_Movement_Pars[,,1:(length(data$years) / 2),7:15,] <- 21:40 # age block 7-15 for years 1 - 31 and sexes
# map_Movement_Pars[,,1:(length(data$years) / 2),16:30,] <- 41:60 # age block 16-30 for years 1 - 31and sexes
# map_Movement_Pars[,,(length(data$years) / 2 + 1):length(data$years),1:6,] <- 61:80 # age block 1-6 for years 1 - 31 and sexes
# map_Movement_Pars[,,(length(data$years) / 2 + 1):length(data$years),7:15,] <- 81:100 # age block 7-15 for years 1 - 31 and sexes
# map_Movement_Pars[,,(length(data$years) / 2 + 1):length(data$years),16:30,] <- 101:120 # age block 16-30 for years 1 - 31and sexes
# mapping$move_pars <- factor(map_Movement_Pars)
# data$map_Movement_Pars <- array(as.numeric(mapping$move_pars), dim = dim(parameters$move_pars))
# Converges and good gradients!

# Exploration 21,
# Exploration 19, but allowing for sex-varying movement
# map_Movement_Pars[,,1:length(data$years),1:6,1] <- 1:20 # age block 1-6 for all years and females
# map_Movement_Pars[,,1:length(data$years),7:15,1] <- 21:40 # age block 7-15 for all years and females
# map_Movement_Pars[,,1:length(data$years),16:30,1] <- 41:60 # age block 16-30 for all years and females
# map_Movement_Pars[,,1:length(data$years),1:6,2] <- 61:80 # age block 1-6 for all years and males
# map_Movement_Pars[,,1:length(data$years),7:15,2] <- 81:100 # age block 7-15 for all years and males
# map_Movement_Pars[,,1:length(data$years),16:30,2] <- 101:120 # age block 16-30 for all years and males
# mapping$move_pars <- factor(map_Movement_Pars)
# data$map_Movement_Pars <- array(as.numeric(mapping$move_pars), dim = dim(parameters$move_pars))
# data$move_sex_blocks <- list(1,2)
# Doesn't converge, and gradients are a bit high ~ 0.01 (negbin parameter can't be estiamted)

# Run Model ---------------------------------------------------------------
# make AD model function
data$UseTagging = 1 # Using tagging data
data$use_fixed_movement = 0
data$Use_Movement_Prior = 1
data$Movement_prior[] = 1.2 # prior to push away from edge
data$Use_TagRep_Prior = 1 # prior to push away from edge

sabie_rtmb_model <- RTMB::MakeADFun(cmb(sabie_RTMB, data), parameters = parameters, map = mapping)
sabie_rtmb_model$rep <- sabie_rtmb_model$report(sabie_rtmb_model$env$last.par.best) # Get report

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

beepr::beep(3)

sabie_rtmb_model$sd_rep$par.fixed[order(sabie_rtmb_model$sd_rep$gradient.fixed, decreasing = T)[1:10]]
rowSums(sabie_rtmb_model$rep$Movement[,,1,1,2])

sabie_rtmb_model$rep$Movement[,,1,1,1]
sabie_rtmb_model$rep$Movement[,,1,1,2]

sabie_rtmb_model$rep$Movement[,,1,13,1]
sabie_rtmb_model$rep$Movement[,,1,13,2]

sabie_rtmb_model$rep$Movement[,,1,30,1]
sabie_rtmb_model$rep$Movement[,,1,30,2]

sabie_rtmb_model$rep$Tag_Reporting
sabie_rtmb_model$rep$R0
sabie_rtmb_model$rep$srv_q
sabie_rtmb_model$optim$par[names(sabie_rtmb_model$optim$par) == 'ln_tag_theta']
sum(sabie_rtmb_model$rep$R0)

sqrt(diag(sabie_rtmb_model$sd_rep$cov.fixed))[1092]

# Plots -------------------------------------------------------------------

ssb_vals <- reshape2::melt(matrix(sabie_rtmb_model$sd_rep$value[names(sabie_rtmb_model$sd_rep$value) == 'SSB'], 5, 62))
sd_vals <- reshape2::melt(matrix(sabie_rtmb_model$sd_rep$sd[names(sabie_rtmb_model$sd_rep$value) == 'SSB'], 5, 62)) %>% 
  rename(sd = value)
ssb_vals <- ssb_vals %>% left_join(sd_vals, by = c("Var1", "Var2")) %>% 
  mutate(lwr = value - 2*sd,
         upr = value + 2*sd,
         value = value)

ssb_vals %>% 
  ggplot(aes(x = Var2, y = value, color = factor(Var1), fill = factor(Var1), ymin = lwr, ymax = upr)) +
  geom_line(lwd = 1) +
  geom_ribbon(alpha = 0.3, color = NA) +
  facet_wrap(~Var1, nrow = 1)

ssb_vals %>% 
  ggplot(aes(x = Var2, y = value, color = factor(Var1), fill = factor(Var1), ymin = lwr, ymax = upr)) +
  geom_line(lwd = 1) +
  geom_ribbon(alpha = 0.3, color = NA) +
  ylim(0,NA)

reshape2::melt(sabie_rtmb_model$rep$SSB) %>% 
  ggplot(aes(x = Var2, y = value, color = factor(Var1))) +
  geom_line(lwd = 1) +
  facet_wrap(~Var1, nrow = 1)
  
reshape2::melt(sabie_rtmb_model$rep$SSB) %>% 
  group_by(Var2) %>% 
  summarize(value = sum(value)) %>% 
  ggplot(aes(x = Var2, y = value)) +
  geom_line(lwd = 1) +
  ylim(0,NA)

reshape2::melt(sabie_rtmb_model$rep$Total_Biom) %>% 
  ggplot(aes(x = Var2, y = value, color = factor(Var1))) +
  geom_line(lwd = 1) +
  facet_wrap(~Var1, nrow = 1)

reshape2::melt(sabie_rtmb_model$rep$Rec) %>% 
  group_by(Var2) %>% 
  summarize(value = sum(value)) %>% 
  ggplot(aes(x = Var2, y = value)) +
  geom_line(lwd = 1)

reshape2::melt(sabie_rtmb_model$rep$Rec) %>% 
  ggplot(aes(x = Var2, y = value, color = factor(Var1))) +
  geom_line(lwd = 1) +
  facet_wrap(~Var1, nrow = 1)

reshape2::melt(sabie_rtmb_model$rep$Rec) %>% 
  ggplot(aes(x = Var2, y = value, color = factor(Var1))) +
  geom_line(lwd = 1) 
  # facet_wrap(~Var1, nrow = 1)

reshape2::melt(sabie_rtmb_model$rep$Fmort) %>% 
  ggplot(aes(x = Var2, y = value, color = factor(Var1))) +
  geom_line(lwd = 1) +
  facet_grid(Var3~Var1, scales = "free_y")

reshape2::melt(sabie_rtmb_model$rep$Fmort) %>% 
  ggplot(aes(x = Var2, y = value, color = factor(Var3))) +
  geom_line(lwd = 1)+
  facet_grid(~Var1, scales = "free_y")

reshape2::melt(sabie_rtmb_model$rep$fish_sel) %>% 
  rename(R = Var1, Y = Var2, A = Var3, S = Var4, Fl = Var5) %>% 
  filter(Y == 1) %>% 
  ggplot(aes(x = A, y = value, color = factor(R))) +
  geom_line() +
  facet_grid(S~Fl, scales = "free_y")

reshape2::melt(sabie_rtmb_model$rep$srv_sel) %>% 
  rename(R = Var1, Y = Var2, A = Var3, S = Var4, Fl = Var5) %>% 
  filter(Y == 50) %>% 
  ggplot(aes(x = A, y = value, color = factor(R))) +
  geom_line() +
  facet_grid(S~Fl, scales = "free_y")

# Catch
plot(data$ObsCatch[1,,1])
lines(sabie_rtmb_model$rep$PredCatch[1,,1])

plot(data$ObsCatch[2,,1])
lines(sabie_rtmb_model$rep$PredCatch[2,,1])

plot(data$ObsCatch[3,,1])
lines(sabie_rtmb_model$rep$PredCatch[3,,1])

plot(data$ObsCatch[4,,1])
lines(sabie_rtmb_model$rep$PredCatch[4,,1])

plot(data$ObsCatch[5,,1])
lines(sabie_rtmb_model$rep$PredCatch[5,,1])

# Trawl catch
plot(data$ObsCatch[1,,2])
lines(sabie_rtmb_model$rep$PredCatch[1,,2])

plot(data$ObsCatch[2,,2])
lines(sabie_rtmb_model$rep$PredCatch[2,,2])

plot(data$ObsCatch[3,,2])
lines(sabie_rtmb_model$rep$PredCatch[3,,2])

plot(data$ObsCatch[4,,2])
lines(sabie_rtmb_model$rep$PredCatch[4,,2])

plot(data$ObsCatch[5,,2])
lines(sabie_rtmb_model$rep$PredCatch[5,,2])


# Coop Survey
plot(data$ObsSrvIdx[1,data$UseSrvIdx[1,,1] == 1,1])
lines(sabie_rtmb_model$rep$PredSrvIdx[1,data$UseSrvIdx[1,,1] == 1,1])

plot(data$ObsSrvIdx[2,data$UseSrvIdx[2,,1] == 1,1])
lines(sabie_rtmb_model$rep$PredSrvIdx[2,data$UseSrvIdx[2,,1] == 1,1])

plot(data$ObsSrvIdx[3,data$UseSrvIdx[3,,1] == 1,1])
lines(sabie_rtmb_model$rep$PredSrvIdx[3,data$UseSrvIdx[3,,1] == 1,1])

plot(data$ObsSrvIdx[4,data$UseSrvIdx[4,,1] == 1,1])
lines(sabie_rtmb_model$rep$PredSrvIdx[4,data$UseSrvIdx[4,,1] == 1,1])

plot(data$ObsSrvIdx[5,data$UseSrvIdx[5,,1] == 1,1])
lines(sabie_rtmb_model$rep$PredSrvIdx[5,data$UseSrvIdx[5,,1] == 1,1])

# domestic survey
plot(data$ObsSrvIdx[1,data$UseSrvIdx[1,,2] == 1,2])
lines(sabie_rtmb_model$rep$PredSrvIdx[1,data$UseSrvIdx[1,,2] == 1,2])

plot(data$ObsSrvIdx[2,data$UseSrvIdx[2,,2] == 1,2])
lines(sabie_rtmb_model$rep$PredSrvIdx[2,data$UseSrvIdx[2,,2] == 1,2])

plot(data$ObsSrvIdx[3,data$UseSrvIdx[3,,2] == 1,2])
lines(sabie_rtmb_model$rep$PredSrvIdx[3,data$UseSrvIdx[3,,2] == 1,2])

plot(data$ObsSrvIdx[4,data$UseSrvIdx[4,,2] == 1,2])
lines(sabie_rtmb_model$rep$PredSrvIdx[4,data$UseSrvIdx[4,,2] == 1,2])

plot(data$ObsSrvIdx[5,data$UseSrvIdx[5,,2] == 1,2])
lines(sabie_rtmb_model$rep$PredSrvIdx[5,data$UseSrvIdx[5,,2] == 1,2])

plot(as.vector(data$ObsFishAgeComps[2,56,,,1]))
lines(as.vector(sabie_rtmb_model$rep$CAA[2,56,,,1] / sum(sabie_rtmb_model$rep$CAA[2,56,,,1])))

plot(rowSums(sabie_rtmb_model$rep$NAA[1,,,1] / rowSums(sabie_rtmb_model$rep$NAA[1,,,1]) * matrix(c(2:31), nrow = 63, ncol = 30, byrow = T)), type = 'l')
plot(rowSums(sabie_rtmb_model$rep$NAA[2,,,1] / rowSums(sabie_rtmb_model$rep$NAA[2,,,1]) * matrix(c(2:31), nrow = 63, ncol = 30, byrow = T)), type = 'l')
lines(rowSums(sabie_rtmb_model$rep$NAA[3,,,1] / rowSums(sabie_rtmb_model$rep$NAA[3,,,1]) * matrix(c(2:31), nrow = 63, ncol = 30, byrow = T)), type = 'l')
lines(rowSums(sabie_rtmb_model$rep$NAA[4,,,1] / rowSums(sabie_rtmb_model$rep$NAA[4,,,1]) * matrix(c(2:31), nrow = 63, ncol = 30, byrow = T)), type = 'l')
lines(rowSums(sabie_rtmb_model$rep$NAA[5,,,1] / rowSums(sabie_rtmb_model$rep$NAA[5,,,1]) * matrix(c(2:31), nrow = 63, ncol = 30, byrow = T)), type = 'l')
