# Purpose: To run TMB model
# Creator: Matthew LH. Cheng
# Date Created: 4/19/24


# Set up ------------------------------------------------------------------

library(here)
library(R2admb)
library(tidyverse)
library(TMB)

source(here("R", "Francis_Reweight.R"))

# Read in data
tem_dat <- dget(here('data', 'tem.rdat'))
ageing_dat <- dget(here('data', 'test.rdat'))
tem_admb_dat <- readLines(here("data", "tem_2023_na_wh.dat"))
tem_rep <- readLines(here('data', 'sable.rep'))
tem_par <- read_pars(here('data', 'tem'))

# Get population dynamics -------------------------------------------------
ages <- 1:30 # ages
years <- rownames(tem_dat$t.series) # years
ssb <- tem_dat$t.series$spbiom # spawning biomass
total_biom <- tem_dat$t.series$totbiom # total biomass
waa_f <- tem_dat$growthmat$wt.f.block1 # female WAA
waa_m <- tem_dat$growthmat$wt.m.block1 # male WAA
mat <- tem_dat$growthmat$mage.block2 # maturity
m_f <- exp(tem_par$coefficients[names(tem_par$coefficients) == "logm"])
m_m <- m_f-0.00819813327864 # male M

### Recruitment -------------------------------------------------------------
rec <- tem_dat$t.series$Recr # total recruitment
init_age <- matrix(cbind(tem_dat$natage.female[1,-1] , tem_dat$natage.male[1,-1] ), ncol = 2) # initial ages
mean_rec <- tem_par$coefficients[names(tem_par$coefficients) == "log_mean_rec"] # mean recruitment
init_rec_devs <- rev(tem_par$coefficients[str_detect(names(tem_par$coefficients), "rec_dev") ][1:28]) # initial devs
rec_devs <- tem_par$coefficients[str_detect(names(tem_par$coefficients), "rec_dev")][c(29:91)] # rec devs

### Fishing Mortality -------------------------------------------------------
total_f <- tem_dat$t.series$fmort # total F
mean_ll_fish <- tem_par$coefficients[names(tem_par$coefficients) == "log_avg_F_fish1"] # Average F LL
devs_ll_fish <- tem_par$coefficients[str_detect(names(tem_par$coefficients), "log_F_devs_fish1")] # Devs LL
mean_ll_trwl <- tem_par$coefficients[names(tem_par$coefficients) == "log_avg_F_fish3"] # Average F Trawl
devs_ll_trwl <- tem_par$coefficients[str_detect(names(tem_par$coefficients), "log_F_devs_fish3")] # Devs Trawl

# Prepare Data ------------------------------------------------------------
data <- list() # make data list

# Set up dimensions
data$ages <- ages # ages
data$lens <- seq(41, 99, 2) # lengths
data$yrs <- as.numeric(years) # years
data$n_sexes <- 2 # number of sexes
data$n_fish_fleets <- 2 # number of fishery fleets (fixed gear, trawl gear)
data$n_srv_fleets <- 3 # number of survey fleets (domestic ll survey, domestic trawl survey, coop jp ll survey)
data$init_F_prop <- 0.1 # initial F proportion
data$do_rec_bias_ramp <- 1 # do bias ramp
data$bias_year <- c(length(1960:1979), length(1960:1989), (length(1960:2023) - 6), length(1960:2022) - 1)
data$sigmaR_switch <- as.integer(length(1960:1975)) - 1 # when to switch sigmaR from 0.4 to 1.04
data$likelihood_type <- 0 # ADMB likelihood
data$share_sel <- 0 # share index selectivity across sexes
data$normalize_fish_sel <- c(0,1) # normalize selectivity to max at 1s
data$normalize_srv_sel <- c(0,1,0) # normalize survey selectivity to max at 1s (0 = dont normalize, 1 = normalize)

# Recruitment stuff
data$sexratio <- as.vector(c(0.5, 0.5)) # recruitment sex ratio

# Biologicals
data$Use_M_prior <- 1 # use natural mortality prior
data$M_prior <- c(0.1, 0.1) # Mean and CV for M prior
data$WAA <- array(NA, dim = c(length(data$yrs), length(data$ages), data$n_sexes)) # weight-at-age array
data$WAA[,,1] <- matrix(waa_f, nrow = length(1:64), ncol = length(1:30), byrow = TRUE) # female weight-at-age
data$WAA[,,2] <- matrix(waa_m, nrow = length(1:64), ncol = length(1:30), byrow = TRUE) # male weight-at-age
data$MatAA <- array(0,  dim = c(length(data$yrs), length(data$ages), data$n_sexes)) # maturity at age
data$MatAA[1:64,,1] <- aperm(replicate(length(1:64), mat), c(2,1))
data$MatAA[1:64,,2] <- aperm(replicate(length(1:64), mat), c(2,1))
data$AgeingError <- as.matrix(ageing_dat$age_error) # ageing error matrix
data$SizeAgeTrans <- array(0, dim = c(length(data$yrs), length(data$lens), length(data$ages), data$n_sexes)) # size age transition matrix
rownames(data$SizeAgeTrans) <- data$yrs
# Size age transition, first time block
data$SizeAgeTrans[1:35,,,1] <- aperm(replicate(length(1:35), tem_dat$sizeage.f.block1), perm = c(3,2,1)) # female
data$SizeAgeTrans[1:35,,,2] <- aperm(replicate(length(1:35), tem_dat$sizeage.m.block1), perm = c(3,2,1)) # male

# Size age transition, second time block
data$SizeAgeTrans[36:64,,,1] <- aperm(replicate(length(36:64), tem_dat$sizeage.f.block2), perm = c(3,2,1)) # female
data$SizeAgeTrans[36:64,,,2] <- aperm(replicate(length(36:64), tem_dat$sizeage.m.block2), perm = c(3,2,1)) # male

## Observations  -----------------------------------------------------------
### Fishery Observations ----------------------------------------------------
# Catches
data$ObsCatch <- matrix(NA, nrow = length(1:64), ncol = data$n_fish_fleets)
rownames(data$ObsCatch) <- data$yrs # define row years
data$ObsCatch[,1] <- as.numeric(strsplit(tem_admb_dat[46], split = " ")[[1]]) # fixed gear catches
data$ObsCatch[,2] <- as.numeric(strsplit(tem_admb_dat[48], split = " ")[[1]]) # trawl gear catches
data$ObsCatch[data$ObsCatch == 0] <- NA # set 0 catches to NA so we aren't fitting
data$UseCatch <- matrix(1, nrow = length(1:64), ncol = data$n_fish_fleets)
data$UseCatch[1:3,2] <- 0 # don't use first 3 years of trawl catches (0)
data$Catch_Constant <- c(0.01, 0.8) # catch constants
data$Catch_sd <- c(0.01,0.01)

# Fishery Indices
data$ObsFishIdx <- matrix(NA, nrow = length(1:64), ncol = data$n_fish_fleets)
data$ObsFishIdx_SE <- matrix(NA, nrow = length(1:64), ncol = data$n_fish_fleets)
data$UseFishIdx <- matrix(1, nrow = length(1:64), ncol = data$n_fish_fleets)
rownames(data$ObsFishIdx) <- data$yrs # define row years
rownames(data$ObsFishIdx_SE) <- data$yrs # define row years

# Fixed gear fishery CPUE
# Domestic LL + Japanese Fishery (pre 1995)
data$ObsFishIdx[rownames(data$ObsFishIdx) %in% rownames(tem_dat$obssrv6),1] <- as.numeric(strsplit(tem_admb_dat[144], split = " ")[[1]])
data$ObsFishIdx_SE[rownames(data$ObsFishIdx_SE) %in% rownames(tem_dat$obssrv6),1] <- as.numeric(strsplit(tem_admb_dat[146], split = " ")[[1]])

# Domestic LL Fishery (after 1995)
data$ObsFishIdx[rownames(data$ObsFishIdx) %in% rownames(tem_dat$obssrv5),1] <- as.numeric(strsplit(tem_admb_dat[127], split = " ")[[1]])
data$ObsFishIdx_SE[rownames(data$ObsFishIdx_SE) %in% rownames(tem_dat$obssrv5),1] <- as.numeric(strsplit(tem_admb_dat[129], split = " ")[[1]])
data$UseFishIdx[is.na(data$ObsFishIdx)] <- 0 # don't fit if missing data

# Fishery Age Comps
# Note that NA is in trawl fishery slot so it doesn't fit
data$ObsFishAgeComps <- array(NA, dim = c(length(data$yrs), length(data$ages), data$n_sexes, data$n_fish_fleets))
rownames(data$ObsFishAgeComps) <- data$yrs # define row years
data$ObsFishAgeComps[rownames(data$ObsFishAgeComps) %in% rownames(tem_dat$oac.fish1),,1,1] <- tem_dat$oac.fish1 # Observed fixed gear fishery age comps (aggregated)
data$UseFishAgeComps <- matrix(0, nrow = length(data$yrs), ncol = data$n_fish_fleets)
rownames(data$UseFishAgeComps) <- data$yrs # define row years
data$UseFishAgeComps[rownames(data$UseFishAgeComps) %in% rownames(tem_dat$oac.fish1),1] <- 1 # only fit if have age comp data
data$AggFishAgeComps <- matrix(0, nrow = length(data$yrs), ncol = data$n_fish_fleets) # aggregating all age comps

# Data weighting for fishery age compositions
data$ISS_FishAgeComps <- array(0, dim = c(length(data$yrs), data$n_sexes, data$n_fish_fleets))
rownames(data$ISS_FishAgeComps) <- data$yrs # define row years
data$ISS_FishAgeComps[rownames(data$ISS_FishAgeComps) %in% rownames(tem_dat$oac.fish1),1,1] <- 20 # Assuming constant ISS of 20 for fixed gear fishery ages (aggregated)
data$Wt_FishAgeComps <- array(NA, dim = c(data$n_sexes, data$n_fish_fleets)) # weights for fishery age comps
data$Wt_FishAgeComps[1,1] <- 0.797868466479416 # Weight for fixed gear age comps

# Fishery Length Comps
data$ObsFishLenComps <- array(NA, dim = c(length(data$yrs), length(data$lens), data$n_sexes, data$n_fish_fleets))
rownames(data$ObsFishLenComps) <- data$yrs # define row years

# observed fixed gear fishery length comps
data$ObsFishLenComps[rownames(data$ObsFishLenComps) %in% rownames(tem_dat$olc.fish1.f),,1,1] <- tem_dat$olc.fish1.f # females
data$ObsFishLenComps[rownames(data$ObsFishLenComps) %in% rownames(tem_dat$olc.fish1.m),,2,1] <- tem_dat$olc.fish1.m # males

# observed trawl gear fishery length comps
data$ObsFishLenComps[rownames(data$ObsFishLenComps) %in% rownames(tem_dat$olc.fish3.f),,1,2] <- tem_dat$olc.fish3.f # females
data$ObsFishLenComps[rownames(data$ObsFishLenComps) %in% rownames(tem_dat$olc.fish3.m),,2,2] <- tem_dat$olc.fish3.m # males

# Use fishery length comp controls
data$UseFishLenComps <- matrix(0, nrow = length(data$yrs), ncol = data$n_fish_fleets)
rownames(data$UseFishLenComps) <- data$yrs # define row years
data$UseFishLenComps[rownames(data$UseFishLenComps) %in% rownames(tem_dat$olc.fish1.f),1] <- 1 # only fit if have len comp data
data$UseFishLenComps[rownames(data$UseFishLenComps) %in% rownames(tem_dat$olc.fish3.m),2] <- 1 # only fit if have len comp data

# Data weighting for fishery length compositions
data$ISS_FishLenComps <- array(0, dim = c(length(data$yrs), data$n_sexes, data$n_fish_fleets))
rownames(data$ISS_FishLenComps) <- data$yrs # define row years
data$ISS_FishLenComps[rownames(data$ISS_FishLenComps) %in% rownames(tem_dat$olc.fish1.f),1,1] <- 20 # Assuming constant ISS of 20 for fixed gear fishery females
data$ISS_FishLenComps[rownames(data$ISS_FishLenComps) %in% rownames(tem_dat$olc.fish1.m),2,1] <- 20 # Assuming constant ISS of 20 for fixed gear fishery males
data$ISS_FishLenComps[rownames(data$ISS_FishLenComps) %in% rownames(tem_dat$olc.fish3.f),1,2] <- 20 # Assuming constant ISS of 20 for trawl gear fishery females
data$ISS_FishLenComps[rownames(data$ISS_FishLenComps) %in% rownames(tem_dat$olc.fish3.m),2,2] <- 20 # Assuming constant ISS of 20 for trawl gear fishery males
data$Wt_FishLenComps <- array(0, dim = c(data$n_sexes, data$n_fish_fleets)) # weights for fishery age comps
data$Wt_FishLenComps[1,1] <- 4.94489032547033 # Weight for fixed gear len comps females
data$Wt_FishLenComps[2,1] <- 5.21629903889062 # Weight for fixed gear len comps males
data$Wt_FishLenComps[1,2] <- 0.35007890900583 # Weight for trawl gear len comps females
data$Wt_FishLenComps[2,2] <- 0.255150139457704 # Weight for trawl gear len comps males

### Survey Observations -----------------------------------------------------

# Survey Indices
data$ObsSrvIdx <- matrix(NA, nrow = length(1:64), ncol = data$n_srv_fleets)
data$ObsSrvIdx_SE <- matrix(NA, nrow = length(1:64), ncol = data$n_srv_fleets)
data$UseSrvIdx <- matrix(1, nrow = length(1:64), ncol = data$n_srv_fleets)
rownames(data$ObsSrvIdx) <- data$yrs # define row years
rownames(data$ObsSrvIdx_SE) <- data$yrs # define row years

# Domestic LL Survey post 1995
data$ObsSrvIdx[rownames(data$ObsSrvIdx) %in% rownames(tem_dat$obssrv3),1] <- as.numeric(strsplit(tem_admb_dat[93], split = " ")[[1]])
data$ObsSrvIdx_SE[rownames(data$ObsSrvIdx_SE) %in% rownames(tem_dat$obssrv3),1] <- as.numeric(strsplit(tem_admb_dat[95], split = " ")[[1]])

# Domestic Trawl Survey
data$ObsSrvIdx[rownames(data$ObsSrvIdx) %in% rownames(tem_dat$obssrv7),2] <- as.numeric(strsplit(tem_admb_dat[161], split = " ")[[1]])
data$ObsSrvIdx_SE[rownames(data$ObsFishIdx_SE) %in% rownames(tem_dat$obssrv7),2] <- as.numeric(strsplit(tem_admb_dat[163], split = " ")[[1]])

# Coop LL Survey pre 1995
data$ObsSrvIdx[rownames(data$ObsSrvIdx) %in% rownames(tem_dat$obssrv4),3] <- as.numeric(strsplit(tem_admb_dat[110], split = " ")[[1]])
data$ObsSrvIdx_SE[rownames(data$ObsSrvIdx_SE) %in% rownames(tem_dat$obssrv4),3] <- as.numeric(strsplit(tem_admb_dat[112], split = " ")[[1]])
data$UseSrvIdx[is.na(data$ObsSrvIdx)] <- 0 # don't fit if missing data

# Survey Age Comps
# Note that NA is in trawl survey slot so it doesn't fit
data$ObsSrvAgeComps <- array(NA, dim = c(length(data$yrs), length(data$ages), data$n_sexes, data$n_srv_fleets))
rownames(data$ObsSrvAgeComps) <- data$yrs # define row years
data$ObsSrvAgeComps[rownames(data$ObsSrvAgeComps) %in% rownames(tem_dat$oac.srv1),,1,1] <- tem_dat$oac.srv1 # Observed domestic ll survey age comps (aggregated)
data$ObsSrvAgeComps[rownames(data$ObsSrvAgeComps) %in% rownames(tem_dat$oac.srv2),,1,3] <- tem_dat$oac.srv2 # Observed coop jp ll survey age comps (aggregated)

data$UseSrvAgeComps <- matrix(0, nrow = length(data$yrs), ncol = data$n_srv_fleets)
rownames(data$UseSrvAgeComps) <- data$yrs # define row years
data$UseSrvAgeComps[rownames(data$UseSrvAgeComps) %in% rownames(tem_dat$oac.srv1),1] <- 1 # only fit if have age comp data
data$UseSrvAgeComps[rownames(data$UseSrvAgeComps) %in% rownames(tem_dat$oac.srv2),3] <- 1 # only fit if have age comp data
data$AggSrvAgeComps <- matrix(0, nrow = length(data$yrs), ncol = data$n_srv_fleets) # aggregating all age comps

# Data weighting for fishery age compositions
data$ISS_SrvAgeComps <- array(0, dim = c(length(data$yrs), data$n_sexes, data$n_srv_fleets))
rownames(data$ISS_SrvAgeComps) <- data$yrs # define row years
data$ISS_SrvAgeComps[rownames(data$ISS_SrvAgeComps) %in% rownames(tem_dat$oac.srv1),1,1] <- 20 # Assuming constant ISS of 20 for domestic ll survey ages (aggregated)
data$ISS_SrvAgeComps[rownames(data$ISS_SrvAgeComps) %in% rownames(tem_dat$oac.srv2),1,3] <- 20 # Assuming constant ISS of 20 for coop jp survey ages (aggregated)

data$Wt_SrvAgeComps <- array(NA, dim = c(data$n_sexes, data$n_srv_fleets)) # weights for survey age comps
data$Wt_SrvAgeComps[1,1] <- 3.72382880611179 # Weight for domestic survey ll gear age comps
data$Wt_SrvAgeComps[1,3] <- 1.27151115308662 # Weight for coop jp survey ll gear age comps

# Survey Length Comps
data$ObsSrvLenComps <- array(NA, dim = c(length(data$yrs), length(data$lens), data$n_sexes, data$n_srv_fleets))
rownames(data$ObsSrvLenComps) <- data$yrs # define row years

# observed domestic survey ll length comps
data$ObsSrvLenComps[rownames(data$ObsSrvLenComps) %in% rownames(tem_dat$olc.srv1.f),,1,1] <- tem_dat$olc.srv1.f # females
data$ObsSrvLenComps[rownames(data$ObsSrvLenComps) %in% rownames(tem_dat$olc.srv1.m),,2,1] <- tem_dat$olc.srv1.m # males

# observed domestic trawl survey length comps
data$ObsSrvLenComps[rownames(data$ObsSrvLenComps) %in% rownames(tem_dat$olc.srv7.f),,1,2] <- tem_dat$olc.srv7.f # females
data$ObsSrvLenComps[rownames(data$ObsSrvLenComps) %in% rownames(tem_dat$olc.srv7.m),,2,2] <- tem_dat$olc.srv7.m # males

# observed coop jp ll survey length comps
srv_trawl_iss <- as.numeric(strsplit(tem_admb_dat[580], split = " ")[[1]]) # get ISS from trawl survey
data$ObsSrvLenComps[rownames(data$ObsSrvLenComps) %in% rownames(tem_dat$olc.srv2.f),,1,3] <- tem_dat$olc.srv2.f # females
data$ObsSrvLenComps[rownames(data$ObsSrvLenComps) %in% rownames(tem_dat$olc.srv2.m),,2,3] <- tem_dat$olc.srv2.m # males

# Use survey length comp controls
data$UseSrvLenComps <- matrix(0, nrow = length(data$yrs), ncol = data$n_srv_fleets)
rownames(data$UseSrvLenComps) <- data$yrs # define row years
data$UseSrvLenComps[rownames(data$UseSrvLenComps) %in% rownames(tem_dat$olc.srv1.f),1] <- 1 # only fit if have len comp data
data$UseSrvLenComps[rownames(data$UseSrvLenComps) %in% rownames(tem_dat$olc.srv7.f),2] <- 1 # only fit if have len comp data
data$UseSrvLenComps[rownames(data$UseSrvLenComps) %in% rownames(tem_dat$olc.srv2.f),3] <- 1 # only fit if have len comp data

# Data weighting for survey length compositions
data$ISS_SrvLenComps <- array(0, dim = c(length(data$yrs), data$n_sexes, data$n_srv_fleets))
rownames(data$ISS_SrvLenComps) <- data$yrs # define row years
data$ISS_SrvLenComps[rownames(data$ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv1.f),1,1] <- 20 # Assuming constant ISS of 20 for domestic ll survey females
data$ISS_SrvLenComps[rownames(data$ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv1.m),2,1] <- 20 # Assuming constant ISS of 20 for domestic ll survey males
data$ISS_SrvLenComps[rownames(data$ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv7.f),1,2] <- srv_trawl_iss # Assuming constant ISS of 20 for domestic trawl survey females
data$ISS_SrvLenComps[rownames(data$ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv7.m),2,2] <- srv_trawl_iss # Assuming constant ISS of 20 for domestic trawl survey males
data$ISS_SrvLenComps[rownames(data$ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv2.f),1,3] <- 20 # Assuming constant ISS of 20 for coop jp ll survey females
data$ISS_SrvLenComps[rownames(data$ISS_SrvLenComps) %in% rownames(tem_dat$olc.srv2.m),2,3] <- 20 # Assuming constant ISS of 20 for coop jp ll survey males

# data weights
data$Wt_SrvLenComps <- array(0, dim = c(data$n_sexes, data$n_srv_fleets)) # weights for survey age comps
data$Wt_SrvLenComps[1,1] <- 1.49981566417443 # Weight for domestic ll survey len comps females
data$Wt_SrvLenComps[2,1] <- 1.11519421735366 # Weight for domestic ll survey len comps males
data$Wt_SrvLenComps[1,2] <- 0.672532119442201 # Weight for domestic trawl survey len comps females
data$Wt_SrvLenComps[2,2] <- 0.450498802823617 # Weight for domestic trawl survey len comps males
data$Wt_SrvLenComps[1,3] <- 1.26841555181968 # Weight for coop jp ll survey len comps females
data$Wt_SrvLenComps[2,3] <- 0.902165573415091 # Weight for coop jp ll survey len comps males

### Fishery Stuff -----------------------------------------------------
# Selectivity
# Time Block Specification
data$fish_sel_blocks <- matrix(NA, nrow = length(1:64), ncol = data$n_fish_fleets)
data$fish_sel_blocks[1:35,1] <- 0 # block one fishery ll selex
data$fish_sel_blocks[36:56,1] <- 1 # block two fishery ll selex
data$fish_sel_blocks[57:64,1] <- 2 # block three fishery ll selex
data$fish_sel_blocks[,2] <- 0 # block two fishery ll selex

# Selectivity Model
data$fish_sel_model <- matrix(nrow = length(data$yrs), ncol = data$n_fish_fleets)
data$fish_sel_model[,1] <- 0 # Logistic selectivity for fixed gear
data$fish_sel_model[,2] <- 1 # Gamma selectivity for trawl fishery

# Catchability
# Time Block Specification
data$fish_q_blocks <- data$fish_sel_blocks # catchability blocks are same as selectivity blocks
data$fish_idx_type <- as.vector(c(1,1)) # fishery index type (0 == abundance, 1 == biomass)

### Survey Selectivity ------------------------------------------------------
# Time Block Specification
data$srv_sel_blocks <- matrix(NA, nrow = length(1:64), ncol = data$n_srv_fleets)
data$srv_sel_blocks[1:56,1] <- 0 # block one survey ll selex
data$srv_sel_blocks[57:64,1] <- 1 # block two survey ll selex
data$srv_sel_blocks[,2] <- 0 # block one trawl survey selex
data$srv_sel_blocks[,3] <- 0 # block one coop jp ll survey selex

# Selectivity Model
data$srv_sel_model <- matrix(nrow = length(data$yrs), ncol = data$n_srv_fleets)
data$srv_sel_model[,1] <- 0 # Logistic selectivity for longline survey
data$srv_sel_model[,2] <- 2 # power function selectivity for trawl survey
data$srv_sel_model[,3] <- 0 # logistic selectivity for jp ll survey

### Survey Catchability ----------------------------------------------------
# Time Block Specification
data$srv_q_blocks <- matrix(0, nrow = length(1:64), ncol = data$n_srv_fleets) # catchability blocks are same as selectivity blocks
data$srv_idx_type <- as.vector(c(0,1,0)) # survey index type (0 == abundance, 1 == biomass)

# Prepare Parameters ------------------------------------------------------
parameters <- list()
parameters$dummy <- 1

### Fishery Stuff ---------------------------------------------------------

# Fishing Mortality
parameters$ln_F_mean <- as.vector(c(mean_ll_fish, mean_ll_trwl)) # mean fishing mortality
parameters$ln_F_devs <- matrix(0, nrow = length(data$yrs), ncol = data$n_fish_fleets) # fishing mortality deviations from mean
parameters$ln_F_devs[,1] <- devs_ll_fish # longline fishery Fs
parameters$ln_F_devs[,2] <- c(rep(0, 3), devs_ll_trwl ) # trwl fishery Fs

# Fixed Gear Fishery three time blocks
max_fish_blks <- 3 # maximum number of fishery blocks for any fleet
max_fish_pars <- 2 # maximum number of fishery fixed parameters for any fleet

parameters$ln_fish_fixed_sel_pars <- array(0, dim = c(max_fish_pars, max_fish_blks, data$n_sexes, data$n_fish_fleets))
parameters$ln_fish_fixed_sel_pars[,1,1,1] <- c(1.4257e+000, -6.9100e-001) # Female longline fishery first time block parameters (a50 and then delta)
parameters$ln_fish_fixed_sel_pars[,2,1,1] <- c(1.2222e+000, 5.3936e-001) # Female longline fishery second time block parameters (a50 and then delta)
parameters$ln_fish_fixed_sel_pars[,3,1,1] <- c(6.5829e-001, 8.1499e-001) # Female longline fishery third time block parameters (a50 and then delta)

parameters$ln_fish_fixed_sel_pars[,1,2,1] <- c(1.9380e+000, -6.9100e-001) # Male longline fishery first time block parameters (a50 and then delta)
parameters$ln_fish_fixed_sel_pars[,2,2,1] <- c(1.4919e+000, -9.7358e-002) # Male longline fishery second time block parameters (a50 and then delta)
parameters$ln_fish_fixed_sel_pars[,3,2,1] <- c(1.0631e+000, -3.7957e-001) # Male longline fishery third time block parameters (a50 and then delta)

# Trawl Fishery kept the same across blocks
parameters$ln_fish_fixed_sel_pars[,1,1,2] <- c(1.7701e+000, 2.3315e+000) # Female trawl fishery first time block parameters (amax and then power)
parameters$ln_fish_fixed_sel_pars[,2,1,2] <- c(1.7701e+000, 2.3315e+000) # Female trawl fishery first time block parameters (amax and then power)
parameters$ln_fish_fixed_sel_pars[,3,1,2] <- c(1.7701e+000, 2.3315e+000) # Female trawl fishery first time block parameters (amax and then power)

# Trawl Fishery kept the same across blocks
parameters$ln_fish_fixed_sel_pars[,1,2,2] <- c(2.1048e+000, 2.3315e+000) # Male trawl fishery first time block parameters (amax and then power)
parameters$ln_fish_fixed_sel_pars[,2,2,2] <- c(2.1048e+000, 2.3315e+000) # Male trawl fishery first time block parameters (amax and then power)
parameters$ln_fish_fixed_sel_pars[,3,2,2] <- c(2.1048e+000, 2.3315e+000) # Male trawl fishery first time block parameters (amax and then power)

# Fishery Catchability
# Fixed Gear Fishery Catchability
parameters$ln_fish_q <- array(0, dim = c(max_fish_blks, data$n_fish_fleets))
parameters$ln_fish_q[1,1] <- 1.5737e+000 # first fixed gear fishery time block
parameters$ln_fish_q[2,1] <- -6.5299e+000 # second fixed gear fishery time block
parameters$ln_fish_q[3,1] <- -7.2336e+000 # third fixed gear fishery time block

# Trawl Gear Fishery Catchability
parameters$ln_fish_q[,2] <- 0 # trawl gear fishery time block (none and no index used)

### Survey Stuff -----------------------------------------------------
# Survey Selectivity
max_srv_blks <- 2 # maximum number of survey blocks for any fleet
max_srv_pars <- 2 # maximum number of survey fixed parameters for any fleet
parameters$ln_srv_fixed_sel_pars <- array(0, dim = c(max_srv_pars, max_srv_blks, data$n_sexes, data$n_srv_fleets))

# Longline Survey - 2 time blocks
parameters$ln_srv_fixed_sel_pars[,1,1,1] <- c(1.0927e+000, 1.0046e+000) # Female longline survey first time block parameters (a50 and then delta)
parameters$ln_srv_fixed_sel_pars[,2,1,1] <- c(6.3434e-001, 1.0046e+000) # Female longline survey second time block parameters (a50 and then delta)
parameters$ln_srv_fixed_sel_pars[,1,2,1] <- c(9.6856e-001, 8.8394e-001) # Male longline survey first time block parameters (a50 and then delta)
parameters$ln_srv_fixed_sel_pars[,2,2,1] <- c(5.4886e-001, 8.8394e-001) # Male longline survey second time block parameters (a50 and then delta)

# Trawl Survey - power function, single time block
parameters$ln_srv_fixed_sel_pars[,,1,2] <- -3.9394e-001 # Female trawl survey first time block parameters (power)
parameters$ln_srv_fixed_sel_pars[,,2,2] <- -1.3181e+000 # Male trawl survey first time block parameters (power)

# Coop JP Survey (Logistic) Single time block
parameters$ln_srv_fixed_sel_pars[,,1,3] <- c(0.953479618491, 1.0046e+000) # Female trawl survey first time block parameters (a50, delta)
parameters$ln_srv_fixed_sel_pars[,,2,3] <- c(1.23386822176, 8.8394e-001) # Male trawl survey first time block parameters (a50, delta)

# Survey Catchability
max_q_srv_blks <- 1 # maximum catchability blocks

# Longline Survey Catchability
parameters$ln_srv_q <- array(0, dim = c(max_q_srv_blks, data$n_srv_fleets))
parameters$ln_srv_q[,1] <- 1.8582e+000 # first longline survey (domestic) time block

# Trawl Survey Catchability
parameters$ln_srv_q[,2] <- -1.5314e-001 # trawl survey time block

# Coop JP Survey Catchability
parameters$ln_srv_q[,3] <- 1.5266e+000

### Natural Mortality -------------------------------------------------------
parameters$ln_M <- tem_par$coefficients[names(tem_par$coefficients) == "logm"] # Female M (base)
parameters$M_offset <- -0.00819813327864 # Male M offset (accidently jittered from OG assessment)

# Recruitment -------------------------------------------------------------
parameters$ln_R0 <- tem_par$coefficients[names(tem_par$coefficients) == "log_mean_rec"] # mean recruitment
parameters$ln_InitDevs <- rev(tem_par$coefficients[str_detect(names(tem_par$coefficients), "rec_dev")][c(1:28)])
parameters$ln_RecDevs <- tem_par$coefficients[str_detect(names(tem_par$coefficients), "rec_dev")][c(29:91)]
parameters$ln_sigmaR_early <- log(0.4) # early sigma R
parameters$ln_sigmaR_late <- log(1.044523) # late sigma R

# Mapping -----------------------------------------------------------------
### Fixing all parameters for checking --------------------------------------
mapping <- list()
mapping$ln_fish_fixed_sel_pars <- factor(rep(NA, length(parameters$ln_fish_fixed_sel_pars))) # fix all fishery selectivity parameters for now
mapping$ln_srv_fixed_sel_pars <- factor(rep(NA, length(parameters$ln_srv_fixed_sel_pars))) # fix all survey selectivity parameters for now
mapping$ln_fish_q <- factor(rep(NA, length(parameters$ln_fish_q))) # fix all fishery catchability parameters for now
mapping$ln_srv_q <- factor(rep(NA, length(parameters$ln_srv_q))) # fix all fishery catchability parameters for now
mapping$ln_F_mean <- factor(rep(NA, length(parameters$ln_F_mean))) # fix all fishery mort parameters for now
mapping$ln_F_devs <- factor(rep(NA, length(parameters$ln_F_devs))) # fix all fishery mort parameters for now
mapping$ln_M <- factor(NA) # fix natural mortality
mapping$ln_R0 <- factor(NA) # fix mean recruitment
mapping$ln_InitDevs <- factor(rep(NA, length(parameters$ln_InitDevs)))
# mapping$ln_RecDevs <- factor(rep(NA, length(parameters$ln_RecDevs)))
mapping$ln_sigmaR_late <- factor(NA) # fix late sigma R
mapping$ln_sigmaR_early <- factor(NA) # fix early sigma R
mapping$M_offset <- factor(NA) # fix natural mortality offset]

### Esimating parameters with sharing ---------------------------------------
# mapping <- list()
mapping$dummy <- factor(NA)
# mapping$ln_F_devs <- factor(c(1:length(data$yrs), rep(NA, 3), 65:(128-3))) # fix all fishery mort parameters for now
# mapping$ln_fish_q <- factor(c(1,2,3,rep(NA,3))) # estimate catchabilities only for first fishery with index
# sharing delta across sexes from early domestic fishery (first time block)
# also fixing parameters so that no time block for trawl fishery
# mapping$ln_fish_fixed_sel_pars <- factor(c(1:7, 2, 8:11, rep(12:13,3), rep(c(14,13),3)))

# ll survey, share delta female (index 2) across time blocks and to the coop jp ll survey delta
# ll survey, share delta male (index 5) across time blocks and to the coop jp ll survey delta
# coop jp survey does not estimate parameters and shares deltas with longline survey
# single time block with trawl survey and only one parameter hence, only one parameter estimated across blocks (indices 7 and 8)
# mapping$ln_srv_fixed_sel_pars <- factor(c(1:3, 2, 4:6, 5,
#                                           rep(7,4), rep(8, 4),
#                                           rep(c(NA,2), 2), rep(c(NA, 5), 2)))

# TMB Model -------------------------------------------------------------------
setwd(here("src"))
compile("SabieTMB.cpp")
dyn.load(dynlib('SabieTMB'))
dyn.unload(dynlib('SabieTMB'))
dyn.load(dynlib('SabieTMB'))

# make AD model function
sabie_model <- MakeADFun(data = data, parameters = parameters,
                       map = mapping, random = NULL,
                       DLL = "SabieTMB")
# Now, optimize the function
sabie_optim <- stats::nlminb(sabie_model$par, sabie_model$fn, sabie_model$gr,
                           control = list(iter.max = 1e5, eval.max = 1e5))

# newton steps
try_improve <- tryCatch(expr =
                         for(i in 1:3) {
                           g = as.numeric(sabie_model$gr(sabie_optim$par))
                           h = optimHess(sabie_optim$par, fn = sabie_model$fn, gr = sabie_model$gr)
                           sabie_optim$par = sabie_optim$par - solve(h,g)
                           sabie_optim$objective = sabie_model$fn(sabie_optim$par)
                         }
                       , error = function(e){e}, warning = function(w){w})

sabie_model$optim <- sabie_optim # Save optimized model results
sabie_model$rep <- sabie_model$report(sabie_model$env$last.par.best) # Get report
sabie_model$sd_rep <- sdreport(sabie_model) # Get sd report


# Likelihoods -------------------------------------------------------------

sum(sabie_model$rep$Catch_nLL) * 50 # 3.75378 (admb)
sum(sabie_model$rep$Init_Rec_nLL, sabie_model$rep$Rec_nLL) * 1.5 * 0.5 # 27.9273 (admb)
sum(sabie_model$rep$Fmort_Pen) * 0.1 # 6.04329 (admb)
sabie_model$rep$M_Pen # 0.766277 (admb)

# Check consistency -------------------------------------------------------

# # Predicted Indices
fish_idx <- sabie_model$rep$PredFishIdx
rownames(fish_idx) <- data$yrs
plot(tem_dat$obssrv6$obssrv6)
lines(fish_idx[rownames(fish_idx) %in% rownames(tem_dat$obssrv6), 1], col = "red") # Domestic + JP Fishery (1960 - 1994)
lines(tem_dat$obssrv6$predsrv6)

plot((fish_idx[rownames(fish_idx) %in% rownames(tem_dat$obssrv6), 1] - tem_dat$obssrv6$predsrv6) / tem_dat$obssrv6$predsrv6 * 100)

plot(sabie_model$rep$PredFishIdx[36:63,1], col = "red", type = 'l') # Fixed Gear Fishery (1995 - terminal)
lines(tem_dat$obssrv5$predsrv5)
points((tem_dat$obssrv5$obssrv5))

plot((sabie_model$rep$PredFishIdx[36:63,1] - tem_dat$obssrv5$predsrv5) / tem_dat$obssrv5$predsrv5 * 100)

# # Survey Index
plot(tem_dat$obssrv3$obssrv3)
lines(sabie_model$rep$PredSrvIdx[31:64,1], col = "red")
lines(tem_dat$obssrv3$predsrv3)
plot((sabie_model$rep$PredSrvIdx[31:64,1] - tem_dat$obssrv3$predsrv3) / tem_dat$obssrv3$predsrv3 * 100)

#
srv <- sabie_model$rep$PredSrvIdx
rownames(srv) <- data$yrs
plot(tem_dat$obssrv4$obssrv4)
lines(srv[rownames(srv) %in% rownames(tem_dat$obssrv4),3], col = "red")
lines(tem_dat$obssrv4$predsrv4)
plot((srv[rownames(srv) %in% rownames(tem_dat$obssrv4),3]- tem_dat$obssrv4$predsrv4) / tem_dat$obssrv4$predsrv4 * 100)


50 * sum(sabie_model$rep$Catch_nLL)

y <- 20

# Check ZAA with actual admb
data$ZAA[y,,1] - sabie_model$rep$natmort[y,,1] - (sabie_model$rep$fish_sel[y,,1,1] * sabie_model$rep$Fmort[y,1] +
                                                  sabie_model$rep$fish_sel[y,,1,2]  * sabie_model$rep$Fmort[y,2])

# Check mortality with ADMB
sabie_model$rep$natmort[y,,1] - exp(parameters$ln_M)

# Check fishing mortality with ADMB (fixed gear)
sabie_model$rep$Fmort[y,1] - exp(parameters$ln_F_mean[1] + parameters$ln_F_devs[y,1])

# Check fishing mortality with ADMB (trawl gear)
sabie_model$rep$Fmort[y,2] - exp(parameters$ln_F_mean[2] + parameters$ln_F_devs[y,2])

# Check selectivity (fixed gear)
sabie_model$rep$fish_sel[1,,1,1] - tem_dat$agesel$fish1sel.f
sabie_model$rep$fish_sel[1,,2,1] - tem_dat$agesel$fish1sel.m

sabie_model$rep$fish_sel[50,,1,1] - tem_dat$agesel$fish4sel.f
sabie_model$rep$fish_sel[50,,2,1] - tem_dat$agesel$fish4sel.m

# Trawl Gear selectivity
sabie_model$rep$fish_sel[y,,1,2] - tem_dat$agesel$fish3sel.f


# Check fishing mortality AA with ADMB FAA (fixed gear)
(sabie_model$rep$fish_sel[y,,1,1] * sabie_model$rep$Fmort[y,1]) - 
(tem_dat$agesel$fish1sel.f * exp(parameters$ln_F_mean[1] + parameters$ln_F_devs[y,1]))

# Check fishing mortality AA with ADMB FAA (trawl gear)
(sabie_model$rep$fish_sel[y,,1,2]  * sabie_model$rep$Fmort[y,2]) - 
(tem_dat$agesel$fish3sel.f * exp(parameters$ln_F_mean[2] + parameters$ln_F_devs[y,2]))











   
                                
(sabie_model$rep$fish_sel[y,,1,1] * sabie_model$rep$Fmort[y,1] +
 sabie_model$rep$fish_sel[y,,1,2]  * sabie_model$rep$Fmort[y,2]) 

(tem_dat$agesel$fish1sel.f * exp(parameters$ln_F_mean[1] + parameters$ln_F_devs[y,1]) +
tem_dat$agesel$fish3sel.f * exp(parameters$ln_F_mean[2] + parameters$ln_F_devs[y,2]))

data$ZAA[y,,1] - sabie_model$rep$natmort[y,,1] - (sabie_model$rep$fish_sel[y,,1,1] * sabie_model$rep$Fmort[y,1] +
                                                  sabie_model$rep$fish_sel[y,,1,2]  * sabie_model$rep$Fmort[y,2])
sabie_model$rep$TotalFAA[y,,1] - 
(sabie_model$rep$fish_sel[y,,1,1] * sabie_model$rep$Fmort[y,1] +
sabie_model$rep$fish_sel[y,,1,2]  * sabie_model$rep$Fmort[y,2])



# Data weighting ----------------------------------------------------------
# Francis reweighting
# for(i in 1:10) {
# 
#   if(i == 1) { # set weights to 1 at first iteration
#     data$Wt_FishAgeComps[!is.na(data$Wt_FishAgeComps)] <- 1
#     data$Wt_FishLenComps[!is.na(data$Wt_FishLenComps)] <- 1
#     data$Wt_SrvAgeComps[!is.na(data$Wt_SrvAgeComps)] <- 1
#     data$Wt_SrvLenComps[!is.na(data$Wt_SrvLenComps)] <- 1
#   } # if i == 1 set weights to 1
# 
#   # make AD model function
#   sabie_model <- MakeADFun(data = data, parameters = parameters,
#                            map = mapping, random = NULL,
#                            DLL = "SabieTMB")
#   # Now, optimize the function
#   sabie_optim <- stats::nlminb(sabie_model$par, sabie_model$fn, sabie_model$gr,
#                                control = list(iter.max = 1e5, eval.max = 1e5))
#   # newton steps
#   try_improve <- tryCatch(expr =
#                             for(i in 1:2) {
#                               g = as.numeric(sabie_model$gr(sabie_optim$par))
#                               h = optimHess(sabie_optim$par, fn = sabie_model$fn, gr = sabie_model$gr)
#                               sabie_optim$par = sabie_optim$par - solve(h,g)
#                               sabie_optim$objective = sabie_model$fn(sabie_optim$par)
#                             }
#                           , error = function(e){e}, warning = function(w){w})
# 
#   sabie_model$optim <- sabie_optim # Save optimized model results
#   sabie_model$rep <- sabie_model$report(sabie_model$env$last.par.best) # Get report
# 
#   # Get new francis weights
#   new_weights <- francis_rwgt(data = data, model = sabie_model)
# 
#   data$Wt_FishAgeComps <- new_weights$new_fish_age_wts
#   data$Wt_FishLenComps <- new_weights$new_fish_len_wts
#   data$Wt_SrvAgeComps <- new_weights$new_srv_age_wts
#   data$Wt_SrvLenComps <- new_weights$new_srv_len_wts
# } # end i loop
# 
# sabie_model$sd_rep <- sdreport(sabie_model)

# Check consistency -------------------------------------------------------

# Get parameters
M_df <- data.frame(Par = "ln_M", 
                   TMB = sabie_model$sd_rep$par.fixed[names(sabie_model$sd_rep$par.fixed) == "ln_M"], 
                   ADMB = parameters$ln_M)

R0_df <- data.frame(Par = "ln_R0", 
                    TMB = sabie_model$sd_rep$par.fixed[names(sabie_model$sd_rep$par.fixed) == "ln_R0"], 
                    ADMB = parameters$ln_R0)

rec_sigmaRlate_df <- data.frame(Par = "ln_recSigma", 
                                TMB = sabie_model$sd_rep$par.fixed[names(sabie_model$sd_rep$par.fixed) == "ln_sigmaR_late"], 
                                ADMB = parameters$ln_sigmaR_late)

srv_q_df <- data.frame(Par = c('ln_srv_domLL_q', 'ln_srv_trwl_q', 'ln_srv_coopLL_q'),
                       TMB = sabie_model$sd_rep$par.fixed[names(sabie_model$sd_rep$par.fixed) == "ln_srv_q"],
                       ADMB = parameters$ln_srv_q[1,])

fish_q_df <- data.frame(Par = c('ln_fish_domLL_q1', 'ln_fish_domLL_q2', 'ln_fish_domLL_q3'),
                        TMB = sabie_model$sd_rep$par.fixed[names(sabie_model$sd_rep$par.fixed) == "ln_fish_q"],
                        ADMB = parameters$ln_fish_q[,1])

par_df <- rbind(M_df, R0_df, rec_sigmaRlate_df, srv_q_df, fish_q_df)

# Get time series
rec_series <- data.frame(Par = "Recruitment", 
                         Year = 1960:2022, 
                         TMB = rowSums(sabie_model$rep$NAA[-c(64,65),1,]), 
                         ADMB = rec[-64])

f_series <- data.frame(Par = "Total F", 
                         Year = 1960:2023, 
                         TMB = rowSums(sabie_model$rep$Fmort), 
                         ADMB = tem_dat$t.series$fmort)

ssb_series <- data.frame(Par = "SSB", 
                         Year = 1960:2023, 
                         TMB = sabie_model$rep$SSB, 
                         ADMB = ssb)

ts_df <- rbind(rec_series, f_series, ssb_series)

# Get selectivities
dom_ll_fish_f1 <- data.frame(Age = 1:30, 
                             TMB = sabie_model$rep$fish_sel[1,,1,1],
                             ADMB = tem_dat$agesel$fish1sel.f,
                             Type = "Domestic LL Fishery Female Block 1")

dom_ll_fish_m1 <- data.frame(Age = 1:30, 
                             TMB = sabie_model$rep$fish_sel[1,,2,1],
                             ADMB = tem_dat$agesel$fish1sel.m,
                             Type = "Domestic LL Fishery Male Block 1")

dom_ll_fish_f2 <- data.frame(Age = 1:30, 
                             TMB = sabie_model$rep$fish_sel[40,,1,1],
                             ADMB = tem_dat$agesel$fish4sel.f,
                             Type = "Domestic LL Fishery Female Block 2")

dom_ll_fish_m2 <- data.frame(Age = 1:30, 
                             TMB = sabie_model$rep$fish_sel[40,,2,1],
                             ADMB = tem_dat$agesel$fish4sel.m,
                             Type = "Domestic LL Fishery Male Block 2")

dom_ll_fish_f3 <- data.frame(Age = 1:30, 
                             TMB = sabie_model$rep$fish_sel[60,,1,1],
                             ADMB = tem_dat$agesel$fish5sel.f,
                             Type = "Domestic LL Fishery Female Block 3")

dom_ll_fish_m3 <- data.frame(Age = 1:30, 
                             TMB = sabie_model$rep$fish_sel[60,,2,1],
                             ADMB = tem_dat$agesel$fish5sel.m,
                             Type = "Domestic LL Fishery Male Block 3")

dom_trwl_fish_f <- data.frame(Age = 1:30, 
                              TMB = sabie_model$rep$fish_sel[1,,1,2],
                              ADMB = tem_dat$agesel$fish3sel.f,
                              Type = "Domestic Trawl Female")

dom_trwl_fish_m <- data.frame(Age = 1:30, 
                              TMB = sabie_model$rep$fish_sel[1,,2,2],
                              ADMB = tem_dat$agesel$fish3sel.m,
                              Type = "Domestic Trawl Male")

dom_ll_srv_f1 <- data.frame(Age = 1:30, 
                            TMB = sabie_model$rep$srv_sel[1,,1,1],
                            ADMB = tem_dat$agesel$srv1sel.f,
                            Type = "Domestic LL Survey Female Block 1")

dom_ll_srv_m1 <- data.frame(Age = 1:30, 
                            TMB = sabie_model$rep$srv_sel[1,,2,1],
                            ADMB = tem_dat$agesel$srv1sel.m,
                            Type = "Domestic LL Survey Male Block 1")

dom_ll_srv_f2 <- data.frame(Age = 1:30, 
                            TMB = sabie_model$rep$srv_sel[60,,1,1],
                            ADMB = tem_dat$agesel$srv10sel.f,
                            Type = "Domestic LL Survey Female Block 2")

dom_ll_srv_m2 <- data.frame(Age = 1:30, 
                            TMB = sabie_model$rep$srv_sel[60,,2,1],
                            ADMB = tem_dat$agesel$srv10sel.m,
                            Type = "Domestic LL Survey Male Block 2")

dom_trwl_srv_f2 <- data.frame(Age = 1:30, 
                              TMB = sabie_model$rep$srv_sel[60,,1,2],
                              ADMB = tem_dat$agesel$srv7sel.f,
                              Type = "Domestic Trawl Survey Female")

dom_trwl_srv_m2 <- data.frame(Age = 1:30, 
                              TMB = sabie_model$rep$srv_sel[60,,2,2],
                              ADMB = tem_dat$agesel$srv7sel.m,
                              Type = "Domestic Trawl Survey Male")

coop_ll_srv_f2 <- data.frame(Age = 1:30, 
                             TMB = sabie_model$rep$srv_sel[60,,1,3],
                             ADMB = tem_dat$agesel$srv2sel.f,
                             Type = "Coop LL Survey Female")

coop_ll_srv_m2 <- data.frame(Age = 1:30, 
                             TMB = sabie_model$rep$srv_sel[60,,2,3],
                             ADMB = tem_dat$agesel$srv2sel.m,
                             Type = "Coop LL Survey Male")

combined_sel <- rbind(
  dom_ll_fish_m1,
  dom_ll_fish_f2,
  dom_ll_fish_m2,
  dom_ll_fish_f3,
  dom_ll_fish_m3,
  dom_trwl_fish_f,
  dom_trwl_fish_m,
  dom_ll_srv_f1,
  dom_ll_srv_m1,
  dom_ll_srv_f2,
  dom_ll_srv_m2,
  dom_trwl_srv_f2,
  dom_trwl_srv_m2,
  coop_ll_srv_f2,
  coop_ll_srv_m2
)

# Plots -------------------------------------------------------------------

pdf(here("figs", "Comparison.pdf"), width = 15)

ggplot() +
  geom_line(ts_df, mapping  = aes(x = Year, y = TMB, color = "TMB"), size = 1.3) +
  geom_line(ts_df, mapping  = aes(x = Year, y = ADMB, color = "ADMB"), size = 1.3) +
  scale_color_manual(labels = c("TMB", "ADMB"), values = c("red", "black")) +
  facet_wrap(~Par, scales = "free") +
  theme_bw()

ggplot(ts_df, aes(x = Year, y = (TMB - ADMB) / ADMB, color = Par)) +
  geom_line(size = 1.3) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) +
  coord_cartesian(ylim = c(-0.05, 0.05)) +
  theme_bw()

ggplot(par_df, aes(x = Par, y = (exp(TMB) - exp(ADMB)) / exp(ADMB), group = Par)) +
  geom_point(size = 5) +
  geom_hline(yintercept = 0, lty = 2, size = 1.3) +
  coord_cartesian(ylim = c(-0.05, 0.05)) +
  labs(y = "(TMB - ADMB) / ADMB") +
  theme_bw()

ggplot() +
  geom_line(combined_sel, mapping = aes(x = Age , y = (TMB - ADMB)), lwd = 1.3) +
  facet_wrap(~Type) +   
  geom_hline(yintercept = 0, lty = 2, size = 1.3) +
  coord_cartesian(ylim = c(-0.05, 0.05)) +
  labs(y = "Selex (TMB - ADMB)") +
  theme_bw()

ggplot() +
  geom_line(combined_sel, mapping = aes(x = Age , y = TMB, color = "TMB")) +
  geom_line(combined_sel, mapping = aes(x = Age , y = ADMB, color = "ADMB")) +
  scale_color_manual(labels = c("TMB", "ADMB"), values = c("red", "black")) +
  facet_wrap(~Type) +
  labs(y = "Selex") +
  theme_bw()
dev.off()

