# Purpose: To simulate data for a spatially structured population
# This set up is a bit of a mess ... sorry
library(here)
# source(here("R", "functions", "Utility_Functions.R"))
# Set up
n_sims <- 1e3
n_yrs <- 20
n_regions <- 2
n_ages <- 5
n_lens <- 1e3
n_sexes <- 1
n_fish_fleets <- 1
n_srv_fleets <- 1
init_iter <- n_ages * 10

# Fishery
sigmaC <- 1e-3
Fmort <- array(0, dim = c(n_yrs, n_regions, n_fish_fleets, n_sims))
F_vec2 <- c(seq(0.20, 0.15, length.out = n_yrs / 2), seq(0.15, 0.15, length.out = (n_yrs / 2)))
F_vec1 <- c(seq(0.20, 0.15, length.out = n_yrs / 2), seq(0.15, 0.15, length.out = (n_yrs / 2)))

fish_sel <- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_fish_fleets, n_sims))

# loop through to propagate fishery stuff
for(y in 1:n_yrs) {
  for(r in 1:n_regions) {
    for(sim in 1:n_sims) {
      Fmort[y,1,1,sim] <- F_vec1[y] * exp(rnorm(1, 0, 0))
      Fmort[y,2,1,sim] <- F_vec2[y] * exp(rnorm(1, 0, 0))
    }
    for(s in 1:n_sexes) {
      for(f in 1:n_fish_fleets) {
        for(sim in 1:n_sims) fish_sel[y,r,,s,f,sim] <- 
        1 / (1 + exp(-1 * ((1:n_ages) - n_ages/2))) 
      }
    }
  }
}

# Biological
do_recruits_move <- 0 # recruits don't move == 0, move == 1
ln_rec_devs <- array(0, dim = c(n_yrs, n_regions, n_sims))
M <- array(0.5, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_sims))
Z <- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_sims))
rec_sexratio <- array(1, dim = c(n_yrs, n_regions, n_sexes, n_sims))
r0 <- array(0, dim = c(n_yrs, n_regions, n_sims))
r0[,1,] <- 150
r0[,2,] <- 150

# loop through to propagate biologicals
WAA <- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_sims))
Maturity_AA <- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_sims))
for(y in 1:n_yrs) {
  for(r in 1:n_regions) {
    for(s in 1:n_sexes) {
      if(s == 1) WAA[y,r,,s,] <- 5 * (1 - exp(-0.1 * 1:n_ages))
      if(s == 2) WAA[y,r,,s,] <- 3 * (1 - exp(-0.1 * 1:n_ages))
      Maturity_AA[y,r,,s,] <- 1 / (1 + exp(-0.3 * 1:n_ages))
    }
  }
}

init_sigmaR <- array(0.5, dim = c(1, n_regions))
sigmaR <- array(0.5, dim = c(n_yrs, n_regions))
recruitment_opt <- 0 # == 0 mean recruitment
recdev_opt <- 0 # == 0 global density dependence, == 1 local density dependence

# Set up movement matrix
movement_matrix <- array(0, dim = c(n_regions, n_regions, n_yrs, n_ages, n_sexes, n_sims)) # From, To 
for(y in 1:n_yrs) {
  for(a in 1:n_ages) {
    for(s in 1:n_sexes) {
      for(sim in 1:n_sims) {
        # diag(movement_matrix[,,y,a,s,sim]) <- 1
        movement_matrix[,,y,a,s,sim] <- c(0.6,0.5,0.4,0.5)
      } 
    } # end s loop
  } # end a loop
} # end y loop

# movement_matrix[] <- 1

# Create containers
Init_NAA <- array(0, dim = c(init_iter, n_regions, n_ages, n_sexes, n_sims))
NAA <- array(0, dim = c(n_yrs+1, n_regions, n_ages, n_sexes, n_sims))
SSB <- array(0, dim = c(n_yrs, n_regions, n_sims))
Total_Biom <- array(0, dim = c(n_yrs, n_regions, n_sims))
Obs_Catch <- array(0, dim = c(n_yrs, n_regions, n_fish_fleets, n_sims))
Obs_FishAgeComps <- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_fish_fleets, n_sims))
CAA <- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_fish_fleets, n_sims))
True_Catch <- array(0, dim = c(n_yrs, n_regions, n_fish_fleets, n_sims))
Obs_Catch <- array(0, dim = c(n_yrs, n_regions, n_fish_fleets, n_sims))

CAL <- array(0, dim = c(n_yrs, n_regions, n_lens, n_sexes, n_fish_fleets, n_sims))
Obs_SrvIdx <- array(0, dim = c(n_yrs, n_regions, n_srv_fleets, n_sims))
Srv_IAA <- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_srv_fleets, n_sims))
Srv_IAL <- array(0, dim = c(n_yrs, n_regions, n_lens, n_sexes, n_srv_fleets, n_sims))
srv_sel <- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_srv_fleets, n_sims))
srv_q <- array(1, dim = c(n_yrs, n_regions, n_srv_fleets, n_sims))
sigmaSrvIdx <- array(0.01, dim = c(n_regions, n_srv_fleets))

# loop through to propagate survey selex
for(y in 1:n_yrs) {
  for(r in 1:n_regions) {
    for(s in 1:n_sexes) {
      for(sf in 1:n_srv_fleets) {
        srv_sel[y,r,,s,sf,] <- 1 / (1 + exp(-1 * ((1:n_ages) - n_ages/2.5)))
      }
    }
  }
}

True_SrvIdx <- array(0, dim = c(n_yrs, n_regions, n_srv_fleets, n_sims))
Obs_SrvIdx <- array(0, dim = c(n_yrs, n_regions, n_srv_fleets, n_sims))
Obs_SrvAgeComps <- array(0, dim = c(n_yrs, n_regions, n_ages, n_sexes, n_srv_fleets, n_sims))

# Tagging stuff
n_tags <- 1e2
max_liberty <- 30
tag_years <- seq(1, n_yrs, 5)
n_tag_yrs <- length(tag_years)
n_tag_rel_events <- n_tag_yrs * n_regions
tag_rel_indicator <- expand.grid(regions = 1:n_regions, tag_yrs = tag_years) # get tag release indicator (by tag years and regions = a tag cohort)
Tag_Reporting <- array(0, dim = c(n_yrs, n_regions, n_sims))
Tag_Reporting[,1,] <- 0.2
Tag_Reporting[,2,] <- 0.2
Tag_Fish <- array(0, dim = c(n_tag_rel_events, n_ages, n_sexes, n_sims))
Tag_Ind_Mort <- array(0, dim = c(n_yrs, n_ages, n_sexes, n_sims))
Tag_Shed <- array(0, dim = c(n_yrs, n_ages, n_sexes, n_sims))
Tag_Releases <- array(0, dim = c(n_tag_rel_events, n_regions, n_ages, n_sexes, n_sims))
Tag_Avail <- array(0, dim = c(max_liberty + 1, n_tag_rel_events, n_regions, n_ages, n_sexes, n_sims))
Pred_Tag_Recap <- array(0, dim = c(max_liberty, n_tag_rel_events, n_regions, n_ages, n_sexes, n_sims))
Obs_Tag_Recap <- array(0, dim = c(max_liberty, n_tag_rel_events, n_regions, n_ages, n_sexes, n_sims))
comp_strc <- 2
# 0 = Split by sex and region
# 1 = split by region but not by sex
# 2 = joint by region sex 
comp_like <- 0 # multinomial likelihood

tag_like <- 2
# 0 = Poisson
# 1 = Negative Binomial
# 2 = Multinomial Release Conditioned
# 3 = Multinomial Recapture Conditioned
tag_nbiom_dispersion <- array(1, dim = c(n_ages, n_sexes)) # tagging dispersion for negative binomial
init_F <- array(0, dim = c(1, n_regions, n_fish_fleets, n_sims))

# Model Structure (Tagged population also follows this structure)
# 1) Initialize population
# 2) Recruitment + Release Tags (tag mortality happens instantaneously here)
# 3) Mortality and Ageing 
# 4) Movement
# 5) Tag Shedding

# set.seed(1)
for(sim in 1:n_sims) {
  for(y in 1:n_yrs) {
    
    # Initialize Age Structure ------------------------------------------------
    if(y == 1) {
      # Iterate equilibrium age structure with mortality, followed by movement
      Init_NAA[1,,1,,sim] = r0[1,,sim] * rec_sexratio[1,,,sim]
      for(i in 2:init_iter) {
        for(r in 1:n_regions) {
          for(s in 1:n_sexes) {
            for(f in 1:n_fish_fleets) {
              # Recruitment
              Init_NAA[i,r,1,s,sim] <- r0[1,r,sim] * rec_sexratio[1,r,s,sim]
              # Iterate from previous "virtual time step"
              Init_NAA[i,r,2:n_ages,s,sim] <- Init_NAA[i - 1,r,1:(n_ages-1),s,sim] *
                                              exp(-(M[1,r,1:(n_ages-1),s,sim] + init_F[1,r,f,sim] * fish_sel[1,r,1:(n_ages-1),s,f,sim])) # mortality
              # Accumulate plus group
              Init_NAA[i,r,n_ages,s,sim] <- Init_NAA[i,r,n_ages,s,sim] + # individuals from equilibrium that get aged a year
                                            Init_NAA[i-1,r,n_ages,s,sim] * # decrement and accumulate plus group
                                            exp(-(M[1,r,n_ages,s,sim] + init_F[1,r,f,sim] * fish_sel[1,r,n_ages,s,f,sim])) # mortality
            } # end f loop
          } # end s loop
        } # end r loop
        
        # Apply movement after iteration
        for(a in 1:n_ages) for(s in 1:n_sexes) Init_NAA[i,,a,s,sim] <- Init_NAA[i,,a,s,sim] %*% movement_matrix[,,1,a,s,sim]
        # Recruits don't move
        if(do_recruits_move == 0) {
          for(r in 1:n_regions) {
            for(s in 1:n_sexes) {
              Init_NAA[i,r,1,s,sim] = r0[1,r,sim] * rec_sexratio[1,r,s,sim]
            } # end s loop
          } # end r loop
        } # end if recruits don't move
      } # end i
      
      # Apply initial age structure deviations here (FLAG: Revise to incorporate more options)
      tmp_ln_init_devs <- NULL # Initialize container vector to allow for global recruitment
      for(r in 1:n_regions) {
        if(recdev_opt == 0 && is.null(tmp_ln_init_devs)) tmp_ln_init_devs <- rnorm(n_ages-2, 0, init_sigmaR[r]) # simulate initial deviations (global density dependence)
        if(recdev_opt == 1) tmp_ln_init_devs <- rnorm(n_ages-2, 0, init_sigmaR[r]) # simulate initial deviations (local density dependence)
        Init_NAA[init_iter,r,2:(n_ages-1),,sim] <- Init_NAA[init_iter,r,2:(n_ages-1),,sim] * rep(exp(tmp_ln_init_devs - init_sigmaR[r]^2/2), n_sexes) # apply deviations
        NAA[1,r,2:n_ages,,sim] <- Init_NAA[init_iter,r,2:n_ages,,sim] # Plug in initial age structure into 1st year (w/o recruitment)
      } # end r loop
      
    } # end initializing age structure
    
  # Run Annual Cycle --------------------------------------------------------
  tmp_ln_rec_devs <- NULL # Initialize container vector to allow for global recruitment (remains NULL within a given year)
  for(r in 1:n_regions) {
    if(recdev_opt == 0 && is.null(tmp_ln_rec_devs)) tmp_ln_rec_devs <- rnorm(1, 0, sigmaR[y,r]) # Get recruitment deviates (global density dependence)
    if(recdev_opt == 1) tmp_ln_rec_devs <- ln_rec_devs[y,r,sim] <- rnorm(1, 0, sigmaR[y,r]) # Get recruitment deviates (local density dependence)
    ln_rec_devs[y,r,sim] = tmp_ln_rec_devs # input vector of temporary rec devs
    for(s in 1:n_sexes) {
      # Recruitment (FLAG: Revise to incorporate more options)
      if(recruitment_opt == 0) NAA[y,r,1,s,sim] <- r0[y,r,sim] * exp(ln_rec_devs[y,r,sim] - sigmaR[y,r]^2/2) * rec_sexratio[y,r,s,sim] 
      # Mortality and Ageing
      for(a in 1:n_ages) {
        # Get total mortality
        Z[y,r,a,s,sim] <- M[y,r,a,s,sim] + sum(Fmort[y,r,,sim] * fish_sel[y,r,a,s,,sim]) # Z = M + Fmort
        if(a < n_ages) { 
          # Exponential mortality for individuals not in plus group (recruits experience mortality )
          NAA[y+1,r,a+1,s,sim] <- NAA[y,r,a,s,sim] * exp(-Z[y,r,a,s,sim])
        } else {
          # Accumulate individuals recently "recruited" into plus group and individuals from previous year
          NAA[y+1,r,n_ages,s,sim] <- NAA[y+1,r,n_ages,s,sim] + NAA[y,r,n_ages,s,sim] * exp(-Z[y,r,a,s,sim])
        } # end else (calculations for plus group)
      } # end a loop
      } # end s loop
    } # end r loop    

  # Movement following Tag Releases, Catch and Survey Accounting ------------
  # Apply movement after ageing processes
  for(a in 1:n_ages) for(s in 1:n_sexes) NAA[y+1,,a,s,sim] <- NAA[y+1,,a,s,sim] %*% movement_matrix[,,y,a,s,sim]
  # Recruits don't move (reapply movement)
  if(do_recruits_move == 0) for(r in 1:n_regions) NAA[y,r,1,,sim] <- r0[y,r,sim] * exp(ln_rec_devs[y,r,sim] - sigmaR[y,r]^2/2) * rec_sexratio[y,r,,sim] 
  } # end y loop
  
  # Get Biomass Calculations
  for(r in 1:n_regions) {
    for(y in 1:n_yrs) {
      Total_Biom[y,r,sim] <- sum(as.vector(NAA[y,r,,,sim]) * as.vector(WAA[y,r,,,sim])) # Total Biomass
      SSB[y,r,sim] <- sum(as.vector(NAA[y,r,,1,sim]) * as.vector(WAA[y,r,,1,sim]) * Maturity_AA[y,r,,1,sim]) # Spawning Stock Biomass 
      if(n_sexes == 1) SSB[y,r,sim] <- SSB[y,r,sim] * 0.5 # If single sex model, multiply SSB calculations by 0.5 
    } # end y loop
  } # end r loop
  
  # Catch Accounting --------------------------------------------------------
  for(y in 1:n_yrs) {
    for(f in 1:n_fish_fleets) {
      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {
          # Baranov's catch equation
          CAA[y,r,,s,f,sim] <- (Fmort[y,r,f,sim] * fish_sel[y,r,,s,f,sim]) / Z[y,r,,s,sim] *  NAA[y,r,,s,sim] * (1 - exp(-Z[y,r,,s,sim]))
          
          # Structuring composition data to be split by region and sex
          if(comp_strc == 0) {
            tmp_FishAgeComps_Prob <- CAA[y,r,,s,f,sim] # Get probabilities for a given region and sex
            if(comp_like == 0) Obs_FishAgeComps[y,r,,s,f,sim] <- array(rmultinom(1, 1e3, tmp_FishAgeComps_Prob / sum(tmp_FishAgeComps_Prob)), dim = c(dim(CAA[y,r,,s,f,sim]))) # simulate multinomial probabilities
          } # end if for "Split" approach for composition data (split by region and sex
          
        } # end s loop
        # Generate Catch Data
        True_Catch[y,r,f,sim] <- sum(CAA[y,r,,,f,sim] * WAA[y,r,,,sim]) # True Catch
        Obs_Catch[y,r,f,sim] <- True_Catch[y,r,f,sim] * exp(rnorm(1, 0, sigmaC)) # Observed Catch w/ lognormal deviations
        
        # Structuring composition data to be split by reigon 
        if(comp_strc == 1) {
          tmp_FishAgeComps_Prob <- CAA[y,r,,,f,sim] # Get probabilities for a given region and sex
          if(comp_like == 0) Obs_FishAgeComps[y,r,,,f,sim] <- array(rmultinom(1, 1e3, tmp_FishAgeComps_Prob / sum(tmp_FishAgeComps_Prob)), dim = c(dim(CAA[y,r,,,f,sim, drop = FALSE]))) # simulate multinomial probabilities
        } # end if for 'Split' approach for composition data split by region by not by sex
        
      } # end r loop
      
      # Structuring composition data to be joint across regions, ages, and sexes
      if(comp_strc == 2) { 
        # Store temporary probabilities ordered by ages, sexes, regions (i.e., age 1-30, sex 1, region 1, age 1-30, sex 2, region 1, 
        # age 1-30, sex 1, region 2, age 1-30, sex 2, region 2 ... )
        tmp_FishAgeComps_Prob <- aperm(CAA[y, , , , f, sim, drop = FALSE], c(3,4,2,1,5,6)) # ordered by ages, sexes, regions, year = y, fishery fleet = f, and sim = sim
        if(comp_like == 0) tmp_sim_FishAgeComps <- array(rmultinom(1, 1e3, tmp_FishAgeComps_Prob / sum(tmp_FishAgeComps_Prob)), dim = dim(tmp_FishAgeComps_Prob)) # Simulate Multinomial samples
        # Inputing simulated data into dataframe while reshaping to correct dimension (revert to year, region, ages, sexes, fleet, sim)
        Obs_FishAgeComps[y,,,,f,sim] <- aperm(tmp_sim_FishAgeComps, c(4,3,1,2,5,6)) 
      } # end if for "Joint" approach for composition data across regions, ages, and sexes
      
    } # end f loop
  } # end y loop
  
  # Survey Accounting -------------------------------------------------------
  for(y in 1:n_yrs) {
    for(sf in 1:n_srv_fleets) {
      for(r in 1:n_regions) {
        for(s in 1:n_sexes) {
          # Survey Ages Indexed
          Srv_IAA[y,r,,s,sf,sim] <- NAA[y,r,,s,sim] * srv_sel[y,r,,s,sf,sim] * exp(-0.5 * Z[y,r,,s,sim])
          
          # Structuring composition data to be split by region and sex
          if(comp_strc == 0) {
            tmp_SrvAgeComps_Prob <- Srv_IAA[y,r,,s,sf,sim] # Get probabilities for a given region and sex
            if(comp_like == 0) Obs_SrvAgeComps[y,r,,s,sf,sim] <- array(rmultinom(1, 1e3, tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)), dim = c(dim(Srv_IAA[y,r,,s,sf,sim]))) # simulate multinomial probabilities
          } # end if for "Split" approach for composition data (split by region and sex)
          
        } # end s loop
        True_SrvIdx[y,r,sf,sim] <- srv_q[y,r,sf,sim] * sum(Srv_IAA[y,r,,,sf,sim]) # True Survey Index
        Obs_SrvIdx[y,r,sf,sim] <- srv_q[y,r,sf,sim] * True_SrvIdx[y,r,sf,sim] * exp(rnorm(1, 0, sigmaSrvIdx[r,sf])) # Observed survey index w/ lognormal deviations
        
        # Structuring composition data to be split by reigon 
        if(comp_strc == 1) {
          tmp_SrvAgeComps_Prob <- Srv_IAA[y,r,,,sf,sim] # Get probabilities for a given region and sex
          if(comp_like == 0) Obs_SrvAgeComps[y,r,,,sf,sim] <- array(rmultinom(1, 1e3, tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)), dim = c(dim(Srv_IAA[y,r,,,sf,sim,drop=FALSE]))) # simulate multinomial probabilities
        } # end if for 'Split' approach for composition data split by region by not by sex
      } # end r loop
      
      # Structuring composition data to be joint across regions, ages, and sexes
      if(comp_strc == 2) { 
        # Store temporary Probabilities ordered by ages, sexes, regions 
        # (i.e., age 1-30, sex 1, region 1, age 1-30, sex 2, region 1, age 1-30, sex 1, region 2, age 1-30, sex 2, region 2 ... )
        tmp_SrvAgeComps_Prob <- aperm(Srv_IAA[y, , , , sf, sim, drop = FALSE], c(3,4,2,1,5,6)) # ordered by ages, sexes, regions, year = y, survey fleet = sf, and sim = sim
        if(comp_like == 0) tmp_sim_SrvAgeComps <- array(rmultinom(1, 1e3, tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)), dim = dim(tmp_SrvAgeComps_Prob)) # Simulate Multinomial samples
        # Inputing simulated data into dataframe while reshaping to correct dimension (revert to year, region, ages, sexes, fleet, sim)
        Obs_SrvAgeComps[y,,,,sf,sim] <- aperm(tmp_sim_SrvAgeComps, c(4,3,1,2,5,6)) 
      } # end if for "Joint" approach for composition data across regions, ages, and sexes
      
    } # end sf loop
  } # end y loop
    
  # Tag Releases and Recoveries ------------------------------------------------------------
  for(tag_rel in 1:n_tag_rel_events) {
    
    # Tag Releases
    tag_rel_yr <- tag_rel_indicator$tag_yrs[tag_rel] # get tag release year
    tag_rel_region <- tag_rel_indicator$regions[tag_rel] # get tag release region
    n_tags_rel <- round(Obs_SrvIdx[tag_rel_yr,,1,sim] / sum(Obs_SrvIdx[tag_rel_yr,,1,sim]) * n_tags) # distribute tags relative to regional survey abundance
    tmp_SrvAgeComps_Prob <- as.vector(Srv_IAA[tag_rel_yr, tag_rel_region, , , 1, sim]) # Use survey proportions to distribute tags
    tagged_fish <- (tmp_SrvAgeComps_Prob / sum(tmp_SrvAgeComps_Prob)) * n_tags_rel[tag_rel_region] # Distribute tags across ages, sexes, for a given release event
    Tag_Fish[tag_rel,,,sim] <- array(tagged_fish, dim = c(n_ages, n_sexes)) # Reshape format by ages, and sexes and input into other array
    
    # Tag Recaptures (do we need to pool individuals into a plus max liberty group or can we just ignore?)
    for(recap_yr in 1:min(max_liberty, n_yrs - tag_rel_yr + 1)) { # recapture year for a given cohort - adding a cut off for time at liberty so we are not tracking fish for 100+ years
      actual_yr <- tag_rel_yr + recap_yr - 1 # Define actual year for indexing purposes
      # Input tagged fish into available tags for recapture and adjust initial number of tagged fish for tag induced mortality (exponential mortality process)
      if(recap_yr == 1) Tag_Avail[1,tag_rel,tag_rel_region,,,sim] <- Tag_Fish[tag_rel,,,sim] * exp(-Tag_Ind_Mort[actual_yr,,,sim])
      
      # Get mortality estimates
      if(n_fish_fleets == 1) tmp_F <- Fmort[actual_yr,,,sim] # temporary fishing mortality variable (uniform sel) 
      if(n_fish_fleets > 1) tmp_F <- rowSums(Fmort[actual_yr,,,sim]) # temporary fihsing mortality variable (uniform sel) 
      tmp_Z <- M[actual_yr,,,,sim, drop = FALSE] + tmp_F  # temporary total mortality variable
      
      # Apply mortality and ageing to tagged fish
      for(a in 1:n_ages) {
        for(s in 1:n_sexes) {
          if(a < n_ages) { # If not in plus group
            # same dynamics as all other fish, but with uniform selex
            Tag_Avail[recap_yr+1,tag_rel,,a+1,s,sim] <- Tag_Avail[recap_yr,tag_rel,,a,s,sim]  * exp(-tmp_Z[1,,a,s,1]) 
          } else{ # Accumulate plus group here
            Tag_Avail[recap_yr+1,tag_rel,,n_ages,s,sim] <- Tag_Avail[recap_yr+1,tag_rel,,n_ages,s,sim] +  
                                                           Tag_Avail[recap_yr,tag_rel,,n_ages,s,sim]  * exp(-tmp_Z[1,,n_ages,s,1])
          } # end else for in plus group
        } # end s loop
      } # end a loop
      
      # Get predicted recaptures (Baranov's catch equation for tagged cohorts * tag reporting rate)
      # Because the dynamics are release tags -> tag induced mortality, ageing + mortality, and then movement,
      # the predicted recaptures for recapture year 1 only provides expected values for the release region. In the next year,
      # because individuals have undergone movement dynamics, recaptures are provided for other regions. (i.e., recapture year 1 only
      # provides information on mortality, while all subsequent years provide information on mortality and movement). So the question is,
      # what to do with recapture year 1? Do we fit these or just ignore?
      Pred_Tag_Recap[recap_yr,tag_rel,,,,sim] <- Tag_Reporting[actual_yr,,sim] * (tmp_F / tmp_Z[1,,,,1]) * 
                                                 Tag_Avail[recap_yr,tag_rel,,,,sim] * (1 - exp(-tmp_Z[1,,,,1])) 
      
      # Move tagged fish around 
      for(a in 1:n_ages) for(s in 1:n_sexes) Tag_Avail[recap_yr+1,tag_rel,,a,s,sim] <- 
                                             Tag_Avail[recap_yr+1,tag_rel,,a,s,sim] %*% movement_matrix[,,actual_yr,a,s,sim]
      
      # Apply tag shedding after movement occurs (just a mortality process)
      for(r in 1:n_regions) Tag_Avail[recap_yr+1,tag_rel,r,,,sim] <- Tag_Avail[recap_yr+1,tag_rel,r,,,sim] * exp(-Tag_Shed[actual_yr,,,sim]) 


      # Tag Recapture Likelihoods -----------------------------------------------
      # Simulate observed tag recoveries
      # Poisson tag recovery
      for(r in 1:n_regions) {
        for(a in 1:n_ages) {
          for(s in 1:n_sexes) {
            if(tag_like == 0) Obs_Tag_Recap[recap_yr,tag_rel,r,a,s,sim] <- rpois(1, Pred_Tag_Recap[recap_yr,tag_rel,r,a,s,sim]) # Poisson tag recovery
            if(tag_like == 1) Obs_Tag_Recap[recap_yr,tag_rel,r,a,s,sim] <- rnbinom(1, mu = Pred_Tag_Recap[recap_yr,tag_rel,r,a,s,sim], size = tag_nbiom_dispersion[a,s]) # Negbin tag recovery
          } # end s loop
        } # end a loop
      } # end r loop
      
      # Multinomial tag recovery (release conditioned)
      if(tag_like == 2) { 
      tmp_n_tags_rel <- sum(Tag_Fish[tag_rel,,,sim]) # Number of initial tags released
      tmp_recap <- aperm(Pred_Tag_Recap[recap_yr,tag_rel,,,,sim, drop = FALSE] / tmp_n_tags_rel, c(4,5,3,1,2,6)) # get recapture probabilities ordered by ages, sexes, regions (i.e., age 1-30, sex 1, region 1, age 1-30, sex 2, region 1, age 1-30, sex 1, region 2, age 1-30, sex 2, region 2 ... ) recapture, tag release, sim ... 
      tmp_non_recap <- (tmp_n_tags_rel - sum(Pred_Tag_Recap[1:recap_yr,tag_rel,,,,sim])) / tmp_n_tags_rel # calculate non recapture probabilities (tags releases - sum of all tags recaptured up to a given recapture year)
      tmp_recap_non_recap_probs <- c(tmp_recap, tmp_non_recap) # concatenate recapture and non recapture probabilities
      tmp_sim_recap <- rmultinom(1, tmp_n_tags_rel, tmp_recap_non_recap_probs) # simulate multinomial draws here
      tmp_sim_recap <- aperm(array(tmp_sim_recap[-length(tmp_sim_recap)], dim(tmp_recap)), c(4,5,3,1,2,6)) # remove last group (not recaptured) and then reshape into correct format
      Obs_Tag_Recap[recap_yr,tag_rel,,,,sim] <- tmp_sim_recap # input recaptures from multinomial into observed array
      } # end if for multinomial likelihood (release conditioned)
      
      # Multinomial tag recovery (recovery conditioned)
      if(tag_like == 3) { # Tag reporting doesn't matter here (if spatially invariant) since it appears in the denominator so it cancels out (when calculating proportion of recaptures)
        tmp_n <- sum(Pred_Tag_Recap[recap_yr,tag_rel,,,,sim]) # Get number of tags to simulate
        tmp_recap_probs <- aperm(Pred_Tag_Recap[recap_yr,tag_rel,,,,sim, drop = FALSE] / sum(Pred_Tag_Recap[recap_yr,tag_rel,,,,sim, drop = FALSE]), c(4,5,3,1,2,6)) # get recapture probabilities ordered by ages, sexes, regions (i.e., age 1-30, sex 1, region 1, age 1-30, sex 2, region 1, age 1-30, sex 1, region 2, age 1-30, sex 2, region 2 ... ) recapture, tag release, sim ... 
        tmp_sim_recap <- rmultinom(1, tmp_n, tmp_recap_probs) # simulate multinomial draws here
        tmp_sim_recap <- aperm(array(tmp_sim_recap, dim(tmp_recap_probs)), c(4,5,3,1,2,6)) # reshape into correct format
        Obs_Tag_Recap[recap_yr,tag_rel,,,,sim] <- tmp_sim_recap # input recaptures from multinomial into observed array
      } # end if for Multinomial likelihood (recovery conditioned)
      
    } # end recap_yr loop
  } # end tag_rel loop
} # end sim loop

# Output simulation outputs as a list
sim_out <- list(init_F = init_F,
     Fmort = Fmort,
     fish_sel = fish_sel,
     ln_rec_devs = ln_rec_devs,
     M = M,
     Z = Z,
     rec_sexratio = rec_sexratio,
     r0 = r0,
     WAA = WAA,
     Maturity_AA = Maturity_AA,
     init_sigmaR = init_sigmaR,
     sigmaR = sigmaR,
     movement_matrix = movement_matrix,
     Init_NAA = Init_NAA,
     NAA = NAA,
     SSB = SSB,
     Total_Biom = Total_Biom,
     True_Catch = True_Catch,
     Obs_Catch = Obs_Catch,
     CAA = CAA,
     Obs_FishAgeComps = Obs_FishAgeComps,
     CAL = CAL,
     Obs_SrvIdx = Obs_SrvIdx,
     True_SrvIdx = True_SrvIdx,
     Srv_IAA = Srv_IAA,
     Srv_IAL = Srv_IAL,
     srv_sel = srv_sel,
     srv_q = srv_q,
     Obs_SrvAgeComps = Obs_SrvAgeComps,
     Tag_Reporting = Tag_Reporting,
     Tag_Fish = Tag_Fish,
     Tag_Ind_Mort = Tag_Ind_Mort,
     Tag_Shed = Tag_Shed,
     Tag_Releases = Tag_Releases,
     Tag_Avail = Tag_Avail,
     Pred_Tag_Recap = Pred_Tag_Recap,
     Obs_Tag_Reca = Obs_Tag_Recap)

saveRDS(sim_out, file = here("sim_out.RDS"))

par(mfrow = c(2,1))
for(i in 1:n_sims) {
  if(i == 1) plot(SSB[,1,i], type = 'l') 
  else lines(SSB[,1,i], type = 'l') 
}
# for(i in 1:n_sims) {
#   if(i == 1) plot(SSB[,2,i], type = 'l') 
#   else lines(SSB[,2,i], type = 'l') 
# }


dev.off()

