# Purpose: To run MCMC on the sablefish stock assessment model using STAN
# Creator: Matthew LH. Cheng (UAF-CFOS)
# Date 5/6/24

# Set up ------------------------------------------------------------------

library(tmbstan)
library(RTMB)
library(here)
library(tidyverse)
library(parallel)
library(shinystan)

# Set up MCMC options 
cores <- parallel::detectCores()-6 # Detect number of cores to use/number of independent chains to run 
iter <- 500 # Number of iterations to run per chain (one million draws)
thin <- 1 # Thinning rate  

# Load in model
sabie_model <- readRDS(here('output', 'Model_23.5', 'Model_23.5.rds')) 

# Run MCMC ----------------------------------------------------------------
options(mc.cores = cores) # Set up parrallelization - doesn't work on multiple cores for
# RTMB if you don't define data and parameters explicitly 
fit <- tmbstan::tmbstan(obj = sabie_rtmb_model, iter = iter, thin = thin, 
                        chains = 2, open_progress=FALSE) # Run MCMC

saveRDS(fit, here('output', 'Model_23.5', 'MCMC_Model_23.5.rds'))

# Launch shiny stan to inspect other MCMC diagnostics
launch_shinystan(fit)

