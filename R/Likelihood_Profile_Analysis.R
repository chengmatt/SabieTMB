# Purpose: To do a likelihood profile analysis
# Creator: Matthew LH. Cheng
# Date 5/13/24

# Set up ------------------------------------------------------------------
library(here)
library(tidyverse)
library(TMB)
library(doSNOW)
library(parallel)

setwd(here("src"))
compile("SabieTMB.cpp")
dyn.load(dynlib('SabieTMB'))

# Load in model
sabie_model <- readRDS(here('output', 'Model_23.5', 'Model_23.5.rds')) 

# Do Profile --------------------------------------------------------------
