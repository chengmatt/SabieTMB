# Purpose: To make plots for sablefish assessment models
# Creator: Matthew LH. Cheng
# Date 5/6/23


# Set up ------------------------------------------------------------------

library(here)
library(R2admb)
library(tidyverse)
library(TMB)
library(rstan)
library(doSNOW)
library(parallel)
library(data.table)
library(reshape2)

theme_set(theme_bw() +
            theme(axis.text = element_text(size = 15, color = 'black'),
                  axis.title = element_text(size = 17),
                  strip.text = element_text(size = 17),
                  legend.title = element_text(size = 17),
                  legend.text = element_text(size = 15)))

# Load in sablefish model
setwd(here("src"))
compile("SabieTMB.cpp")
dyn.load(dynlib('SabieTMB'))

# set up cores to extract MCMC results
ncores <- detectCores() 
# Register cluster here
cl <- makeCluster(ncores - 2)
registerDoSNOW(cl)

model_res <- readRDS(here('output', 'Model_23.5', 'Model_23.5.rds')) # load in model
mcmc_res <- readRDS(here('output', 'Model_23.5', 'MCMC_Model_23.5.rds')) # load MCMC

# Abundance at Age and Year -----------------------------------------------

# Extract matrix and name columns
naa_df <- reshape2::melt(model_res$rep$NAA)
names(naa_df) <- c("Year", "Age", "Sex", "Numbers")
naa_df <- naa_df %>% mutate(Sex = ifelse(Sex == 1, "Female", "Male"), # differentiate sexes
                            Year = Year + 1959) 

ggplot(naa_df, aes(x = Year, y = Numbers, color = factor(Sex), lty = Sex)) +
  geom_line(size = 1, alpha = 0.85) +
  facet_wrap(~Age, scales = "free") +
  labs(x = "Year", y = "Numbers-at-age (millions)", color = "Sex") +
  theme(legend.position = 'top')

# Proportions at Age, Sex, and Year ---------------------------------------

# Summarize to get proportions
naa_prop_df <- naa_df %>% 
  group_by(Year) %>% 
  mutate(sum_prop = sum(Numbers)) %>% 
  ungroup() %>% 
  mutate(prop = Numbers / sum_prop)

ggplot(naa_prop_df %>% filter(Year %in% c(2005:2024)), aes(x = Age, y = prop,
                                                           color = Sex, fill = Sex, lty = Sex)) +
  geom_col(position = position_dodge()) +
  facet_wrap(~Year) +
  theme(legend.position = 'top') +
  labs(x = "Age", y = "Proportion", fill = "Sex")

# MCMC Extraction ---------------------------------------------------------

# Get MCMC derived posteriors
posterior <- as.matrix(mcmc_res) # posteriors
n_iters <- nrow(posterior) # number of iterations to loop through

# set up progress bar
pb <- txtProgressBar(max = n_iters, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Loop through each posterior draw, feed into the TMB report function, and
# regenerate report file from a given posterior draw to get derived quantities
derived_post <- foreach(iter = 1:n_iters, .packages = c("TMB", "here"), .options.snow = opts) %dopar% {
  dyn.load(dynlib('SabieTMB')) # Make sure sablefish model is availiable in parallel loop
  post_rep <- model_res$report(posterior[iter,-ncol(posterior)]) # removing last col because its the log posterior
  # output dataframe
  ssb_rec_post_df <- data.frame(Year = 1:length(post_rep$SSB), 
                                SSB = post_rep$SSB, Rec = post_rep$Rec,
                                iter = iter) 
  } # end for each loop
stopCluster(cl)

# Bind the parrallelized list together
derived_post_df <- data.table::rbindlist(derived_post)

# Summarize recruitment and ssb quantiles
derived_post_df <- derived_post_df %>% 
  mutate(Year = Year + 1959) %>%  # Adding years back
  group_by(Year) %>% 
  summarize(mean_SSB = mean(SSB),
            lwr_95_SSB = quantile(SSB, 0.025),
            upr_95_SSB = quantile(SSB, 0.975),
            mean_Rec = mean(Rec),
            lwr_95_Rec = quantile(Rec, 0.025),
            upr_95_Rec = quantile(Rec, 0.975))

# MCMC SSB ----------------------------------------------------------------

ggplot(derived_post_df, aes(x = Year, y = mean_SSB, ymin = lwr_95_SSB, ymax = upr_95_SSB)) +
  geom_line(size = 1.3) +
  geom_ribbon(alpha = 0.3) +
  labs(x = "Year", y = "Spawning Biomass (kt)")

# MCMC Recruitment --------------------------------------------------------

ggplot(derived_post_df, aes(x = Year, y = mean_Rec, ymin = lwr_95_Rec, ymax = upr_95_Rec)) +
  geom_line(size = 1.3) +
  geom_ribbon(alpha = 0.3) +
  labs(x = "Year", y = "Recruitment (millions of fish)")

