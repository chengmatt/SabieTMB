# Purpose: To do a jitter analysis
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
sabie_model <- readRDS(here('output', 'Model_23.5', 'Model_23.5.RDS')) 

sd <- 0.4 # Specify CV/sd value from an rnorm
n_jitter <- 100 # number of jitter runsa
n_newton <- 4 # number of newton steps to take (extr)

# progress bar for parraellelization
pb <- txtProgressBar(max = n_jitter, style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# set up cores
ncores <- detectCores() 
# Register cluster here
cl <- makeCluster(ncores - 4)
registerDoSNOW(cl)

# Jitter Model ------------------------------------------------------------
jitter_all <- data.frame()
for(jit in 1:n_jitter) {
  
  # Jitter parameters
  jittered_pars <- sabie_model$par * exp(rnorm(length(sabie_model$par), 0, sd))
                     
  # Now, optimize the function
  sabie_optim <- stats::nlminb(jittered_pars,
                               sabie_model$fn, sabie_model$gr, 
                               control = list(iter.max = 1e5, eval.max = 1e5))
  # newton steps
  try_improve <- tryCatch(expr =
                            for(i in 1:n_newton) {
                              g = as.numeric(sabie_model$gr(sabie_optim$par))
                              h = optimHess(sabie_optim$par, fn = sabie_model$fn, gr = sabie_model$gr)
                              sabie_optim$par = sabie_optim$par - solve(h,g)
                              sabie_optim$objective = sabie_model$fn(sabie_optim$par)
                            }
                          , error = function(e){e}, warning = function(w){w})
  
  sabie_model$optim <- sabie_optim # Save optimized model results
  sabie_model$rep <- sabie_model$report(sabie_model$env$last.par.best) # Get report
  sabie_model$sd_rep <- sdreport(sabie_model) # Get sd report
  
  # put jitter results into a dataframe
  jitter_ts_df <- data.frame(Year = 1:length(sabie_model$rep$SSB), 
                          SSB = sabie_model$rep$SSB, Rec = sabie_model$rep$Rec, 
                          jitter = jit, Hessian = sabie_model$sd_rep$pdHess,
                          jnLL = sabie_model$rep$jnLL,
                          Max_Gradient = max(abs(sabie_model$sd_rep$gradient.fixed)))
  
  # # Jitter parameters
  jitter_all <- rbind(jitter_ts_df, jitter_all)
}


# Diagnostic Plots -------------------------------------------------------------------


# SSB Trends
ggplot(jitter_all, aes(x = Year + 1960, y = SSB, group = jitter)) +
  geom_line(size = 1.3, color = "grey21") +
  geom_line() +
  labs(x = "Year", y = "SSB (kt)", title = "Jitter Analysis SSB Trends") +
  theme_sablefish() 

ggsave(filename = here("figs", "Diagnostics", "Jitter_SSB.png"))

# Recruitment
ggplot(jitter_all, aes(x = Year + 1960, y = Rec, group = jitter)) +
  geom_line(size = 1.3, color = "grey21") +
  geom_line() +
  labs(x = "Year", y = "Age 2 Recruitment (millions)", title = "Jitter Analysis SSB Trends") +
  theme_sablefish() 

ggsave(filename = here("figs", "Diagnostics", "Jitter_Rec.png"))

# jnLL Plot
ggplot(jitter_all, aes(x = jitter, y = jnLL, color = Max_Gradient, shape = Hessian)) +
  geom_point(size = 5, alpha = 0.3) +
  geom_hline(yintercept = min(jitter_all$jnLL), lty = 2, size = 2, color = "blue") +
  facet_wrap(~Hessian, scales = 'free') +
  scale_color_viridis_c() +
  theme_sablefish() +
  theme(legend.position = "bottom") +
  guides(color = guide_colorbar(barwidth = 15, barheight = 0.5)) +
  labs(x = "Iteration", y = "Joint Negative Log-Liklihood", 
       color = "Maximum Gradient", title = "Jitter Analysis Joint Negative Log-Liklihood Values",
       shape = "Positive Definite Hessian")
  
ggsave(filename = here("figs", "Diagnostics", "Jitter_jnLL.png"), width = 11)
