# Purpose: To visualize the general dynamics of the genearlized spatial stock assessment model
# Creator: Matthew LH. Cheng
# Date 2/4/25


# Set up ------------------------------------------------------------------
source(here("R", "model", "dbeta_symmetric.R")) # symmetric beta
source(here("R", "model", "get_beta_scaled_pars.R")) # scaled

library(here)
library(tidyverse)
library(ggtern)
library(MCMCpack)


# Steepness (Prior) -------------------------------------------------------
h_mu <- c(0.25,0.35, 0.5, 0.95)
h_sd <- c(0.1, rep(0.1, 3))
vals <- seq(0.01,0.99,0.001) # values to profile across
scaled_beta_df <- data.frame()
scaled_beta_sim <- data.frame()

for(i in 1:length(h_sd)) {
  tmp_pars <- get_beta_scaled_pars(0.2, 1, h_mu[i], h_sd[i])
  tmp_ll <- dbeta((vals - 0.2) / (1 - 0.2),tmp_pars[1],tmp_pars[2],log=TRUE) # get beta likelihood
  scaled_beta_df <- rbind(scaled_beta_df, data.frame(x = vals, vals = tmp_ll, pars = paste("mu = ", h_mu[i], "sd = ", h_sd[i]))) # get scaled beta density
  tmp_sim <- rbeta(5e5,tmp_pars[1],tmp_pars[2]) * tmp_pars[4] + tmp_pars[3]
  scaled_beta_sim <- rbind(scaled_beta_sim, data.frame(vals = tmp_sim, pars = paste("mu = ", h_mu[i], "sd = ", h_sd[i])))
}

# histogram
scaled_beta_hist <- ggplot(scaled_beta_sim, aes(x = vals, color = factor(pars))) +
  geom_histogram() +
  labs(x = "x", y = "Count", color = "Parameter Combination") +
  facet_wrap(~pars, scales = 'free_y') +
  # xlim(0,1) +
  theme_bw(base_size = 15) +
  theme(legend.position = "top")

# simulated values
png(here("figs", "Model Structure", "h Beta Sim.png"), width = 800, height = 800)
scaled_beta_hist
dev.off()

# density
scaled_beta_dens <- ggplot(scaled_beta_df, aes(x = x, y = vals, color = factor(pars))) +
  geom_line(size = 2) +
  labs(x = "x", y = "log-likelihood", color = "Parameter Combination") +
  facet_wrap(~pars, scales = 'free') +
  theme_bw(base_size = 15) +
  theme(legend.position = "top")

png(here("figs", "Model Structure", "h Beta Density.png"), width = 800, height = 800)
scaled_beta_dens
dev.off()

# Movement (Prior) --------------------------------------------------------
n <- 10000
samples2 <- rdirichlet(n, c(2, 2, 2))  # slightly away from edges
samples1 <- rdirichlet(n, c(1.0, 1.0, 1.0))  # uniform
samples3 <- rdirichlet(n, c(0.5, 0.5, 0.5))  # concentrates on edges
samples4 <- rdirichlet(n, c(3, 3, 3))  # concentrates in center

# combine plots
df_all <- rbind(
  data.frame(samples3, alpha = "alpha = 0.5 (near edges)"),
  data.frame(samples1, alpha = "alpha = 1.0 (uniform)"),
  data.frame(samples2, alpha = "alpha = 2 (away from edges)"),
  data.frame(samples4, alpha = "alpha = 3 (closer to middle)")
) 

# plot dirichlet ternary
move_dir_prior_plot <- ggtern(df_all, aes(X1, X2, X3)) +
  geom_point(alpha = 0.1, size = 0.5) +
  facet_wrap(~alpha) +
  theme_bw(base_size = 15)

png(here("figs", "Model Structure", "Movement Dirichlet Prior.png"), width = 800, height = 800)
move_dir_prior_plot
dev.off()


# Tag Reporting (Prior) ---------------------------------------------------

# Symmetric Beta prior
vals <- seq(0,1,0.001) # values to profile across
beta_sd <- c(0.05, 0.1, 0.3, 0.5, 1) # sd to profile across
symm_beta_df <- data.frame()
for(i in 1:length(beta_sd)) {
  tmp <- dbeta_symmetric(vals,1,0,beta_sd[i],log=TRUE) # get beta likelihood
  symm_beta_df <- rbind(symm_beta_df, data.frame(x = vals, vals = tmp, sd = beta_sd[i]))
}

tagrep_beta_prior_plot <- ggplot(symm_beta_df, aes(x = x, y = vals, color = factor(sd))) +
  geom_line(size = 2) +
  labs(x = "x", y = "log-likelihood", color = "SD") +
  theme_bw(base_size = 15) +
  theme(legend.position = "top")

png(here("figs", "Model Structure", "TagRep Sym Beta Prior.png"), width = 800, height = 800)
tagrep_beta_prior_plot
dev.off()

# Full beta prior
vals <- seq(0.01,0.99,0.001) # values to profile across
beta_mu <- c(0.1, 0.1, 0.1, 0.1) # mu to proilfe across
beta_sd <- c(0.05, 0.1, 0.3, 0.5) # sd to profile across
full_beta_df <- data.frame()
for(i in 1:length(beta_sd)) {
  a = beta_mu[i] / beta_sd[i]^2 # alpha parameter
  b = (1 - beta_mu[i]) / beta_sd[i]^2 # beta parameter
  tmp <- dbeta(vals,a,b,log=TRUE) # get beta likelihood
  full_beta_df <- rbind(full_beta_df, data.frame(x = vals, vals = tmp, pars = paste("mu = ", beta_mu[i], "sd = ", beta_sd[i])))
}

tagrep_full_beta_prior_plot <- ggplot(full_beta_df, aes(x = x, y = vals, color = factor(pars))) +
  geom_line(size = 2) +
  labs(x = "x", y = "log-likelihood", color = "Parameter Combination") +
  facet_wrap(~pars) +
  theme_bw(base_size = 15) +
  theme(legend.position = "top")

png(here("figs", "Model Structure", "TagRep Full Beta Prior.png"), width = 800, height = 800)
tagrep_full_beta_prior_plot
dev.off()

