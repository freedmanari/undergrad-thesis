# stan attempt (Sampling Through Adaptive Neighborhoods)

library(rstan)
library(tidyverse)
library(shinystan)
library(bayesplot)
library(coda)
library(loo)
library(pracma)

### load data
setwd("~/Desktop/school/dwyer_research/undergrad-thesis")
SOK_data <- read.csv(file="SOK_reordered.csv",
                     header=TRUE, sep=",", stringsAsFactors=FALSE)
y <- SOK_data$dose_response

### load fits
informative_priors <- readRDS("../stan_fits/informative_priors.rds")
no_means <- readRDS("../stan_fits/no_means.rds")

fit_hier <- informative_priors

# pulling out the fits
# NUTS (the No-U-Turn Sampler) optimizes HMC adaptively
posterior_fit_hier <- as.array(fit_hier)
nuts_params_fit_hier <- nuts_params(fit_hier)
log_posterior_fit_hier <- log_posterior(fit_hier)

summary(fit_hier, pars = c("sigma_alpha","sigma_beta"))$summary
head(log_posterior_fit_hier)

# trace plots for the runs
# to make sure nothing is being wonky
mcmc_trace(posterior_fit_hier, regex="alpha",np=nuts_params_fit_hier)
mcmc_trace(posterior_fit_hier, regex="beta",np=nuts_params_fit_hier)
mcmc_dens(posterior_fit_hier, regex="alpha")
mcmc_dens(posterior_fit_hier, regex="beta")

# looks at the divergence in the plots
color_scheme_set("red")
mcmc_nuts_divergence(nuts_params_fit_hier, log_posterior_fit_hier)

# calculates the rhats for the fit
# if they are above 1.05 that is bad
# still working on getting the Gelman-Rubin statistic
# that is just a split r-hat (I believe)
# last line plots the rhats
rhats <- rhat(fit_hier)
# print(rhats_simple)
mcmc_rhat(rhats)

# plots the fit and deviation of the parameters
plot(fit_hier, pars = "alpha")
plot(fit_hier, pars = "beta")
plot(fit_hier, pars = c("theta"), fill_color = "red", show_density=TRUE) +
  geom_point(data = data.frame(x_coords=logit(y),y_coords=48:1),
             mapping = aes(x = x_coords, y = y_coords),
             colour = "blue",
             shape = 4,
             size = 3)
plot(fit_hier, pars = c("inv_logit_theta"), fill_color = "red", show_density=TRUE) +
  geom_point(data = data.frame(x_coords=y,y_coords=48:1),
             mapping = aes(x = x_coords, y = y_coords),
             colour = "blue",
             shape = 4,
             size = 3) +
  geom_hline(yintercept=25)

plot(fit_hier, pars = "sigma_alpha", show_density=TRUE)
plot(fit_hier, pars = "sigma_beta", show_density=TRUE)

# gets the log likelihood matrix out and computes the waic
log_lik_informative_priors <- extract_log_lik(informative_priors, parameter_name = "log_lik", merge_chains=TRUE)
log_lik_no_means <- extract_log_lik(no_means, parameter_name = "log_lik", merge_chains=TRUE)
waic_informative_prios <- waic(log_lik_informative_priors)
waic_no_means <- waic(log_lik_no_means)

