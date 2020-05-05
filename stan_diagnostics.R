# stan attempt (Sampling Through Adaptive Neighborhoods)

library(rstan)
library(tidyverse)
library(shinystan)
library(bayesplot)
library(coda)
library(loo)
library(pracma)

### load data
SOK_data <- read.csv(file="SOK_reordered.csv",
                     header=TRUE, sep=",", stringsAsFactors=FALSE)
y <- SOK_data$dose_response

### load fits
tree_and_strain <- readRDS("../stan_fits/tree_and_strain.rds")
tree_only <- readRDS("../stan_fits/tree_only.rds")
strain_only <- readRDS("../stan_fits/strain_only.rds")
no_tree_or_strain <- readRDS("../stan_fits/no_tree_or_strain.rds")
pooled_means <- readRDS("../stan_fits/pooled_means.rds")
no_hierarchy <- readRDS("../stan_fits/no_hierarchy.rds")

bhms <- c(tree_and_strain,no_hierarchy,pooled_means,strain_only,tree_only,no_tree_or_strain)
name <- c("B4","B5","B6","B2","B3","B1")
looic <- sapply(bhms,function(model) loo(model)$looic)
delta_looic <- looic - min(looic)
weight <- exp(-looic/2) / sum(exp(-looic/2))
p_looic <- sapply(bhms,function(model) loo(model)$p_loo)
lhood <- (p_looic-looic)/2
looic_table <- data.frame(name,lhood,looic,delta_looic,weight,p_looic)


fit_hier <- tree_and_strain

# pulling out the fits
# NUTS (the No-U-Turn Sampler) optimizes HMC adaptively
posterior_fit_hier <- as.array(fit_hier)
nuts_params_fit_hier <- nuts_params(fit_hier)
log_posterior_fit_hier <- log_posterior(fit_hier)

#summary(fit_hier, pars = c("sigma_alpha","sigma_beta"))$summary
#head(log_posterior_fit_hier)

# trace plots for the runs
# to make sure nothing is being wonky
mcmc_trace(posterior_fit_hier, regex="alpha",np=nuts_params_fit_hier)
mcmc_trace(posterior_fit_hier, regex="beta",np=nuts_params_fit_hier)



mcmc_dens(posterior_fit_hier, regex="^alpha.{2,}1",facet_args=list(ncol=8)) +
  scale_x_continuous(limits=c(-6,2),expand = c(0,0)) +
  theme(text=element_text(size=11,family="Palatino"),
        panel.spacing= unit(1.2, "lines")) +
  facet_text(FALSE)
mcmc_dens(posterior_fit_hier, regex="^alpha.{2,}2",facet_args=list(ncol=8)) +
  scale_x_continuous(limits=c(-6,2),expand = c(0,0)) +
  theme(text=element_text(size=11,family="Palatino"),
        panel.spacing= unit(1.2, "lines")) +
  facet_text(FALSE)

mcmc_dens(posterior_fit_hier, regex="^beta.{2,}1",facet_args=list(ncol=8)) +
  scale_x_continuous(limits=c(-.0004,.0006),expand = c(0,0),
                     breaks=c(-.0004,-.0002,0,.0002,.0004,.0006),
                     labels=c("-.0004","","0","",".0004","")) +
  theme(text=element_text(size=11,family="Palatino"),
        strip.background = element_blank(),
        panel.spacing= unit(1.2, "lines")) +
  facet_text(FALSE)
mcmc_dens(posterior_fit_hier, regex="^beta.{2,}2",facet_args=list(ncol=8)) +
  scale_x_continuous(limits=c(-.0004,.0006),expand = c(0,0),
                     breaks=c(-.0004,-.0002,0,.0002,.0004,.0006),
                     labels=c("-.0004","","0","",".0004","")) +
  theme(text=element_text(size=11,family="Palatino"),
        strip.background = element_blank(),
        panel.spacing= unit(1.2, "lines")) +
  facet_text(FALSE)


mcmc_dens(posterior_fit_hier, regex="sigma_alpha") +
  scale_x_continuous(limits=c(0,3),expand = c(0,0)) +
  theme(text=element_text(size=11,family="Palatino"),
        panel.spacing= unit(1.2, "lines"),
        plot.margin=margin(0,20,0,10)) +
  facet_text(FALSE)
mcmc_dens(posterior_fit_hier, regex="sigma_beta") +
  scale_x_continuous(limits=c(0,.0006),breaks=c(0,.0002,.0004,.0006),
                     labels=c("0",".0002",".0004",".0006"),expand = c(0,0)) +
  theme(text=element_text(size=11,family="Palatino"),
        panel.spacing= unit(1.2, "lines"),
        plot.margin=margin(0,20,0,10)) +
  facet_text(FALSE)



labs <- c("α_COL,GR","α_CUB,GR","α_LOV,GR","α_LST,GR",
               "α_DRY,GR","α_KLP,GR","α_TAM,GR","α_TMB,GR",
               "α_COL,DO","α_CUB,DO","α_LOV,DO","α_LST,DO",
               "α_DRY,DO","α_KLP,DO","α_TAM,DO","α_TMB,DO")
names(labs) <- c("alpha[1,1]","alpha[2,1]","alpha[3,1]","alpha[4,1]",
                      "alpha[5,1]","alpha[6,1]","alpha[7,1]","alpha[8,1]",
                      "alpha[1,2]","alpha[2,2]","alpha[3,2]","alpha[4,2]",
                      "alpha[5,2]","alpha[6,2]","alpha[7,2]","alpha[8,2]")



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

f <- function(i) return(paste("inv_logit_theta[",i,"]",sep=""))
plot(fit_hier, pars = sapply(48:1,f), fill_color = "red", show_density=TRUE) +
  geom_point(data = data.frame(x_coords=y,y_coords=1:48),
             mapping = aes(x = x_coords, y = y_coords),
             colour = "blue",
             shape = 4,
             size = 3) +
  geom_hline(yintercept=25) +
  geom_hline(yintercept=c(7,13,19,31,37,43),linetype="dashed")

plot(fit_hier, pars = "sigma_alpha", show_density=TRUE)
plot(fit_hier, pars = "sigma_beta", show_density=TRUE)

# gets the log likelihood matrix out and computes the waic
log_lik_informative_priors <- extract_log_lik(informative_priors, parameter_name = "log_lik", merge_chains=TRUE)
log_lik_no_means <- extract_log_lik(no_means, parameter_name = "log_lik", merge_chains=TRUE)
waic_informative_prios <- waic(log_lik_informative_priors)
waic_no_means <- waic(log_lik_no_means)



# K-S tests
sigma_alpha_SNPV <- c(sapply(1:4,function(i) tree_and_strain@sim$samples[[1]]$`sigma_alpha[1]`))
sigma_alpha_MNPV <- c(sapply(1:4,function(i) tree_and_strain@sim$samples[[1]]$`sigma_alpha[2]`))
sigma_beta_SNPV <- c(sapply(1:4,function(i) tree_and_strain@sim$samples[[1]]$`sigma_beta[1]`))
sigma_beta_MNPV <- c(sapply(1:4,function(i) tree_and_strain@sim$samples[[1]]$`sigma_beta[2]`))

ks.test(sigma_alpha_SNPV,sigma_alpha_MNPV,alternative="greater")
ks.test(sigma_beta_SNPV,sigma_beta_MNPV,alternative="greater")

t.test(sigma_alpha_SNPV,sigma_alpha_MNPV)
t.test(sigma_beta_SNPV,sigma_beta_MNPV)

