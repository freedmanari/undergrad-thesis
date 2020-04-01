library(tidyverse)
library(rstan)
library(pracma)
library(loo)

set.seed(682932)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

setwd("~/Desktop/school/dwyer_research/undergrad-thesis")
SOK_data <- read.csv(file="SOK_reordered.csv",
                               header=TRUE, sep=",", stringsAsFactors=FALSE)

factor_to_int <- function(factor, col_name) {
  factors <- t(unique(SOK_data[col_name]))
  dict <- as.list(1:length(factors))
  names(dict) <- factors
  return(dict[[factor]])
}

strain_factor_to_int = function(i) factor_to_int(i,"strain")
tree_factor_to_int = function(i) factor_to_int(i,"tree_sp")



## initialize stan data

N <- nrow(SOK_data) # number of treatments = 48
I <- length(unique(SOK_data$tree_sp)) # number of tree species = 2
J <- length(unique(SOK_data$strain)) # number of strains = 8
K <- length(unique(SOK_data$capsid)) # number of capsids = 2

tid <- sapply(SOK_data$tree_sp, tree_factor_to_int) # GR = 1, DO = 2
sid <- sapply(SOK_data$strain, strain_factor_to_int) # COL=1,CUB=2,LOV=3,LST=4,DRY=5,KLP=6,TAM=7,TMB=8
cid <- c(rep(1,4), rep(2,4)) # SNPV = 1, MNPV = 2

x <- SOK_data$ob_count
y <- SOK_data$total_virus
total <- SOK_data$total_n


## get priors from glm with lowest AIC

data_MNPV_first <- SOK_data
data_SNPV_first <- SOK_data
data_SNPV_first$capsid <- factor(data_MNPV_first$capsid,levels=c("SNPV","MNPV"))

glm_MNPV_first <- glm(dose_response ~ ob_count * capsid + tree_sp * capsid,
                      family = "binomial", data = data_MNPV_first, weights=total_n)
glm_SNPV_first <- glm(dose_response ~ ob_count * capsid + tree_sp * capsid,
                      family = "binomial", data = data_SNPV_first, weights=total_n)

prior_mu_alpha <- matrix(nrow = J, ncol = I)
prior_mu_beta <- matrix(nrow = J, ncol = I)

p <- predict(glm_MNPV_first)
for (i in seq(1,N,3)) {
  prior_mu_beta[sid[i],tid[i]] <- (p[i+1] - p[i]) / x[i]
  prior_mu_alpha[sid[i],tid[i]] <- p[i] - prior_mu_beta[sid[i],tid[i]] * x[i]
}

prior_sigma_alpha <- c(Norm(coef(summary(glm_SNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]),
                       Norm(coef(summary(glm_MNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]))
prior_sigma_beta <- c(coef(summary(glm_SNPV_first))["ob_count", "Std. Error"],
                      coef(summary(glm_MNPV_first))["ob_count", "Std. Error"])



## run stan model

fit_hier <- stan(file="hier_model.stan",
                 data=list(N=N,I=I,J=J,K=K,
                           tid=tid,sid=sid,cid=cid,
                           x=x,y=y,total=total,
                           prior_mu_alpha=prior_mu_alpha,prior_mu_beta=prior_mu_beta,
                           prior_sigma_alpha=prior_sigma_alpha,prior_sigma_beta=prior_sigma_beta),
                 chains=4,
                 iter=10000,
                 control = list(adapt_delta=0.99, max_treedepth=25))

### save fits
fit@stanmodel@dso <- new("cxxdso")
saveFit <- function(fit_obj, fit_name) {
  saveRDS(fit_obj, file=fit_name)
  fit <- fit_obj
  rm(fit_obj)
  return(fit)
}
fit <- saveFit(fit_hier, "../stan_fits/fit_hier.rds")
