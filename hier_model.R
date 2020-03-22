library(tidyverse)
library(rstan)
library(bayesplot)
library(coda)
library(loo)

setwd("~/Desktop/school/dwyer_research/undergrad-thesis")
SOK_data_w_control <- read.csv(file="DoseR_SOK_reformatted.csv",
                               header=TRUE, sep=",", stringsAsFactors=FALSE)
SOK_data <- SOK_data_w_control %>% filter(capsid != "none")

factor_to_int <- function(factor, col_name) {
  factors <- t(unique(SOK_data[col_name]))
  dict <- as.list(1:length(factors))
  names(dict) <- factors
  return(dict[[factor]])
}

strain_factor_to_int = function(i) factor_to_int(i,"strain")
tree_factor_to_int = function(i) factor_to_int(i,"tree_sp")


## run stan model

N <- nrow(SOK_data) # number of treatments = 48
I <- length(unique(SOK_data$tree_sp)) # number of tree species = 2
J <- length(unique(SOK_data$strain)) # number of strains = 8
K <- length(unique(SOK_data$capsid)) # number of capsids = 2

tid <- sapply(SOK_data$tree_sp, tree_factor_to_int)
sid <- sapply(SOK_data$strain, strain_factor_to_int)
cid <- c(rep(1,4), rep(2,4))

x <- SOK_data$ob_count
y <- SOK_data$total_virus
total <- SOK_data$total_n

m_hier <- stan(file="hier_model.stan",
               data=list(N=N,I=I,J=J,K=K,
                         tid=tid,sid=sid,cid=cid,
                         x=x,y=y,total=total),
               chains=4,
               iter=2000,
               control = list(adapt_delta=0.90, max_treedepth=15))

### save fit
fit@stanmodel@dso <- new("cxxdso")
saveFit <- function(fit_obj, fit_name) {
  saveRDS(fit_obj, file=fit_name)
  fit <- fit_obj
  rm(fit_obj)
  return(fit)
}
fit <- saveFit(m_hier, "../stan_fits/fit_hier.rds")

### load fit
m_hier <- readRDS("../stan_fits/fit_hier.rds")
