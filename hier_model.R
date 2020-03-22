library(tidyverse)
library(rstan)
library(bayesplot)
library(coda)
library(loo)

setwd("~/Desktop/school/dwyer_research/SOK")
SOK_data_w_control <- read.csv(file="DoseR_SOK_reformatted.csv",
                               header=TRUE, sep=",", stringsAsFactors=FALSE)
SOK_data <- SOK_data_w_control %>% filter(capsid != "none")

## Data reformatting to each row representing an individual caterpillar,
## with death outcome y = 1 if virus-killed and 0 otherwise

M <- nrow(SOK_data)    #number of distinct treatments (without control)
hier_data <- data.frame(capsid=integer(), # SNPV=1, MNPV=2
                 strain=integer(),        # COL=1, CUB=2, LOV=3, LST=4, DRY=5, KLP=6, TAM=7, TMB=8
                 tree_sp=integer(),       # GR=0, DO=1
                 ob_count=integer(),
                 y=integer())             # virus-killed = 1, otherwise = 0

factor_to_int <- function(row, col_name) {
  factors <- t(unique(SOK_data[col_name]))
  dict <- as.list(1:length(factors))
  names(dict) <- factors
  return(dict[[SOK_data[row, col_name]]])
}

for (i in 1:M) {
  virus <- SOK_data[i, "total_virus"]
  total <- SOK_data[i, "total_n"]
  if (virus > 0) {
    for (j in 1:virus) {
      hier_data <- add_row(hier_data,
                           capsid=factor_to_int(i, "capsid"),
                           strain=factor_to_int(i, "strain"),
                           tree_sp=factor_to_int(i, "tree_sp") - 1, # -1 for 0-based indexing
                           ob_count=SOK_data$ob_count[i],
                           y=1)
    }
  }
  for (j in 1:(total-virus)) {
    hier_data <- add_row(hier_data,
                         capsid=factor_to_int(i, "capsid"),
                         strain=factor_to_int(i, "strain"),
                         tree_sp=factor_to_int(i, "tree_sp") - 1, # -1 for 0-based indexing
                         ob_count=SOK_data$ob_count[i],
                         y=0)
  }
}


## run stan model

N <- nrow(hier_data) # number of caterpillars = 1319
J <- length(unique(hier_data$strain)) # number of strains = 8
K <- length(unique(hier_data$capsid)) # number of capsids = 2

sid <- hier_data$strain
cid <- c(rep(1,4), rep(2,4))

x_obs <- hier_data$ob_count
x_tree <- hier_data$tree_sp
y <- hier_data$y

m_hier <- stan(file="hier_model.stan",
               data=list(N=N,J=J,K=K,
                         sid=sid,cid=cid,
                         x_obs=x_obs,x_tree=x_tree,y=y),
               chains=4,
               iter=5000,
               control = list(adapt_delta=0.95, max_treedepth=15))

### save fit
fit@stanmodel@dso <- new("cxxdso")
saveFit <- function(fit_obj, fit_name) {
  saveRDS(fit_obj, file=fit_name)
  fit <- fit_obj
  rm(fit_obj)
  return(fit)
}
fit <- saveFit(m_hier, "fit_hier.rds")

### load fit
m_hier <- readRDS("fit_hier.rds")
