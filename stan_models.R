library(tidyverse)
library(rstan)
library(pracma)
library(loo)

set.seed(682932)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

SOK_data <- read.csv(file="SOK_reordered.csv",
                               header=TRUE, sep=",", stringsAsFactors=FALSE)

data_MNPV_first <- SOK_data
data_SNPV_first <- SOK_data
data_SNPV_first$capsid <- factor(data_MNPV_first$capsid,levels=c("SNPV","MNPV"))

factor_to_int <- function(factor, col_name) {
  factors <- t(unique(SOK_data[col_name]))
  dict <- as.list(1:length(factors))
  names(dict) <- factors
  return(dict[[factor]])
}

strain_factor_to_int = function(i) factor_to_int(i,"strain")
tree_factor_to_int = function(i) factor_to_int(i,"tree_sp")

N <- nrow(SOK_data) # number of treatments = 48
H <- length(unique(SOK_data$capsid)) # number of capsids = 2
I <- length(unique(SOK_data$strain)) # number of strains = 8
J <- length(unique(SOK_data$tree_sp)) # number of tree species = 2

cid <- c(rep(1,4),rep(2,4)) # SNPV = 1, MNPV = 2
sid <- sapply(SOK_data$strain, strain_factor_to_int) # COL=1,CUB=2,LOV=3,LST=4,DRY=5,KLP=6,TAM=7,TMB=8
tid <- sapply(SOK_data$tree_sp, tree_factor_to_int) # GR = 1, DO = 2
  
x <- SOK_data$ob_count
y <- SOK_data$total_virus
total <- SOK_data$total_n



## tree_and_strain

glm_SNPV_first <- glm(dose_response ~ ob_count * capsid * tree_sp,
                      family = "binomial", data = data_SNPV_first, weights=total_n)
glm_MNPV_first <- glm(dose_response ~ ob_count * capsid * tree_sp,
                      family = "binomial", data = data_MNPV_first, weights=total_n)

prior_mu_alpha <- matrix(nrow = I, ncol = J)
prior_mu_beta <- matrix(nrow = I, ncol = J)
p <- as.vector(predict(glm_SNPV_first))
for (n in seq(1,N,3)) {
  prior_mu_beta[sid[n],tid[n]] <- (p[n+1] - p[n]) / (x[n+1] - x[n])
  prior_mu_alpha[sid[n],tid[n]] <- p[n] - prior_mu_beta[sid[n],tid[n]] * x[n]
}

prior_sigma_alpha <- c(Norm(coef(summary(glm_SNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]),
                       Norm(coef(summary(glm_MNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]))
prior_sigma_beta <- c(Norm(coef(summary(glm_SNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"]),
                      Norm(coef(summary(glm_MNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"]))

tree_and_strain <- stan(file="tree_and_strain.stan",
                        data=list(N=N,H=H,I=I,J=J,
                                  cid=cid,sid=sid,tid=tid,
                                  x=x,y=y,total=total,
                                  prior_mu_alpha=prior_mu_alpha,prior_mu_beta=prior_mu_beta,
                                  prior_sigma_alpha=prior_sigma_alpha,prior_sigma_beta=prior_sigma_beta),
                        chains=4,
                        iter=5000,
                        control = list(adapt_delta=0.99, max_treedepth=25))
waic_tree_and_strain <- waic(extract_log_lik(tree_and_strain))




## tree_only

glm_SNPV_first <- glm(dose_response ~ ob_count * tree_sp,
                      family = "binomial", data = data_SNPV_first, weights=total_n)

prior_mu_alpha <- rep(0,J)
prior_mu_beta <- rep(0,J)
p <- as.vector(predict(glm_SNPV_first))
for (n in c(1,4)) {
  prior_mu_beta[tid[n]] <- (p[n+1] - p[n]) / (x[n+1] - x[n])
  prior_mu_alpha[tid[n]] <- p[n] - prior_mu_beta[tid[n]] * x[n]
}

prior_sigma_alpha <- Norm(coef(summary(glm_SNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"])
prior_sigma_beta <- Norm(coef(summary(glm_SNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"])

tree_only <- stan(file="tree_only.stan",
                  data=list(N=N,J=J,
                            tid=tid,
                            x=x,y=y,total=total,
                            prior_mu_alpha=prior_mu_alpha,prior_mu_beta=prior_mu_beta,
                            prior_sigma_alpha=prior_sigma_alpha,prior_sigma_beta=prior_sigma_beta),
                  chains=4,
                  iter=5000,
                  control = list(adapt_delta=0.99, max_treedepth=25))
waic_tree_only <- waic(extract_log_lik(tree_only))



## strain_only

glm_SNPV_first <- glm(dose_response ~ ob_count * capsid,
                      family = "binomial", data = data_SNPV_first, weights=total_n)
glm_MNPV_first <- glm(dose_response ~ ob_count * capsid,
                      family = "binomial", data = data_MNPV_first, weights=total_n)

prior_mu_alpha <- rep(0,I)
prior_mu_beta <- rep(0,I)
p <- as.vector(predict(glm_SNPV_first))
for (n in seq(1,N,6)) {
  prior_mu_beta[sid[n]] <- (p[n+1] - p[n]) / (x[n+1] - x[n])
  prior_mu_alpha[sid[n]] <- p[n] - prior_mu_beta[sid[n]] * x[n]
}

prior_sigma_alpha <- c(coef(summary(glm_SNPV_first))["(Intercept)", "Std. Error"],
                       coef(summary(glm_MNPV_first))["(Intercept)", "Std. Error"])
prior_sigma_beta <- c(coef(summary(glm_SNPV_first))["ob_count", "Std. Error"],
                      coef(summary(glm_MNPV_first))["ob_count", "Std. Error"])

strain_only <- stan(file="strain_only.stan",
                    data=list(N=N,H=H,I=I,
                              cid=cid,sid=sid,
                              x=x,y=y,total=total,
                              prior_mu_alpha=prior_mu_alpha,prior_mu_beta=prior_mu_beta,
                              prior_sigma_alpha=prior_sigma_alpha,prior_sigma_beta=prior_sigma_beta),
                    chains=4,
                    iter=5000,
                    control = list(adapt_delta=0.99, max_treedepth=25))
waic_strain_only <- waic(extract_log_lik(strain_only))


## no_tree_or_strain

glm_SNPV_first <- glm(dose_response ~ ob_count,
                      family = "binomial", data = data_SNPV_first, weights=total_n)

prior_mu_alpha <- coef(summary(glm_SNPV_first))["(Intercept)", "Estimate"]
prior_mu_beta <- coef(summary(glm_SNPV_first))["ob_count", "Estimate"]

prior_sigma_alpha <- coef(summary(glm_SNPV_first))["(Intercept)", "Std. Error"]
prior_sigma_beta <- coef(summary(glm_SNPV_first))["ob_count", "Std. Error"]

no_tree_or_strain <- stan(file="no_tree_or_strain.stan",
                          data=list(N=N,
                                    x=x,y=y,total=total,
                                    prior_mu_alpha=prior_mu_alpha,prior_mu_beta=prior_mu_beta,
                                    prior_sigma_alpha=prior_sigma_alpha,prior_sigma_beta=prior_sigma_beta),
                          chains=4,
                          iter=5000,
                          control = list(adapt_delta=0.99, max_treedepth=25))
waic_no_tree_or_strain <- waic(extract_log_lik(no_tree_or_strain))



## pooled_means

glm_SNPV_first <- glm(dose_response ~ ob_count * capsid * tree_sp,
                      family = "binomial", data = data_SNPV_first, weights=total_n)
glm_MNPV_first <- glm(dose_response ~ ob_count * capsid * tree_sp,
                      family = "binomial", data = data_MNPV_first, weights=total_n)

prior_mu_alpha <- matrix(nrow = H, ncol = J)
prior_mu_beta <- matrix(nrow = H, ncol = J)
p <- as.vector(predict(glm_SNPV_first))
for (n in c(1,4,25,28)) {
  prior_mu_beta[cid[sid[n]],tid[n]] <- (p[n+1] - p[n]) / (x[n+1] - x[n])
  prior_mu_alpha[cid[sid[n]],tid[n]] <- p[n] - prior_mu_beta[cid[sid[n]],tid[n]] * x[n]
}

prior_sigma_alpha <- c(Norm(coef(summary(glm_SNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]),
                       Norm(coef(summary(glm_MNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]))
prior_sigma_beta <- c(Norm(coef(summary(glm_SNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"]),
                      Norm(coef(summary(glm_MNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"]))

pooled_means <- stan(file="pooled_means.stan",
                        data=list(N=N,H=H,I=I,J=J,
                                  cid=cid,sid=sid,tid=tid,
                                  x=x,y=y,total=total,
                                  prior_mu_alpha=prior_mu_alpha,prior_mu_beta=prior_mu_beta,
                                  prior_sigma_alpha=prior_sigma_alpha,prior_sigma_beta=prior_sigma_beta),
                        chains=4,
                        iter=5000,
                        control = list(adapt_delta=0.99, max_treedepth=25))
waic_pooled_means <- waic(extract_log_lik(pooled_means))


## no_hierarchy

glm_SNPV_first <- glm(dose_response ~ ob_count * capsid * tree_sp,
                      family = "binomial", data = data_SNPV_first, weights=total_n)
glm_MNPV_first <- glm(dose_response ~ ob_count * capsid * tree_sp,
                      family = "binomial", data = data_MNPV_first, weights=total_n)

prior_mu_alpha <- matrix(nrow = I, ncol = J)
prior_mu_beta <- matrix(nrow = I, ncol = J)
p <- as.vector(predict(glm_SNPV_first))
for (n in seq(1,N,3)) {
  prior_mu_beta[sid[n],tid[n]] <- (p[n+1] - p[n]) / (x[n+1] - x[n])
  prior_mu_alpha[sid[n],tid[n]] <- p[n] - prior_mu_beta[sid[n],tid[n]] * x[n]
}

prior_sigma_alpha <- c(rep(Norm(coef(summary(glm_SNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]),4),
                       rep(Norm(coef(summary(glm_MNPV_first))[c("(Intercept)","tree_spGR"), "Std. Error"]),4))
prior_sigma_beta <- c(rep(Norm(coef(summary(glm_SNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"]),4),
                      rep(Norm(coef(summary(glm_MNPV_first))[c("ob_count","ob_count:tree_spGR"), "Std. Error"]),4))

no_hierarchy <- stan(file="no_hierarchy.stan",
                        data=list(N=N,I=I,J=J,
                                  sid=sid,tid=tid,
                                  x=x,y=y,total=total,
                                  prior_mu_alpha=prior_mu_alpha,prior_mu_beta=prior_mu_beta,
                                  prior_sigma_alpha=prior_sigma_alpha,prior_sigma_beta=prior_sigma_beta),
                        chains=4,
                        iter=5000,
                        control = list(adapt_delta=0.99, max_treedepth=25))
waic_no_hierarchy <- waic(extract_log_lik(no_hierarchy))



### save fits
fit@stanmodel@dso <- new("cxxdso")
saveFit <- function(fit_obj, fit_name) {
  saveRDS(fit_obj, file=fit_name)
  fit <- fit_obj
  rm(fit_obj)
  return(fit)
}

fit <- saveFit(tree_and_strain, "../stan_fits/tree_and_strain.rds")
fit <- saveFit(tree_only, "../stan_fits/tree_only.rds")
fit <- saveFit(strain_only, "../stan_fits/strain_only.rds")
fit <- saveFit(no_tree_or_strain, "../stan_fits/no_tree_or_strain.rds")
fit <- saveFit(pooled_means, "../stan_fits/pooled_means.rds")
fit <- saveFit(no_hierarchy, "../stan_fits/no_hierarchy.rds")



###
xs <- seq(0,4,.001)
df <- data.frame(x=xs,
                 y1=dnorm(xs,prior_sigma_alpha[1],3*prior_sigma_alpha[1]),
                 y2=dnorm(xs,prior_sigma_alpha[2],3*prior_sigma_alpha[2]))
ggplot(df) +
  geom_area(aes(x=x,y=y2),fill="#AA4400",alpha=.8) +
  geom_line(aes(x=x,y=y2)) +
  scale_x_continuous(expand=c(0,0),breaks = NULL) +
  scale_y_continuous(limits=c(0,.25),expand=c(0,0),breaks = NULL) +
  geom_segment(linetype="dashed",
               aes(x=prior_sigma_alpha[2],y=0,
                   xend=prior_sigma_alpha[2],
                   yend=dnorm(prior_sigma_alpha[2],prior_sigma_alpha[2],3*prior_sigma_alpha[2]))) +
  xlab("") +
  ylab("")



