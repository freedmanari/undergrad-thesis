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
SOK_data$strain <- factor(SOK_data$strain,levels=c("COL","CUB","LOV","LST","DRY","KLP","TAM","TMB"))
SOK_data$capsid <- factor(SOK_data$capsid,levels=c("SNPV","MNPV"))
SOK_data$tree_sp <- factor(SOK_data$tree_sp,levels=c("GR","DO"))
y <- SOK_data$dose_response

### load fits
morphotype_and_tree <- readRDS("../stan_fits/morphotype_and_tree.rds")
tree_only <- readRDS("../stan_fits/tree_only.rds")
morphotype_only <- readRDS("../stan_fits/morphotype_only.rds")
neither_morphotype_nor_tree <- readRDS("../stan_fits/neither_morphotype_nor_tree.rds")
complete_hierarchy <- readRDS("../stan_fits/complete_hierarchy.rds")
no_hierarchy <- readRDS("../stan_fits/no_hierarchy.rds")

bhms <- c(morphotype_and_tree,
          tree_only,morphotype_only,neither_morphotype_nor_tree,no_hierarchy,complete_hierarchy)
name <- c("M and T","T only","M only","neither M nor T","No hierarchy","Complete hierarchy")
looic <- sapply(bhms,function(model) loo(model)$looic)
delta_looic <- looic - min(looic)
weight <- exp(-looic/2) / sum(exp(-looic/2))
p_looic <- sapply(bhms,function(model) loo(model)$p_loo)
lhood <- (p_looic-looic)/2
looic_table <- data.frame(name,lhood,looic,delta_looic,weight,p_looic)


fit_hier <- morphotype_and_tree

### plot of best model

factor_to_int <- function(factor, col_name) {
  factors <- t(unique(SOK_data[col_name]))
  dict <- as.list(1:length(factors))
  names(dict) <- factors
  return(dict[[factor]])
}

strain_factor_to_int = function(i) factor_to_int(i,"strain")
tree_factor_to_int = function(i) factor_to_int(i,"tree_sp")

estimates <- summary(fit_hier)$summary[,"mean"]
alphas <- estimates[str_detect(names(estimates),pattern="^alpha")]
betas <- estimates[str_detect(names(estimates),pattern="^beta")]

doses <- seq(0,6500,100)
B4_predictions <- data.frame(dose=numeric(16*length(doses)),
                             strain=character(16*length(doses)),
                             tree_sp=character(16*length(doses)),
                             prediction=numeric(16*length(doses)),
                             stringsAsFactors = FALSE)
i <- 1
for (s in unique(SOK_data[,"strain"])) {
  s_int <- strain_factor_to_int(s)
  for (t in unique(SOK_data[,"tree_sp"])) {
    t_int <- tree_factor_to_int(t)
    for (d in doses) {
      B4_predictions[i,] <- list(dose=d,
                                 strain=s,
                                 tree_sp=t,
                                 prediction = alphas[paste("alpha[",s_int,",",t_int,"]",sep="")] +
                                              betas[paste("beta[",s_int,",",t_int,"]",sep="")] * d)
      i <- i+1
    }
  }
}

B4_predictions$strain <- factor(B4_predictions$strain, levels = unique(SOK_data[,"strain"]))
B4_predictions$tree_sp <- factor(B4_predictions$tree_sp, levels = unique(SOK_data[,"tree_sp"]))


SOK_data_grouped_isolate <- SOK_data %>%
  group_by(strain,density,tree_sp) %>%
  summarise(total_virus=sum(total_virus) + .5 * (sum(total_virus)==0),
            total_n=sum(total_n),
            dose_response=logit(sum(total_virus)/sum(total_n)),
            dose_response_lower=logit(binom.confint(total_virus,total_n,method="wilson")$lower),
            dose_response_upper=logit(binom.confint(total_virus,total_n,method="wilson")$upper),
            dose_var=mean(ob_count))

ggplot(SOK_data_grouped_isolate) +
  geom_point(aes(x=dose_var/1000,y=dose_response,color=tree_sp)) +
  geom_errorbar(aes(x=dose_var/1000,ymin=dose_response_lower,ymax=dose_response_upper,
                    color=tree_sp),width=.30) +
  geom_line(data=B4_predictions,aes(x=dose/1000,y=prediction,color=tree_sp)) +
  facet_wrap(~strain,nrow=2,scales="free_x") +
  geom_hline(yintercept=0, linetype="dashed",size=.3) +
  scale_x_continuous(limits=c(0,6.5),expand = c(0,0)) +
  scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6),limits=c(-6.5,6),expand = c(0,0)) +
  xlab("Dose (thousands of occlusion bodies)") + ylab("logit (Proportion virus-killed)") +
  scale_color_discrete(name = "Tree",labels = c("Grand fir", "Douglas fir")) +
  theme(text=element_text(size=12,family="Palatino"),
        strip.background = element_blank(),
        panel.spacing= unit(1.5, "lines"),
        plot.margin = margin(0,20,0,10))




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
sigma_alpha_SNPV <- c(sapply(1:4,function(i) morphotype_and_tree@sim$samples[[1]]$`sigma_alpha[1]`))
sigma_alpha_MNPV <- c(sapply(1:4,function(i) morphotype_and_tree@sim$samples[[1]]$`sigma_alpha[2]`))
sigma_beta_SNPV <- c(sapply(1:4,function(i) morphotype_and_tree@sim$samples[[1]]$`sigma_beta[1]`))
sigma_beta_MNPV <- c(sapply(1:4,function(i) morphotype_and_tree@sim$samples[[1]]$`sigma_beta[2]`))

ks.test(sigma_alpha_SNPV,sigma_alpha_MNPV,alternative="greater")
ks.test(sigma_beta_SNPV,sigma_beta_MNPV,alternative="greater")

t.test(sigma_alpha_SNPV,sigma_alpha_MNPV)
t.test(sigma_beta_SNPV,sigma_beta_MNPV)

