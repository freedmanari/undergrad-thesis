library(tidyverse)
#setwd("~/Desktop/school/dwyer_research/SOK")
SOK_data_w_control <- read.csv(file="DoseR_SOK_reformatted.csv",header=TRUE,sep=",")
SOK_data <- SOK_data_w_control %>% filter(capsid != "none")


#tidy the data

days <- colnames(SOK_data)[9:31]

day_num <- function(str) {
  return(substring(str,5,) %>% strtoi())
}

tidy <- SOK_data %>%
  pivot_longer(days, names_to="day")

tidy <- tidy %>% cbind(numeric_day = unlist(lapply(tidy$day, day_num)))
#write.csv(tidy, "tidy_SOK.csv", row.names = FALSE)

tidy %>% ggplot() +
  aes(x=numeric_day, y=value) + geom_col() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(capsid~tree_sp)



MNPV_DO <- filter(tidy, capsid == "MNPV", tree_sp == "DO")
MNPV_GR <- filter(tidy, capsid == "MNPV", tree_sp == "GR")
SNPV_DO <- filter(tidy, capsid == "SNPV", tree_sp == "DO")
SNPV_GR <- filter(tidy, capsid == "SNPV", tree_sp == "GR")

MNPV_DO_days <- rep(MNPV_DO$numeric_day,MNPV_DO$value)
MNPV_GR_days <- rep(MNPV_GR$numeric_day,MNPV_GR$value)
SNPV_DO_days <- rep(SNPV_DO$numeric_day,SNPV_DO$value)
SNPV_GR_days <- rep(SNPV_GR$numeric_day,SNPV_GR$value)



#likelihood function

logLHoodHtg<-function(data, par){
  
  alpha = par[1];
  beta  = par[2];
  
  #data = MNPV_DO_days
  
  logLHood = dgamma(data, shape=alpha, scale=beta, log=TRUE); #2. Type in the log-likelihood function for the pure death model
  return(-sum(logLHood));  #We are looking for the negative sum of the log likelihoods across experimental units, because optim is a minimizer
  
}


#MNPV_DO

hist(MNPV_DO_days, freq=FALSE, breaks = 15)

MNPV_DO_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = MNPV_DO_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = MNPV_DO_OptOut$par[1], scale = MNPV_DO_OptOut$par[2])

lines(xnums,ynums)

MNPV_DO_LHood <- MNPV_DO_OptOut$value


#MNPV_GR

hist(MNPV_GR_days, freq=FALSE, breaks = 15)

MNPV_GR_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = MNPV_GR_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = MNPV_GR_OptOut$par[1], scale = MNPV_GR_OptOut$par[2])

lines(xnums,ynums)

MNPV_GR_LHood <- MNPV_GR_OptOut$value



#SNPV_DO

hist(SNPV_DO_days, freq=FALSE, breaks = 15)

SNPV_DO_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = SNPV_DO_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = SNPV_DO_OptOut$par[1], scale = SNPV_DO_OptOut$par[2])

lines(xnums,ynums)

SNPV_DO_LHood <- SNPV_DO_OptOut$value



#SNPV_GR

hist(SNPV_GR_days, freq=FALSE, breaks = 15)

SNPV_GR_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = SNPV_GR_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = SNPV_GR_OptOut$par[1], scale = SNPV_GR_OptOut$par[2])

lines(xnums,ynums)

SNPV_GR_LHood <- SNPV_GR_OptOut$value







##grouping by morphotype

MNPV <- filter(tidy, capsid == "MNPV")
SNPV <- filter(tidy, capsid == "SNPV")

MNPV_days <- rep(MNPV$numeric_day, MNPV$value)
SNPV_days <- rep(SNPV$numeric_day, SNPV$value)


#MNPV

hist(MNPV_days, freq=FALSE, breaks = 15)

MNPV_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = MNPV_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = MNPV_OptOut$par[1], scale = MNPV_OptOut$par[2])

lines(xnums,ynums)

MNPV_LHood <- MNPV_OptOut$value



#SNPV

hist(SNPV_days, freq=FALSE, breaks = 15)

SNPV_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = SNPV_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = SNPV_OptOut$par[1], scale = SNPV_OptOut$par[2])

lines(xnums,ynums)

SNPV_LHood <- SNPV_OptOut$value



##grouping by tree

DO <- filter(tidy, tree_sp == "DO")
GR <- filter(tidy, tree_sp == "GR")

DO_days <- rep(DO$numeric_day, DO$value)
GR_days <- rep(GR$numeric_day, GR$value)


#DO

hist(DO_days, freq=FALSE, breaks = 15)

DO_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = DO_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = DO_OptOut$par[1], scale = DO_OptOut$par[2])

lines(xnums,ynums)

DO_LHood <- DO_OptOut$value



#GR

hist(GR_days, freq=FALSE, breaks = 15)

GR_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = GR_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = GR_OptOut$par[1], scale = GR_OptOut$par[2])

lines(xnums,ynums)

GR_LHood <- GR_OptOut$value



##grouping by everything

SOK_days <- rep(tidy$numeric_day, tidy$value)

hist(SOK_days, freq=FALSE, breaks = 15)

SOK_OptOut = optim(par=c(5,0.5),fn=logLHoodHtg, data = SOK_days);

xnums <- seq(0,28,0.1)
ynums <- dgamma(xnums, shape = SOK_OptOut$par[1], scale = SOK_OptOut$par[2])

lines(xnums,ynums)

SOK_LHood <- SOK_OptOut$value



##AICs
AICCalc <- function(negLogLHood, K) {
  2*negLogLHood + 2*K
}


#tree and capsid
AICCalc(MNPV_DO_LHood + MNPV_GR_LHood + SNPV_DO_LHood + SNPV_GR_LHood, 2*4)

AICCalc(MNPV_LHood + SNPV_LHood, 2*2)

AICCalc(DO_LHood + GR_LHood, 2*2)

AICCalc(SOK_LHood, 2*1)
