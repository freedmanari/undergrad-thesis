library(tidyverse)
setwd("~/Desktop/school/dwyer_research/SOK")
SOK_data <- read.csv(file="SOK_reordered.csv",header=TRUE,sep=",")

dose_var <- SOK_data$ob_count


##cumulating death curves

cumulate_cols <- function(arr) {
  for (i in 2:ncol(arr)) {
    arr[,i] <- arr[,i] + arr[,i-1]
  }
  return(arr)
}

cumulated_SOK <- SOK_data
cumulated_SOK[,9:31] <- cumulate_cols(SOK_data[,9:31])



## SOK histograms

days <- colnames(SOK_data)[9:31]

day_num <- function(str) {
  return(substring(str,5,) %>% strtoi())
}

tidy <- SOK_data %>%
  pivot_longer(days, names_to="day")

tidy <- tidy %>% cbind(numeric_day = unlist(lapply(tidy$day, day_num)))

tidy %>% ggplot() +
  aes(x=numeric_day, y=value) + geom_col() +
  theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(capsid~tree_sp)




##average SOK

weighted_SOK <- SOK_data
for (i in 9:31) {
  weighted_SOK[i] <- (i - 3) * SOK_data[i]
}

weighted_SOK$weighted_sums <- rowSums(weighted_SOK[9:31])
grouped_SOK <- weighted_SOK %>%
  group_by(capsid,tree_sp,numeric_dose) %>%
  summarise(SOK=sum(weighted_sums)/sum(total_virus))

ggplot(data = grouped_SOK) + facet_wrap(~capsid) +
  aes(x=numeric_dose,y=SOK, color=tree_sp) + geom_line() + ylim(c(0,17))

ggplot(data = SOK_data) + facet_wrap(~strain,nrow=3) +
  aes(x=numeric_dose,y=SOK, color=tree_sp) + geom_line()



## glms for all data



null <- glm(dose_response ~ 1,
            family = "binomial", data = SOK_data, weights=total_n)

summary(glm(dose_response ~ dose_var, 
            family = "binomial", data = SOK_data, weights=total_n))$aic
summary(glm(dose_response ~ capsid,
            family = "binomial", data = SOK_data, weights=total_n))$aic
summary(glm(dose_response ~ tree_sp,
            family = "binomial", data = SOK_data, weights=total_n))$aic

summary(glm(dose_response ~ dose_var + capsid, 
            family = "binomial", data = SOK_data, weights=total_n))$aic
summary(glm(dose_response ~ dose_var + tree_sp,
            family = "binomial", data = SOK_data, weights=total_n))$aic
summary(glm(dose_response ~ capsid + tree_sp,
            family = "binomial", data = SOK_data, weights=total_n))$aic

summary(glm(dose_response ~ dose_var + capsid + tree_sp,
            family = "binomial", data = SOK_data, weights=total_n))$aic

best_model <- glm(dose_response ~ dose_var * capsid + tree_sp * capsid,
            family = "binomial", data = SOK_data, weights=total_n)

summary(glm(dose_response ~ dose_var * capsid + dose_var * tree_sp,
            family = "binomial", data = SOK_data, weights=total_n))$aic
summary(glm(dose_response ~ tree_sp * capsid + dose_var * tree_sp,
            family = "binomial", data = SOK_data, weights=total_n))$aic
summary(glm(dose_response ~ dose_var * capsid + tree_sp * capsid + dose_var * tree_sp,
            family = "binomial", data = SOK_data, weights=total_n))$aic

summary(glm(dose_response ~ dose_var * capsid * tree_sp,
            family = "binomial", data = SOK_data, weights=total_n))$aic


ggplot(data = SOK_data) + facet_wrap(~capsid) +
  aes(x=dose_var,y=dose_response, color=tree_sp) +
  geom_point() +
  geom_smooth(method="glm",method.args=c(family="quasibinomial")) +
  xlab("Dose (µL virus)") +
  ylab("Proportion virus-killed")


## glms for just MNPV

SNPVs <- filter(SOK_data, capsid == "SNPV")
MNPVs <- filter(SOK_data, capsid == "MNPV")

summary(glm(dose_response ~ 1,
            family = "binomial", data = MNPVs, weights=total_n))$aic

summary(glm(dose_response ~ ob_count, 
            family = "binomial", data = MNPVs, weights=total_n))$aic
summary(glm(dose_response ~ strain,
            family = "binomial", data = MNPVs, weights=total_n))$aic
summary(glm(dose_response ~ tree_sp,
            family = "binomial", data = MNPVs, weights=total_n))$aic

summary(glm(dose_response ~ ob_count + strain, 
            family = "binomial", data = MNPVs, weights=total_n))#$aic
summary(glm(dose_response ~ ob_count + tree_sp,
            family = "binomial", data = MNPVs, weights=total_n))#$aic
summary(glm(dose_response ~ strain + tree_sp,
            family = "binomial", data = MNPVs, weights=total_n))#$aic

summary(glm(dose_response ~ ob_count + strain + tree_sp,
            family = "binomial", data = MNPVs, weights=total_n))#$aic

ggplot(data = MNPVs) + facet_wrap(~strain) +
  aes(x=ob_count,y=dose_response, color=tree_sp) +
  geom_point() +
  geom_smooth(method="glm",method.args=c(family="quasibinomial")) +
  xlab("Dose (µL virus)") +
  ylab("Proportion virus-killed")




## more graph

#summary of tranmission rates
ggplot(data = SOK_data_w_control) + facet_wrap(~strain,nrow=3) +
  aes(x=ob_count,y=dose_response, color=tree_sp) + geom_line()




