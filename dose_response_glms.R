library(tidyverse)
library(AICcmodavg)
library(boot)
library(binom)
library(gmodels)
library(ggstance)
library(scales)
library(ggpubr)
theme_set(theme_pubr())
SOK_data <- read.csv(file="SOK_reordered.csv",header=TRUE,sep=",")
SOK_data_MNPV_first <- SOK_data
SOK_data$strain <- factor(SOK_data$strain,levels=c("COL","CUB","LOV","LST","DRY","KLP","TAM","TMB"))
SOK_data$capsid <- factor(SOK_data$capsid,levels=c("SNPV","MNPV"))
SOK_data$tree_sp <- factor(SOK_data$tree_sp,levels=c("GR","DO"))

dose_var <- SOK_data$ob_count


# ##cumulating death curves
# 
# cumulate_cols <- function(arr) {
#   for (i in 2:ncol(arr)) {
#     arr[,i] <- arr[,i] + arr[,i-1]
#   }
#   return(arr)
# }
# 
# cumulated_SOK <- SOK_data
# cumulated_SOK[,9:31] <- cumulate_cols(SOK_data[,9:31])

# ## SOK histograms
# 
# days <- colnames(SOK_data)[9:31]
# 
# day_num <- function(str) {
#   return(substring(str,5,) %>% strtoi())
# }
# 
# tidy <- SOK_data %>%
#   pivot_longer(days, names_to="day")
# 
# tidy <- tidy %>% cbind(numeric_day = unlist(lapply(tidy$day, day_num)))
# 
# tidy %>% ggplot() +
#   aes(x=numeric_day, y=value) + geom_col() +
#   theme(axis.text.x = element_text(angle = 90)) +
#   facet_grid(capsid~tree_sp)




# ##average SOK
# 
# weighted_SOK <- SOK_data
# for (i in 9:31) {
#   weighted_SOK[i] <- (i - 3) * SOK_data[i]
# }
# 
# weighted_SOK$weighted_sums <- rowSums(weighted_SOK[9:31])
# grouped_SOK <- weighted_SOK %>%
#   group_by(capsid,tree_sp,numeric_dose) %>%
#   summarise(SOK=sum(weighted_sums)/sum(total_virus))
# 
# ggplot(data = grouped_SOK) + facet_wrap(~capsid) +
#   aes(x=numeric_dose,y=SOK, color=tree_sp) + geom_line() + ylim(c(0,17))
# 
# ggplot(data = SOK_data) + facet_wrap(~strain,nrow=2) +
#   aes(x=numeric_dose,y=SOK, color=tree_sp) + geom_line()



## plot of average mortality rate for morphotype-tree combos
grouped <- SOK_data %>%
  group_by(capsid, tree_sp) %>%
  summarise(y = sum(total_virus)/sum(total_n),
            ymin=binom.confint(sum(total_virus),sum(total_n),method="wilson")$lower,
            ymax=binom.confint(sum(total_virus),sum(total_n),method="wilson")$upper) %>% 
  as.data.frame()

ggplot(data = grouped) +
  geom_line(aes(x=capsid,y=y,color=tree_sp,group=tree_sp)) +
  geom_point(aes(x=capsid,y=y,color=tree_sp)) +
  geom_errorbar(width=.25,
                aes(x=capsid,ymin=ymin,ymax=ymax,color=tree_sp))+
  scale_y_continuous(limits = c(0,.8),expand = c(0,0)) +
  scale_color_discrete(name = "Tree",labels = c("GF", "DF")) +
  xlab("Morphotype") +
  ylab("Proportion virus-killed") +
  theme(text=element_text(size=12,family="Palatino"))



## glms for all data



A1 <- glm(dose_response ~ 1,
            family = "binomial", data = SOK_data, weights=total_n)

A2 <- glm(dose_response ~ dose_var, 
            family = "binomial", data = SOK_data, weights=total_n)
A3 <- glm(dose_response ~ capsid,
            family = "binomial", data = SOK_data, weights=total_n)
A4 <- glm(dose_response ~ tree_sp,
            family = "binomial", data = SOK_data, weights=total_n)

A5 <- glm(dose_response ~ dose_var + capsid, 
            family = "binomial", data = SOK_data, weights=total_n)
A6 <- glm(dose_response ~ dose_var + tree_sp,
            family = "binomial", data = SOK_data, weights=total_n)
A7 <- glm(dose_response ~ capsid + tree_sp,
            family = "binomial", data = SOK_data, weights=total_n)

A8 <- glm(dose_response ~ dose_var + capsid + tree_sp,
            family = "binomial", data = SOK_data, weights=total_n)

A9 <- glm(dose_response ~ dose_var * capsid,
                  family = "binomial", data = SOK_data, weights=total_n)
A10 <- glm(dose_response ~ dose_var * tree_sp,
          family = "binomial", data = SOK_data, weights=total_n)
A11 <- glm(dose_response ~ capsid * tree_sp,
          family = "binomial", data = SOK_data, weights=total_n)

A12 <- glm(dose_response ~ dose_var * capsid + dose_var * tree_sp,
            family = "binomial", data = SOK_data, weights=total_n)
A13 <- glm(dose_response ~ dose_var * capsid + tree_sp * capsid,
            family = "binomial", data = SOK_data, weights=total_n) # best_model
A14 <- glm(dose_response ~ capsid * tree_sp + dose_var * tree_sp,
            family = "binomial", data = SOK_data, weights=total_n)

A15 <- glm(dose_response ~ dose_var * capsid + dose_var * tree_sp + capsid * tree_sp,
            family = "binomial", data = SOK_data, weights=total_n)

A16 <- glm(dose_response ~ dose_var * capsid * tree_sp,
            family = "binomial", data = SOK_data, weights=total_n)


glms <- list(A1=A1,A2=A2,A3=A3,A4=A4,
             A5=A5,A6=A6,A7=A7,A8=A8,
             A9=A9,A10=A10,A11=A11,A12=A12,
             A13=A13,A14=A14,A15=A15,A16=A16)
formulas <- sapply(glms, function(g) str_remove_all(str_replace_all(str_replace_all(str_replace_all(
  as.character(formula(g))[3],"tree_sp","T"),"dose_var","D"),"capsid","M")," "))
aicc_table <- aictab(glms,second.ord=TRUE)
aicc_table <- cbind(formula=formulas[aicc_table$Modnames %>% substring(,first=2) %>% strtoi],aicc_table)
rownames(aicc_table) <- c()


## graph of best model
#alpha+beta1*D+beta2*M+beta3*T+beta4*D*M+beta5*T*M
doses <- seq(0,6000,1)
A13_predictions <- data.frame(dose=numeric(4*length(doses)),
                              capsid=character(4*length(doses)),
                              tree_sp=character(4*length(doses)),
                              prediction=numeric(4*length(doses)),
                              stringsAsFactors = FALSE)
i <- 1
for (c in c("SNPV","MNPV")) {
  c_int <- c == "MNPV"
  for (t in c("GR","DO")) {
    t_int <- t == "DO"
    for (d in doses) {
      A13_predictions[i,] <- list(dose=d,
                                  capsid=c,
                                  tree_sp=t,
                                  prediction=inv.logit(coef(A13)[["(Intercept)"]] +
                                                       coef(A13)[["dose_var"]]*d +
                                                       coef(A13)[["capsidMNPV"]]*c_int +
                                                       coef(A13)[["tree_spDO"]]*t_int +
                                                       coef(A13)[["dose_var:capsidMNPV"]]*d*c_int +
                                                       coef(A13)[["capsidMNPV:tree_spDO"]]*t_int*c_int))
      i <- i+1
    }
  }
}

A13_predictions$capsid <- factor(A13_predictions$capsid, levels = c("SNPV","MNPV"))
A13_predictions$tree_sp <- factor(A13_predictions$tree_sp, levels = c("GR","DO"))

ggplot() +
  geom_line(data = A13_predictions, aes(x=dose, y=prediction, color=tree_sp)) +
  geom_point(data = SOK_data, aes(x=dose_var,y=dose_response, color=tree_sp)) +
  facet_wrap(~capsid) + xlab("Dose (µL virus)") + ylab("Proportion virus-killed")

SOK_data_grouped <- SOK_data %>%
                    group_by(capsid,density,tree_sp) %>%
                    summarise(total_virus=sum(total_virus),
                              total_n=sum(total_n),
                              dose_response=sum(total_virus)/sum(total_n),
                              dose_response_lower=binom.confint(total_virus,total_n,method="wilson")$lower,
                              dose_response_upper=binom.confint(total_virus,total_n,method="wilson")$upper,
                              dose_var=mean(ob_count),
                              dose_lower=mean_cl_normal(ob_count)$ymin,
                              dose_upper=mean_cl_normal(ob_count)$ymax)


ggplot() +
  geom_line(data = A13_predictions, aes(x=dose/1000, y=prediction, color=tree_sp)) +
  geom_point(data = SOK_data_grouped, aes(x=dose_var/1000,y=dose_response, color=tree_sp)) +
  geom_errorbar(data=SOK_data_grouped,
                width=.15,
                aes(x=dose_var/1000,ymin=dose_response_lower,ymax=dose_response_upper,color=tree_sp))+
  geom_errorbarh(data=SOK_data_grouped,
                 width=.02,aes(y=dose_response,xmin=dose_lower/1000,xmax=dose_upper/1000,color=tree_sp))+
  facet_wrap(~capsid) + xlab("Dose (thousands of occlusion bodies)") + ylab("Proportion virus-killed") +
  scale_x_continuous(limits = c(0,6),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1),expand = c(0,0)) +
  scale_color_discrete(name = "Tree",labels = c("Grand fir", "Douglas fir")) +
  theme(text=element_text(size=12,family="Palatino"),
        strip.background = element_blank(),
        panel.spacing= unit(1.5, "lines"),
        plot.margin = margin(0,20,0,10),
        legend.margin=margin(0,0,-10,0))



## more graph

ggplot(data = SOK_data) + facet_wrap(~capsid) +
  aes(x=dose_var,y=dose_response, color=tree_sp) +
  geom_point() +
  geom_smooth(method="glm",method.args=c(family="quasibinomial")) +
  xlab("Dose (µL virus)") +
  ylab("Proportion virus-killed")

#summary of tranmission rates
SOK_data_grouped_isolate <- SOK_data %>%
  group_by(strain,density,tree_sp) %>%
  summarise(total_virus=sum(total_virus),
            total_n=sum(total_n),
            dose_response=sum(total_virus)/sum(total_n),
            dose_response_lower=binom.confint(total_virus,total_n,method="wilson")$lower,
            dose_response_upper=binom.confint(total_virus,total_n,method="wilson")$upper,
            dose_var=mean(ob_count))


ggplot(data = SOK_data_grouped_isolate) +
  facet_wrap(~strain,nrow=2,scales="free_x") +
  aes(x=dose_var/1000,y=dose_response, color=tree_sp) +
  geom_errorbar(width=.3,
                aes(x=dose_var/1000,ymin=dose_response_lower,ymax=dose_response_upper,color=tree_sp))+
  geom_line() + geom_point() +
  scale_color_discrete(name = "Tree",labels = c("Grand fir", "Douglas fir")) +
  scale_x_continuous(limits=c(0,6),expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1),expand = c(0,0)) +
  xlab("Dose (thousands of occlusion bodies)") + ylab("Proportion virus-killed") +
  theme(text=element_text(size=12,family="Palatino"),
        strip.background = element_blank(),
        panel.spacing= unit(1.5, "lines"),
        plot.margin = margin(0,20,0,10))



