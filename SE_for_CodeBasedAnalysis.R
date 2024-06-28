library(readxl)
library(tidyverse)
library(xlsx)
library(glue)
source("fn_fisher_robust_handcalc.R")


filepath = "ColDuc.xlsx"
num_GM = 55 # How many GM were you using for NTHR


data <- read_excel(filepath)
sitename <- 'lat37.8lon-122.25'
Colduc <- data[data[,1] == sitename,]
Colduc <- unname(unlist(Colduc %>% dplyr::select(-1)))

# Retrieve EDP values from the 1st ground motion until the specified num_GM
Colduc_pernumGM = Colduc[1:num_GM]

# First we estimate the median and log-std from their unbiased estimator
meanlog_hat_mom = mean(log(Colduc_pernumGM));
median_hat_mom = exp(mean(log(Colduc_pernumGM)));
sigma_hat_mom = sqrt(sum((log(Colduc_pernumGM)-mean(log(Colduc_pernumGM)))^2)/length(log(Colduc_pernumGM)));
glue('The median and log-std calculated estimated using their unbiased estimator are: {format(median_hat_mom, digits = 4)} and {format(sigma_hat_mom,digits = 4)}.')

# Using the Robust SE to estimate the SE of median and log-std
Vhat <- vfun(c(log(median_hat_mom),sigma_hat_mom), x=log(Colduc_pernumGM))
Jhat <- jfun(c(log(median_hat_mom),sigma_hat_mom), x=log(Colduc_pernumGM))
Mhat <- solve(Jhat) %*% Vhat %*% solve(Jhat) # the sandwich asymptotic variance
sd_sigma_SW = sqrt(Mhat[2,2])
sd_meanlog_SW = sqrt(Mhat[1,1])
sd_median_SW <- sqrt((exp(meanlog_hat_mom))^2 * sd_meanlog_SW^2)

IF <- ginv(Jhat)
sd_meanlog_IF = sqrt(IF[1,1])
sd_sigma_IF = sqrt(IF[2,2])
sd_median_IF <- sqrt((exp(meanlog_hat_mom))^2 * sd_meanlog_IF^2)

sd_median_larger = max(sd_median_IF,sd_median_SW)
sd_sigma_larger = max(sd_sigma_IF,sd_sigma_SW)

glue('The SE for median and log-std calculated from the classical Fisher-based method are: {format(sd_median_IF, digits = 4)} and {format(sd_sigma_IF,digits = 4)}.')
glue('The SE for median and log-std calculated from the robust estimator are: {format(sd_median_larger, digits = 4)} and {format(sd_sigma_larger,digits = 4)}.')
