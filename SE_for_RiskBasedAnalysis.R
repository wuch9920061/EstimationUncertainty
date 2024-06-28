rm(list = ls())
library(VGAM)
library(sandwich)
library(lmtest)
library(readxl)
library(xlsx)
library(glue)
#library(tidyverse)
source("fn_RiskCalu.R")

# ---- Start User Input  -------#
filepath = "ExceedanceCount.xlsx"
link_fit = 'probit'
sheetname = "Site1"
# ---- End of User Input ---------# 

# When using probit, convert the moments of beta0 and beta1 to lognormal parameters (median and eta)
mean_betaln = function(mean_beta1) {1/mean_beta1}
var_betaln = function(mean_beta1, var_beta1) (1/(mean_beta1^4)) * var_beta1
mean_thetaln = function(mean_beta0, mean_beta1) exp(-mean_beta0/mean_beta1)
var_thetaln = function(mean_beta0,mean_beta1,vcovbeta01) t(c(-1/mean_beta1*exp(-mean_beta0/mean_beta1), mean_beta0/mean_beta1^2*exp(-mean_beta0/mean_beta1))) %*% vcovbeta01 %*% c(-1/mean_beta1*exp(-mean_beta0/mean_beta1), mean_beta0/mean_beta1^2*exp(-mean_beta0/mean_beta1))

# Read sheets name from xlsx file
excelpath = filepath
data = read_excel(excelpath, sheet=sheetname)
logIM = log(data$IM)
trial = cbind(data$Exceedance, data$`Non-exceedance`)


# Using classical Fisher-based estimator
glm.mod = glm(formula = trial ~ logIM, family = binomial(link=link_fit))
glm.estimate = coeftest(glm.mod)[,1]

beta_estimate = glm.estimate #estimates for GLM coefficients
glm.vcov= vcov(glm.mod)
beta_se_if = sqrt(diag(glm.vcov)) # SE estimates for GLM coefficients (Inverse Fisher)
beta_cov_if = glm.vcov[1,2]  

# Using Robust SE estimator
sw.vcov <-  vcovHC(glm.mod)
beta_se_sw = sqrt(diag(sw.vcov)) # SE estimates for GLM coefficients (Robust)
beta_cov_sw = sw.vcov[1,2] # covariance estimates for GLM coefficients (Robust)

# Convert the GLM coefficients into log-normal theta and sigma
sigma_estimate = sapply(X = glm.estimate[2],  FUN = mean_betaln)
theta_estimate = sapply(X = glm.estimate[1],  FUN = mean_thetaln, mean_beta1 = glm.estimate[2])
se_sigma_if = sqrt(sapply(X = glm.estimate[2],  FUN = var_betaln, var_beta1 = glm.vcov[2,2]))
se_sigma_sw = sqrt(sapply(X = glm.estimate[2],  FUN = var_betaln, var_beta1 = sw.vcov[2,2]))
se_theta_if = sqrt(sapply(X = glm.estimate[1],  FUN = var_thetaln, mean_beta1 = glm.estimate[2], vcovbeta01 = glm.vcov))
se_theta_sw = sqrt(sapply(X = glm.estimate[1],  FUN = var_thetaln, mean_beta1 = glm.estimate[2], vcovbeta01 = sw.vcov))
# The Robust SE is taken as the maximum of IF and Robust
se_theta_sw_larger = max(se_theta_sw,se_theta_if)
se_sigma_sw_larger = max(se_sigma_sw,se_sigma_if)

# Calculate the risk
IMvals = data$IM
lambdaIM = c(1/50,1/225,1/475,1/975,1/2475,1/10000)
dlambdaIM <- c( abs(diff(lambdaIM)), lambdaIM[length(lambdaIM)])
risk_if = fn_RiskCalu(IMvals, link_fit, beta_estimate, beta_se_if, beta_cov_if,dlambdaIM)
risk_mean_if <- risk_if$Risk.mean
risk_se_if <- risk_if$Risk.se
risk_CV_if <- risk_if$Risk.CV
risk_sw = fn_RiskCalu(IMvals, link_fit, beta_estimate, beta_se_sw, beta_cov_sw,dlambdaIM)
risk_mean_sw <- risk_sw$Risk.mean
risk_se_sw <- risk_sw$Risk.se
risk_CV_sw <- risk_sw$Risk.CV

if ((risk_CV_if <= .001) || (risk_CV_sw <= .001)) { # in this case the estimation is useless
  risk_CV_if = NaN
  risk_CV_sw = NaN
} else if ((risk_CV_if >= 2) || (risk_CV_sw >= 2)) { # in this case the estimation is useless
  risk_CV_if = NaN
  risk_CV_sw = NaN
}    

glue('The SE of lambda from the classical Fisher-based method is: {format(risk_CV_if, digits = 4)}.')
glue('The SE of lambda from the robust estimator is: {format(risk_CV_sw, digits = 4)}.')
