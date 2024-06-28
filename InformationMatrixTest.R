rm(list = ls())
library(VGAM)
library(sandwich)
library(lmtest)
library(readxl)
library(tidyverse)
library(RobustSE)

#--------- Start User Input ------ #
filepath = "ExceedanceCount.xlsx"
link_fit = 'probit'
sheetname = "Site1"
B = 300 # The first bootstrap number in the nested simulation
B2 = 200# The second bootstrap number in the nested simulation
#-------- End of User Input-------- #

data_perstruct = read_excel(filepath, sheetname)
logIM = log(data_perstruct$IM)
trial = cbind(data_perstruct$Exceedance, data_perstruct$`Non-exceedance`)
# convert the two-column integer matrix formula into factor ( see https://stat.ethz.ch/R-manual/R-devel/library/stats/html/family.html)
trial_new = c()
logIM_new = c()
for (row in 1:length(logIM)) {
  logIM_new = c(logIM_new, rep(logIM[row], trial[row,1] + trial[row,2]))
  trial_new = c(trial_new, rep(1, trial[row,1]), rep(0, trial[row,2]))
}

# Fit the generalized linear model
glm.mod = glm(formula = trial_new ~ logIM_new, family = binomial(link=link_fit))
#a_full = GIM(glm.mod, full = TRUE, B = 50, B2 = 50)
#p_val = (a_full$`GIM pval`)

# The Information Matrix test

#--- Calculate the V (Meat) and I (Fisher Info) matrices
grad = sandwich::estfun(glm.mod) # the score function
V_mat = t(grad) %*% grad
I_mat = MASS::ginv(vcov(glm.mod)) # the fisher info matrix is the inverse of the parameter variance matrix

glm.mod.list = model.matrix(glm.mod)
glm.mod.coefficients = glm.mod$coefficients
sample.size = nrow(glm.mod.list)

D_hat = sample.size^(-1/2)*(V_mat - I_mat)
d_hat = diag(D_hat)

T_list= list()
d_hat_b_list = list()
d_b_avg = rep(0,length(d_hat))
p_predicted = pnorm(glm.mod.list %*% glm.mod.coefficients)

for (b in 1:B){
  yb = rbinom(sample.size, 1, prob=p_predicted)
  glm.mod.b = glm(yb ~ logIM_new, family=binomial("probit"))
  grad.b = sandwich::estfun(glm.mod.b) # the score function
  V_mat.b = t(grad.b) %*% grad.b
  I_mat.b = MASS::ginv(vcov(glm.mod.b))
  d_hat.b = diag(sample.size^(-1/2) * (V_mat.b - I_mat.b))
  d_hat_b_list[[b]] = d_hat.b
  d_b_avg = d_b_avg + d_hat.b
  
  d_hat_b2_list = list()
  d_b2_avg = rep(0,length(d_hat))
  
  for (b2 in 1:B2){
    yb2 =  rbinom(sample.size, 1, prob=p_predicted)
    glm.mod.b2 = glm(yb2 ~ logIM_new, family=binomial("probit"))
    grad.b2 = sandwich::estfun(glm.mod.b2) # the score function
    V_mat.b2 = t(grad.b2) %*% grad.b2
    I_mat.b2 = MASS::ginv(vcov(glm.mod.b2))
    d_hat.b2 = diag(sample.size^(-1/2) * (V_mat.b2 - I_mat.b2))
    d_hat_b2_list[[b2]] = d_hat.b2
    d_b2_avg = d_b2_avg + d_hat.b2
  }# end B2
  d_b2_avg = d_b2_avg/B2
  
  V_d.b.sum = matrix(0, nrow=length(d_b2_avg), ncol=length(d_b2_avg))
  for(b2 in 1:B2){
    V_d.b.sum <- V_d.b.sum + (d_hat_b2_list[[b2]] - d_b2_avg)%*%t(d_hat_b2_list[[b2]] - d_b2_avg)
  }
  Q_hat.b = V_d.b.sum/(B2-1)
  omega_hat.b = t(d_hat.b)%*%MASS::ginv(Q_hat.b)%*%d_hat.b
  T[b] = omega_hat.b
} # end B
d_b_avg = d_b_avg/B

V_d.sum = matrix(0, nrow=length(d_b_avg), ncol=length(d_b_avg))
for(b in 1:B){
  V_d.sum <- V_d.sum + (d_hat_b_list[[b]] - d_b_avg)%*%t(d_hat_b_list[[b]] - d_b_avg)
}
Q_hat = V_d.sum/(B-1)
omega_hat = t(d_hat)%*%MASS::ginv(Q_hat)%*%d_hat
p_val = sum(T > as.numeric(omega_hat))/(B)
print(p_val)

