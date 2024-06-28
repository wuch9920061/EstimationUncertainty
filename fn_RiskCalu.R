fn_RiskCalu <- function (IMvals, link_fit, beta_estimate, beta_se, beta_cov,dlambdaIM) {
  
  # Inputs: 
  # ---- IMvals: IM levels of interest (not-log);
  # ---- link_fit: link function
  # ---- beta_estimate: GLM coefficients estimator; row1 for beta_0, row2 for beta_1
  # ---- beta_se: standard error for GLM coefficients, either from inverse Fisher or from the sandwich
  # ---- beta_cov: covariance for beta
  # ---- dlambdaIM: absolute differences of hazard curve (with length of number of IMvals)
  
  # Outputs:
  # ----- Mean, SE, and CV of risk exceeding a certain LS

  inverselink = switch(link_fit,
            "logit" = expression(1/(1+exp(-val))),
            "c-loglog" = expression(1-exp(-exp(val))),
            "cauchit" = expression(1/pi*atan(eta) + 1/2),
            "probit" = expression(pnorm(q = val))
)
  inverselink_diff = D(inverselink,"val")
  
  mean_CR=0; var_CR=0;
  
  #Calculate Collapse Risk - Mean
  for (IM_idx in 1:length(IMvals)) {
    IM_c = IMvals[IM_idx]
    eta_estimate_perIM = beta_estimate[1] + beta_estimate[2]*log(IM_c)
    val = eta_estimate_perIM
    mean_CR = mean_CR + eval(inverselink) * dlambdaIM[IM_idx]
  }
  
  #Calculate Collapse Risk - Variance
  for (IM_idxi in 1:length(IMvals)) {
    dlambdaIM_idxi = dlambdaIM[IM_idxi]
    IM_i = IMvals[IM_idxi]
    eta_estimate_idxi = beta_estimate[1] + beta_estimate[2]*log(IM_i)
    eta_var_idxi = beta_se[1]^2 + beta_se[2]^2*log(IM_i)^2 + 2*beta_cov*log(IM_i)
    val = eta_estimate_idxi
    p_mean_idxi = eval(inverselink)
    p_var_idxi = eval(inverselink_diff)^2*eta_var_idxi
    if (is.nan(p_var_idxi)) {p_var_idxi=0}
    if (p_var_idxi >=0) {
      p_se_idxi = sqrt(p_var_idxi)
    } else {
      p_se_idxi = 1e-10
    }
    
    
    for (IM_idxj in 1:length(IMvals)) {
      dlambdaIM_idxj = dlambdaIM[IM_idxj]
      IM_j = IMvals[IM_idxj]
      eta_estimate_idxj = beta_estimate[1] + beta_estimate[2]*log(IM_j)
      eta_var_idxj = beta_se[1]^2 + beta_se[2]^2*log(IM_j)^2 + 2*beta_cov*log(IM_j)
      val = eta_estimate_idxj
      p_mean_idxj = eval(inverselink)
      p_var_idxj = eval(inverselink_diff)^2*eta_var_idxj
      
      if (is.nan(p_var_idxj)) {p_var_idxj=0}
      
      if (p_var_idxj >=0) {
      p_se_idxj = sqrt(p_var_idxj)
      } else {
        p_se_idxj = 1e-10
      }
      
      var_CR = var_CR + dlambdaIM_idxi*dlambdaIM_idxj*p_se_idxi*p_se_idxj
    }
  }
  
  
  
  return(list(Risk.mean = mean_CR, Risk.se = sqrt(var_CR), Risk.CV = sqrt(var_CR)/mean_CR))
}