This project aims to quantify estimation uncertainty (a part of epistemic uncertainty that is due to the use of a limited number of ground motion sets). The source codes are programmed in R.

The theoretical background is documented in the paper, C. Wu and H.V.Burton. Effects of probability model misspecification on the number of ground motions required for seismic performance assessment. Earthquake Spectra. DOI: 10.1177/87552930241262044


### The functionality of the main script

`SE_for_CodeBasedAnalysis.R` is to calculate the standard errors of the two parameters in a lognormally distributed probabilistic seismic demand model (PSDM). The parametric PSDM has the functional form:
$$P(EDP>edp|IM) = 1- \Phi\left(\frac{\log(edp/m_{EDP|IM})}{\sigma_{\log (EDP)|IM}}\right),$$
where $m_{EDP|IM}$ and $\sigma_{\log(EDP)|IM}$ are the conditional median and dispersion for the EDP of interest. 
This script is used to estimate the standard errors for these two parameters, using either the inverse-Fisher method or the robust method. 

`SE_for_RiskBasedAnalysis.R` is to calculate the standard errors of the  mean annual frequency of exceeding a given limit state, $\lambda_{ls}$, calculated through integrating a mean hazard curve and a lognormally distributed fragility function:
$$\lambda_{ls} &= \int_{im}^{} P(LS >ls|IM=im)|d\lambda(im)| \approx \sum_{i}^{} P(LS>ls|IM=im) |\lambda(im_{i+1}) - \lambda(im_{i})|.$$

