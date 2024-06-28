library(MASS)

loglik <- function(parvec,x) (-sum(dnorm((x), mean=parvec[1], sd=parvec[2], log = T)))

# Now we calculate the bread and meat
logl <- expression(-0.5*log(2*pi) -0.5* log(sigma^2) - 0.5*(x - mu)^2/sigma^2)
#score.mu <- D(logl,'mu')
score.mu <- expression(0.5 * (2 * (x - mu))/sigma^2)
#score.sigma <- D(logl,'sigma')
score.sigma <- expression(-(0.5 * (2 * sigma/sigma^2) - 0.5 * (x - mu)^2 * (2 * sigma)/(sigma^2)^2))
#hessian.mu.mu <- D(score.mu , 'mu')
hessian.mu.mu <- expression(-(0.5 * 2/sigma^2) +x-x)
#hessian.mu.sigma <- D(score.mu , 'sigma')
hessian.mu.sigma <- expression(-(0.5 * (2 * (x - mu)) * (2 * sigma)/(sigma^2)^2))
#hessian.sigma.sigma <- D(score.sigma , 'sigma')
hessian.sigma.sigma <- expression(-(0.5 * (2/sigma^2 - 2 * sigma * (2 * sigma)/(sigma^2)^2) - (0.5 * 
                                                                                                 (x - mu)^2 * 2/(sigma^2)^2 - 0.5 * (x - mu)^2 * (2 * sigma) * 
                                                                                                 (2 * (2 * sigma * (sigma^2)))/((sigma^2)^2)^2)))

# calculate loglik values
loglfun <- function(theta,x) {
  mu <- theta[1]
  sigma <- theta[2]
  sum(eval(logl))
}

gradfun <- function(theta,x) {
  mu <- theta[1]
  sigma <- theta[2]
  grad.mu <- sum(eval(score.mu))
  grad.sigma <- sum(eval(score.sigma))
  c(grad.mu,grad.sigma)
}

vfun <- function(theta, x) {
  mu <- theta[1]
  sigma <- theta[2]
  s.mu <- eval(score.mu)
  s.sigma <- eval(score.sigma)
  v.mu.mu <- sum(s.mu^2)
  v.mu.sigma <- sum(s.mu * s.sigma)
  v.sigma.sigma <- sum(s.sigma^2)
  matrix(c(v.mu.mu, v.mu.sigma, v.mu.sigma, v.sigma.sigma), 2, 2)
}

jfun <- function(theta, x) {
  mu <- theta[1]
  sigma <- theta[2]
  j.mu.mu <- sum(- eval(hessian.mu.mu))
  j.mu.sigma <- sum(- eval(hessian.mu.sigma))
  j.sigma.sigma <- sum(- eval(hessian.sigma.sigma))
  matrix(c(j.mu.mu, j.mu.sigma, j.mu.sigma, j.sigma.sigma), 2, 2)
}
