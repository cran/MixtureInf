## Generate a random sample from a mixture of binomial distributions
##
## n:      sample size
## alpha:  component mixing proportions
## prob:   vector of probability of success of all components
## size:   size parameter for Binomial mixture.
rbinommix <- function(n, alpha, prob, size) {
  m <- length(alpha)                    # number of components
  
  z <- sample(m, n, replace=TRUE, prob=alpha)
  rbinom(n, size, prob[z])
}


## Generate a random sample from a mixture of exponential distributions
## 
## n:      sample size
## alpha:  component mixing proportions
## scale:  vector of scales (inverse rates) of all components
rexpmix <- function (n, alpha, scale) {
  m <- length(alpha)                    # number of components
  
  z <- sample(m, n, replace=TRUE, prob=alpha)
  rexp(n, 1/scale[z])
}


## This function computes the cdf of Exponential mixture minus a given probability alp0.
## 
## x:      data.
## alpha:  mixing proportions.
## theta:  mixing parameters.
## alp0:   a given probability.
pexpmix <- function(x, alpha, theta, alp0) {
  sum( alpha * pexp(x, rate=1/theta) )-alp0
}


## This function computes the alp0-quantile of Exponential mixture.
## 
## alpha:  mixing proportions.
## theta:  mixing parameters.
## alo0:   a given probability.
qexpmix <- function(alpha, theta, alp0) {
  uniroot(pexpmix, c(0, qexp(alp0, rate=1/max(theta))), alpha=alpha, theta=theta, alp0=alp0)$root
}


## Generate a random sample from a mixture of normal distributions
##
## n:      sample size
## alpha:  component mixing proportions
## mean:   vector of means of all components
## sigma:  vector of standard deviations of each component
rnormmix <- function (n, alpha, mean, sd=rep(1, length(alpha))) {
  m <- length(alpha)                    # number of components
  
  z <- sample(m, n, replace=TRUE, prob=alpha)
  rnorm(n, mean[z], sd[z])
}


## This function computes the cdf of Normal mixture minus a given probability alp0.
##
## x:      data.
## alpha:  mixing proportions.
## theta:  mixing parameters.
## alp0:   a given probability.
pnormmix <- function(x, alpha, theta, alp0) {	
  sum(alpha*pnorm(x, theta, 1)) - alp0
}


## This function computes the alp0-quantile of Normal mixture.
##
## alpha:  mixing proportions.
## theta:  mixing parameters.
## alp0:   a given probability.
qnormmix <- function(alpha, theta, alp0) {	
  uniroot(pnormmix, c(qnorm(alp0, min(theta)), qnorm(alp0, max(theta))), alpha=alpha, theta=theta, alp0=alp0)$root
}


## Generate a random sample from a mixture of Poisson distributions
## 
## n:      sample size
## alpha:  component mixing proportions
## lambda: vector of means of all components
rpoismix <- function(n, alpha, lambda) {
  m <- length(alpha)                    # number of components
  
  z <- sample(m, n, replace=TRUE, prob=alpha)
  rpois(n, lambda[z])
}
