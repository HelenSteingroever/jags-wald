#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# bern.R
#
# (c) 2013 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2013-03-19
# last mod 2015-02-14 18:16 DW
#

rm(list=ls())

# necessary libs
library(rjags)
library(statmod)
load.module("wald")

# model
my.model <- function() {
  # fixed parameters
  tau <- .5
  kappa <- 9

  # priors
  alpha ~ dunif(0.01,2.1) 

  for (i in 1:N) {
    RT[i] ~ dwald_gamma(alpha, tau, kappa)
  }
}

# inits  
inits1 <- list(alpha=1.2)
inits2 <- list(alpha=1.2)
inits3 <- list(alpha=1.2)
inits <- list(inits1,inits2,inits3)
  
# Parameters to be monitored  
params <- c("alpha")

# Genarate data
N <- 10000   # Number of reaction times

# (1) Generate drift rates
kappa <- 9   # shape: small values: right-skewed curve.
             # hihger: symmetric, shifted to the right
tau <- .5    # scale (strechting; the higher, the more stretching)
nu <- rgamma(N, shape=kappa, scale=tau)
plot(density(nu))

# (2) Generate RTs
alpha <- 1.2
RT <- rinvgauss(N, mean=alpha/nu, shape=alpha^2)
plot(density(RT))

dat <- list(RT=RT, N=N)

# sample
jags(dat, inits=NULL, params,  # inits=myinits
	 			  model.file =my.model, n.chains=3, n.iter=2000, 
          n.burnin=1000, n.thin=1, DIC=T)  


j.model <- jags.model(my.model, dat, inits=NULL, n.chains=3, n.adapt=100)
j.samples <- coda.samples(j.model, params, n.iter=400, thin=3)

# plot
par(mfrow=c(3,4))
plot(j.samples)

