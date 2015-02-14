#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# bern.R
#
# (c) 2013 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2013-03-19
# last mod 2014-09-18 08:53 DW
#

rm(list=ls())

# necessary libs
library(rjags)
library(statmod)
load.module("wald")

# model
mf <- textConnection("model {
  # priors
  lambda <- 1

  alpha ~ dunif(0.7,1.2) 
  tau ~ dunif(0.7,1.2)
  gamma ~ dunif(0.7,1.2)

  for (i in 1:N) {
    x[i] ~ dwald_gamma(alpha, tau, gamma)
  }
}")

# inits  
inits1 <- list(alpha=0.7, tau=0.7, gamma=0.9)
inits2 <- list(alpha=1.0, tau=0.9, gamma=1.0)
inits3 <- list(alpha=1.1, tau=1.1, gamma=1.1)
inits <- list(inits1,inits2,inits3)
  
# Parameters to be monitored  
params <- c("alpha", "tau", "gamma")

# Genarate data
# (1) Generate drift rates
N <- 100
d <- 1
v <- 1
nu <- rgamma(N, d, v)
# (2) Generate RTs
alpha=lambda=1
RT=rinvgauss(N, mean=alpha/nu, shape=lambda*alpha^2)

dat <- list(x=RT, N=N)

# sample
j.model <- jags.model(mf, dat, inits, n.chains=3, n.adapt=100)
j.samples <- coda.samples(j.model, params, n.iter=400, thin=3)

# plot
par(mfrow=c(3,4))
plot(j.samples)

