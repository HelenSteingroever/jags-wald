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
mf <- textConnection("model {
  # fixed parameters
  tau <- 1
  kappa <- 1

  # priors
  alpha ~ dunif(0.1,2.1) 

  for (i in 1:N) {
    x[i] ~ dwald_gamma(alpha, tau, kappa)
  }
}")

# inits  
inits1 <- list(alpha=1.0)
inits2 <- list(alpha=1.0)
inits3 <- list(alpha=1.0)
inits <- list(inits1,inits2,inits3)
  
# Parameters to be monitored  
params <- c("alpha")

# Genarate data
# (1) Generate drift rates
N <- 500
tau <- 1
kappa <- 1
nu <- rgamma(N, tau, kappa)
# (2) Generate RTs
alpha=1
RT=rinvgauss(N, mean=alpha/nu, shape=alpha^2)

dat <- list(x=RT, N=N)

# sample
j.model <- jags.model(mf, dat, inits, n.chains=3, n.adapt=100)
j.samples <- coda.samples(j.model, params, n.iter=400, thin=3)

# plot
par(mfrow=c(3,4))
plot(j.samples)

