#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# bern.R
#
# (c) 2013 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2013-03-19
# last mod 2013-04-29 18:34 DW
#

# necessary libs
library(rjags)
library(truncnorm)
library(statmod)
load.module("wald")

# model
mf <- textConnection("model {
  # priors
  lambda ~ dunif(0.01,1.5)
  alpha ~ dunif(0.01,1.5) 
  v ~ dunif(0.01,1.5)
  d ~ dunif(0.01,1.5)

  for (i in 1:N) {
    x[i] ~ dwald_trunc(lambda, alpha, v, d)
  }

}")

# Genarate data
# (1) Generate drift rates
N <- 1000
nu <- rtruncnorm(N, a=0, b=Inf, mean=1, sd=1/sqrt(1))  
# (2) Generate RTs
alpha=lambda=1
RT=rinvgauss(N, mu=alpha/nu, lambda*alpha^2)
plot(density(RT[RT<=3]))

x <- RT[RT<=3]
N <- length(x)
dat <- list(x=x, N=N)

# inits
inits1 <- list(lambda=.7, alpha=.7, v=.7, d=.7)
inits2 <- list(lambda=.9, alpha=.9, v=.9, d=.9)
inits3 <- list(lambda=1.1, alpha=1.1, v=1.1, d=1.1)
inits <- list(inits1,inits2,inits3)

# sample
j.model <- jags.model(mf, dat, inits, n.chains=3, n.adapt=0)
j.samples <- coda.samples(j.model, c("alpha","lambda", "v", "d"), n.iter=2000, thin=1)

# plot
plot(j.samples[,1])
