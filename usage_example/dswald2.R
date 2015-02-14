#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# dswald2.R
#
# (c) 2015 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2015-02-14
# last mod 2015-02-14 17:48 DW
#

# necessary libs
library(rjags)
load.module("wald")
library(statmod)

# model
mf <- textConnection("model {
  # fixed parameters
  theta <- 0

  # priors
  alpha ~ dunif(0,2) 
  nu ~ dunif(0,2)

  for (i in 1:N) {
    x[i] ~ dswald(alpha, nu, theta)
  }
}")

inits1 <- list(alpha=0.7, nu=0.8)
inits2 <- list(alpha=0.9, nu=1.2)
inits3 <- list(alpha=1.1, nu=1.1)
inits <- list(inits1,inits2,inits3)

params <- c("alpha", "nu")

# Genarate data
N <- 1000

alpha <- 1
lambda <- 1
nu <- 1

RT <- rinvgauss(N, mean=alpha/nu, shape=lambda*alpha^2)

dat <- list(x=RT, N=N)

# sample
j.model <- jags.model(mf, dat, inits, n.chains=3, n.adapt=100)
j.samples <- coda.samples(j.model, params, n.iter=400, thin=3)

# plot
plot(j.samples)
