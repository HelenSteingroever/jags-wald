#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# dswald_miller.R
#
# (c) 2015 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2015-02-14
# last mod 2015-02-16 19:56 DW
#

# nu = delta
# alpha = gamma
# theta = theta

# necessary libs
library(rjags)
load.module("wald")
library(RWiener)

# model
mf <- textConnection("model {
  # priors
  alpha ~ dunif(0,10) 
  nu~ dunif(0,10)
  theta ~ dunif(0,1.0)

  for (i in 1:length(rt_hit)) {
    rt_hit[i] ~ dswald(alpha, nu, theta)
  }
  for (i in 1:length(rt_miss)) {
    rt_miss[i] ~ pswald_upper(alpha, nu, theta)
  }
}")

inits1 <- list(alpha=2.7, nu=0.8, theta=0.1)
inits2 <- list(alpha=0.9, nu=1.2, theta=0.01)
inits3 <- list(alpha=1.1, nu=1.1, theta=0.05)
inits <- list(inits1,inits2,inits3)

params <- c("alpha", "nu", "theta")

# Genarate data
N <- 1000

alpha <- 1
nu <- 1
theta <- 0.3

rswald <- function(N, alpha, nu, theta) {
  beta.par <- 0.5
  x <- rwiener(N, alpha/(1-beta.par), theta, beta.par, nu)
  return(x)
}

RT <- rswald(N, alpha, nu, theta)

dat <- list(rt_hit=RT$q[RT$resp=="upper"], 
            rt_miss=RT$q[RT$resp=="lower"])

# sample
j.model <- jags.model(mf, dat, inits, n.chains=3, n.adapt=100)
j.samples <- coda.samples(j.model, params, n.iter=500, thin=3)

# plot
plot(j.samples)
