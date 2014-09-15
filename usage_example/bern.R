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
load.module("bernoulli")

# model
mf <- textConnection("model {
  # priors
  p ~ dunif(0,1)

  for (i in 1:N) {
    x[i] ~ dbern2(p)
  }

  y[1] <- logbern(0,.5)
  y[2] <- logbern(1,.7)

}")

# data
x <- rbinom(100,1,.7)
dat <- list(x=x, N=100)

# inits
inits1 <- list(p=.5)
inits2 <- list(p=.9)
inits3 <- list(p=.1)
inits <- list(inits1,inits2,inits3)

# sample
j.model <- jags.model(mf, dat, inits, n.chains=3, n.adapt=0)
j.samples <- coda.samples(j.model, c("p","y"), n.iter=2000, thin=1)

# plot
plot(j.samples[,1])
