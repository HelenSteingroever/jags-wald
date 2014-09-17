#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# bern.R
#
# (c) 2013 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2013-03-19
# last mod 2014-09-17 17:06 DW
#

rm(list=ls())

# necessary libs
library(rjags)
library(truncnorm)
library(statmod)
load.module("wald")

# Should lambda be fixed to 1?
lambda_fixed <- F

# model
if (lambda_fixed == T){
  mf <- textConnection("model {
    # priors
    lambda <- 1
    alpha ~ dunif(0.01,2) 
    v ~ dunif(0.01,2)
    d ~ dunif(0.01,2)

    for (i in 1:N) {
      x[i] ~ dwald_gamma(lambda, alpha, v, d)
      #x[i] ~ dwald_trunc(lambda, alpha, v, d)
    }
  }")
  
  inits1 <- list(alpha=.7, v=.7, d=.7)
  inits2 <- list(alpha=.9, v=.9, d=.9)
  inits3 <- list(alpha=1.1, v=1.1, d=1.1)
  
  params <- c("alpha", "v", "d")
} else {
  mf <- textConnection("model {
    # priors
    lambda ~ dunif(0.01,2) 
    alpha ~ dunif(0.01,2) 
    v ~ dunif(0.01,2)
    d ~ dunif(0.01,2)

    for (i in 1:N) {
      x[i] ~ dwald_gamma(lambda, alpha, v, d)
    }
  }")
  inits1 <- list(lambda=.7, alpha=.7, v=.7, d=.7)
  inits2 <- list(lambda=.9, alpha=.9, v=.9, d=.9)
  inits3 <- list(lambda=1.1, alpha=1.1, v=1.1, d=1.1)
  
  params <- c("lambda", "alpha", "v", "d")
}

# Genarate data
# (1) Generate drift rates
N <- 10000
d <- 1
v <- 1
nu <- rtruncnorm(N, a=0, b=Inf, mean=d, sd=sqrt(v))  
# (2) Generate RTs
alpha=lambda=1
RT=rinvgauss(N, mean=alpha/nu, shape=lambda*alpha^2)
plot(density(RT[RT<=3]))
x <- RT
N <- length(x)
dat <- list(x=x, N=N)

# inits
inits <- list(inits1,inits2,inits3)

# sample
j.model <- jags.model(mf, dat, inits, n.chains=3, n.adapt=1000)
j.samples <- coda.samples(j.model, params, n.iter=4000, thin=3)

# plot
par(mfrow=c(3,4))
plot(j.samples)


save.image("~/Dropbox/2014/Wald/TruncNorm/Results/samples_lambda_free.rdata")


########
alpha <- .965
d <- 1.34
v <- .015
lambda <- 1
nu <- rtruncnorm(N, a=0, b=Inf, mean=d, sd=sqrt(v))  
RT2 <- rinvgauss(N, mean=alpha/nu, shape=lambda*alpha^2)

plot(density(RT[RT<=3]))
lines(density(RT2[RT2<=3]))

