#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# bern.R
#
# (c) 2013 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2013-03-19
# last mod 2015-02-14 17:33 DW
#

rm(list=ls())

# necessary libs
library(rjags)
library(truncnorm)
library(statmod)
load.module("wald")


# model
my.model <- function() {
  # priors
  lambda <- 1  # is fixed!

  alpha ~ dunif(0, 10) 
  v ~ dunif(0, 10)
  d ~ dunif(0, 10)

  for (i in 1:N) {
    x[i] ~ dwald_trunc(lambda, alpha, v, d)
  }
}

# inits  
inits1 <- list(alpha=0.7, v=0.7, d=0.9)
inits2 <- list(alpha=0.9, v=0.9, d=1.0)
inits3 <- list(alpha=1.1, v=1.1, d=1.1)
inits <- list(inits1,inits2,inits3)

# Parameters to be monitored  
params <- c("alpha", "d", "v")

# Genarate data
# (1) Generate drift rates
N <- 10000
d <- 9  # mean of the truncated normal distribution
v <- (.5)^2  # nu; std of the trunc. normal distribution 
nu <- rtruncnorm(N, a=0, b=Inf, mean=d, sd=sqrt(v))  
plot(density(nu))

# (2) Generate RTs
# alpha: boundary separation
# lambda: diffusion coefficient
alpha=1.2
lambda=1 
tr.params <- c(alpha, d, v)
RT=rinvgauss(N, mean=alpha/nu, shape=lambda*alpha^2)
plot(density(RT))

dat <- list(x=RT, N=N)

# sample
samples <- jags(dat, inits=inits, params,  # inits=NULL
	 			model.file=my.model, n.chains=3, n.iter=2000, 
                n.burnin=1000, n.thin=1, DIC=T)  
samples
j.samples <- as.mcmc(samples)

index <- c(1, 2, 4)
tr.params <- c(alpha, d, v)
layout(matrix(1:3, nrow=1))
for (i in 1:3) {  # loop over chains
  h <- index[i]
  samples <- c(j.samples[[1]][,h], j.samples[[2]][,h], j.samples[[3]][,h])
  plot(density(samples), main=params[i])
  lines(c(tr.params[i],tr.params[i]), c(0, max(density(samples)$y)))  	
} 
