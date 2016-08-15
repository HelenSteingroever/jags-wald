################################################################################
##### Fit individual-level IG-GAM mixture
##### Sample-based solution
##### Helen Steingroever
##### Last updated: August 2016
################################################################################

rm(list=ls())
# Load required libraries
#library(rjags)
library(R2jags)
library(statmod)     # for rinvgauss 
#library(SuppDists)  # for rinvGauss
#library(msm)	     # for rtnorm
#library(RWald)
load.module("wald")  # for the JAGS distribution dswald()
                     # To activate the extension, modules have to be loaded
                     # explicetly. 
################################################################################

################################################################################
##### Genarate data
################################################################################
T <- 120      # number of trials

# (1) Generate drift rates
mu  <- 9.42  # drift mean
var <- 4.62  # drift variance
kappa <- mu^2 / var  # shape
tau <- mu  / var     # rate           
xi <- rgamma(T, shape=kappa, rate=tau)

# (2) Generate RTs
alpha <- 0.90  # threshold  
theta <- 0.10  # shift
RT <- rinvgauss(T, mean=alpha/xi, shape=alpha^2) + theta
minRT <- min(RT)
################################################################################

################################################################################ 
##### JAGS
################################################################################
# Parameters to be monitored:	
params <- c("alpha", "theta", "mu_xi", "sig2_xi")

# Data to be passed on to JAGS
dat <- list(RT=RT, T=T, minRT=minRT)

# Collect samples from posterior distributions
# Change the working directory to the folder containing the model file
setwd("~/Dropbox/2014/Wald/SWmixturesCode")
samples <- jags(dat, inits=NULL, params,  # inits=NULL
	 			model.file="Model_IG_GAM_ind.txt", n.chains=3, n.iter=31000, 
	 			n.burnin=1000, n.thin=15, DIC=T)  
# Print a summary of the posterior samples
samples
################################################################################  
  