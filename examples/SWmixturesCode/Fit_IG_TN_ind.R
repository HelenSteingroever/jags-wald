################################################################################
##### Fit individual-level IG-TN mixture
##### Analytical solution
##### Helen Steingroever
##### Last updated: August 2016
################################################################################

rm(list=ls())
# Load required libraries
library(R2jags)
library(statmod)     # for rinvgauss 
library(msm)	     # for rtnorm
load.module("wald")  # for the JAGS distribution dwald_trunc()
                     # To activate the extension, modules have to be loaded
                     # explicetly. 
################################################################################

################################################################################
##### Genarate data
################################################################################
T <- 120  # number of trials

# (1) Generate drift rates
mu_xi   <- 8.93  # drift mean
sig2_xi <- 4.05  # drift variance
xi <- rtnorm(T, mean=mu_xi, sd=sqrt(sig2_xi), lower=0)

# (2) Generate RTs
alpha <- 0.71  # threshold  
theta <- 0.12  # shift
RT <- rinvgauss(T, mean=alpha/xi, shape=alpha^2) + theta
minRT <- min(RT)
################################################################################

################################################################################ 
##### JAGS
################################################################################
# Parameters to be monitored  
params <- c("alpha", "theta", "mu_xi", "sig2_xi")

# Data to be passed on to JAGS
dat <- list(RT=RT, T=T, minRT=minRT)

# Change the working directory to the folder containing the model file
setwd("examples/SWmixturesCode")

# Collect samples from posterior distributions
samples <- jags(dat, inits=NULL, params,  # inits=NULL
	 			model.file="Model_IG_TN_ind.txt", n.chains=3, n.iter=31000, 
                n.burnin=1000, n.thin=15, DIC=T) 
# Print a summary of the posterior samples
samples                  
################################################################################  
   
