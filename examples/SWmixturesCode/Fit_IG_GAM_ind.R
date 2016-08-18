################################################################################
##### Fit individual-level IG-GAM mixture
##### Sample-based solution
##### Helen Steingroever
##### Last updated: August 2016
################################################################################

rm(list=ls())
# Load required libraries
library(R2jags)
library(statmod)     # for rinvgauss 
load.module("wald")  # for the JAGS distribution dswald()
                     # To activate the extension, modules have to be loaded
                     # explicetly. 
################################################################################

################################################################################
##### Genarate data
################################################################################
set.seed(3)
T <- 120      # number of trials

# (1) Generate drift rates
mu_xi   <- 10.04  # drift mean
sig2_xi <-  4.54  # drift variance
kappa <- mu_xi^2 / sig2_xi  # shape
tau <- mu_xi / sig2_xi    # rate           
xi <- rgamma(T, shape=kappa, rate=tau)

# (2) Generate RTs
alpha <- 0.98  # threshold  
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

# Change the working directory to the folder containing the model file
setwd("examples/SWmixturesCode")

# Collect samples from posterior distributions
samples <- jags(dat, inits=NULL, params,  # inits=NULL
	 			model.file="Model_IG_GAM_ind.txt", n.chains=3, n.iter=41000, 
	 			n.burnin=1200, n.thin=10, DIC=T)  
# Print a summary of the posterior samples
samples
################################################################################  
  
