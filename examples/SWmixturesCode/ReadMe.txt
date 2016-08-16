################################################################################
# Helen Steingroever, Last updated August 2016
# This text file lists all files that are necessary to fit the SW-TN and SW-GAM 
# mixtures as described in Steingroever, Wabersich, & Wagenmakers (submitted).
# The programs R and JAGS (a recent version) are needed to fit the models. In
# addition, the jags-wald module has to be build (see Wabersich & 
# Vandekerckhove, 2014; Behav Res; they also describe how the installation of
# the module can be tested).
# Please cite "Steingroever, Wabersich & Wagenmakers (submitted). Modeling 
# Across-Trial Variability in the Wald Drift Rate Parameter" when using this code.
################################################################################

1. Fit_IG_TN_ind.r:  R code to generate data and to fit the IG-TN mixture.
2. Fit_IG_GAM_ind.r: R code to generate data and to fit the IG-GAM mixture.
3. Model_IG_TN_ind.txt:  JAGS model file of the IG-TN mixture. This file is
   needed to run the Fit_IG_TN_ind.r code.
4. Model_IG_GAM_ind.txt: JAGS model file of the IG-GAM mixture. This file is 
   needed to run the Fit_IG_GAM_ind.r code.
