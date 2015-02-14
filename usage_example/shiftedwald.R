#!/usr/bin/R --silent -f
# -*- encoding: utf-8 -*-
# shiftedwald.R
#
# (c) 2015 Dominik Wabersich <dominik.wabersich [aet] gmail.com>
# GPL 3.0+ or (cc) by-sa (http://creativecommons.org/licenses/by-sa/3.0/)
#
# created 2015-02-14
# last mod 2015-02-14 17:49 DW
#

# Genarate data
N <- 100

alpha <- 1
nu <- 1
theta <- 0.3

rswald <- function(N, alpha, nu, theta) {
  beta.par <- 0.5
  x <- rwiener(N, alpha/(1-beta.par), theta, beta.par, nu)
  return(x)
}

RT <- rswald(N, alpha, nu, theta)
