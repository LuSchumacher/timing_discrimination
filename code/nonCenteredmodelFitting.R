library(tidyverse)
library(brms)
library(rstan)
library(magrittr)
library(rstan)
library(bayesplot)

# read data
setwd("/users/lukas/documents/github/timing_discrimination/data")
df <- read_csv("exp1_cond_C.csv")


##------------------------------------------------
## FITTING PREPARATION
##------------------------------------------------
# for quicker estimation
df %<>% filter(sub<10)

# number of participants and trials
S <- length(unique(df$sub))
`T` <- nrow(df)

# create stan data
stan_data = list(
  `T`     = `T`,
  S       = S,
  sub     = df$sub,
  d1      = df$d1-500,
  d2      = df$d2-500,
  resp    = df$resp,
  rt      = df$RT / 1000
)

init = function(chains=4) {
  L = list()
  for (c in 1:chains) {
    L[[c]]=list()
    
    L[[c]]$mu_a       = 0.5
    L[[c]]$mu_ndt     = -2.5
    # L[[c]]$mu_ndt_var = -4
    L[[c]]$mu_z0      = rnorm(1,0,1)
    L[[c]]$mu_v0      = rnorm(1,0,1)
    L[[c]]$mu_b1v     = rnorm(1,0,1)
    L[[c]]$mu_b2v     = rnorm(1,0,1)
    L[[c]]$ndt_var    = 0.001
    L[[c]]$eps        = rep(0.001,`T`)
    

    L[[c]]$sd_a       = 0.001
    L[[c]]$sd_ndt     = 0.001
    # L[[c]]$sd_ndt_var = 0.001
    L[[c]]$sd_z0      = 0.001
    L[[c]]$sd_v0      = 0.001
    L[[c]]$sd_b1v     = 0.001
    L[[c]]$sd_b2v     = 0.001

    L[[c]]$z_a        = abs(rnorm(S,0,1))
    L[[c]]$z_ndt      = rnorm(S,0,1)
    # L[[c]]$z_ndt_var  = rnorm(1,0,1)
    L[[c]]$z_z0       = rnorm(S,0,1)
    L[[c]]$z_v0       = rnorm(S,0,1)
    L[[c]]$z_b1v      = rnorm(S,0,1)
    L[[c]]$z_b2v      = rnorm(S,0,1)
  }
  return (L)
}

##------------------------------------------------
## FIT MODEL
##------------------------------------------------
setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/code")
nonC_m3 <-  stan("m3_nonCentered.stan",
                     init=init(4),
                     data=stan_data,
                     chains=4,
                     iter = 500,
                     cores=parallel::detectCores(),
                     control = list(adapt_delta=0.95))

params <- c("mu_a","mu_ndt","mu_z0","mu_v0","mu_b1v","mu_b2v")
mcmc_pairs(nonC_m3,pars=params)
