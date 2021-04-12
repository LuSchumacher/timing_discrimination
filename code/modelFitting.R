library(tidyverse)
library(brms)
library(rstan)
library(magrittr)

# read data
setwd("/users/lukas/documents/github/timing_discrimination/data")
df <- read_csv("exp1_cond_C.csv")

##------------------------------------------------
## FITTING PREPARATION
##------------------------------------------------

# number of participants
S = length(unique(df$sub))
`T` = nrow(df)

# create stan data
stan_data = list(
  `T`  = nrow(df),
  S    = S,
  sub  = df$sub,
  d1   = df$d1-500,
  d2   = df$d2-500,
  resp = df$resp,
  rt   = df$RT / 1000
)

# set initial values
init = function(chains=4) {
  L = list()
  for (c in 1:chains) {
    L[[c]]=list()
    
    L[[c]]$mu_a   = runif(1,0.5,2)
    L[[c]]$mu_ndt = runif(1,0.08,0.1)
    L[[c]]$mu_z0  = rnorm(1,0.5,0.1)
    L[[c]]$mu_bz  = 0.0
    L[[c]]$mu_v0  = runif(1,0.5,2)
    L[[c]]$mu_bv  = 0.0
    L[[c]]$mu_b1v = 0.0
    L[[c]]$mu_b2v = 0.0
    L[[c]]$mu_g   = rnorm(1,0.5,0.1)
    L[[c]]$mu_w   = rnorm(1,0.5,0.1)
    # L[[c]]$ndt_sd_mu   = 0.01
    L[[c]]$ndt_var   = runif(`T`,0.08,0.1)
    
    L[[c]]$sd_a   = 0.001
    L[[c]]$sd_ndt = 0.001
    L[[c]]$sd_z0  = 0.001
    L[[c]]$sd_bz  = 0.001
    L[[c]]$sd_v0  = 0.001
    L[[c]]$sd_bv  = 0.001
    L[[c]]$sd_b1v = 0.001
    L[[c]]$sd_b2v = 0.001
    L[[c]]$sd_g   = 0.001
    L[[c]]$sd_w   = 0.001
    L[[c]]$sd_ndt_sd   = 0.001
    
    L[[c]]$a   = runif(S,0.5,2)
    L[[c]]$ndt = runif(S,0.08,0.1)
    L[[c]]$z0  = rnorm(S,0.5,0.1)
    L[[c]]$bz  = rep(0.0,S)
    L[[c]]$v0  = runif(S,0.5,2)
    L[[c]]$bv  = rep(0.0,S)
    L[[c]]$b1v = rep(0.0,S)
    L[[c]]$b2v = rep(0.0,S)
    L[[c]]$g   = rnorm(S,0.5,0.1)
    L[[c]]$w   = rnorm(S,0.5,0.1)
    # L[[c]]$ndt_sd   = rep(0.001,S)
  }
  return (L)
}

##------------------------------------------------
## FIT MODEL
##------------------------------------------------
setwd("/users/lukas/documents/github/timing_discrimination/models")
fit_tbtVarm7 <-  stan("trial-by-trial_var_m7.stan",
                         init=init(4),
                         data=stan_data,
                         chains=4,
                         iter = 500,
                         cores=parallel::detectCores(),
                         control = list(adapt_delta=0.95))

# saveRDS(fit_m11,"/users/lukas/documents/UniHeidel/Project_Discrimination/fits/fit_m11_new.rds")


fit_tbtVarm7
# ndt_sd_mu
params <- c("mu_a","mu_ndt", "mu_z0","mu_bz","mu_v0","mu_b1v","mu_b2v","ndt_sd_mu")
mcmc_pairs(fit_tbtVarm7,pars=params)



