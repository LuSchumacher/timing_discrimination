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
    
    L[[c]]$mu_a   = 1.0
    L[[c]]$mu_ndt = 0.1
    L[[c]]$mu_z0  = 0.5
    L[[c]]$mu_bz  = 0.0
    L[[c]]$mu_v0  = 1.0
    L[[c]]$mu_bv  = 0.0
    L[[c]]$mu_b1v = 0.0
    L[[c]]$mu_b2v = 0.0
    L[[c]]$mu_g   = 0.5
    L[[c]]$mu_w   = 0.5
    
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
    
    L[[c]]$a   = rep(1,S)
    L[[c]]$z0  = rep(0.5,S)
    L[[c]]$bz  = rep(0.0,S)
    L[[c]]$v0  = rep(1.0,S)
    L[[c]]$bv  = rep(0.0,S)
    L[[c]]$b1v = rep(0.0,S)
    L[[c]]$b2v = rep(0.0,S)
    L[[c]]$ndt = rep(0.1,S)
    L[[c]]$g   = rep(0.5,S)
    L[[c]]$w   = rep(0.5,S)
  }
  return (L)
}

##------------------------------------------------
## FIT MODEL
##------------------------------------------------
setwd("/users/lukas/documents/github/timing_discrimination/models")
fit_cond_C_m6 <-  stan("hierarchical_m6.stan",
                         init=init(2),
                         data=stan_data,
                         chains=2,
                         iter = 500,
                         cores=parallel::detectCores(),
                         control = list(adapt_delta=0.95))

loo_fit_cond_C_m6 <- loo(fit_cond_C_m6)
# saveRDS(fit_cond_C_m4,"/users/lukas/documents/UniHeidel/Project_Discrimination/fits/fit_cond_C_m4.rds")




