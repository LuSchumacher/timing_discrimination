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
    
    L[[c]]$mu_a     = rnorm(1,0.5,0.1)
    L[[c]]$mu_ndt   = rnorm(1,-2.5,0.1)
    L[[c]]$mu_z0    = rnorm(1,0.5,0.1)
    L[[c]]$mu_bz    = 0
    L[[c]]$mu_v0    = rnorm(1,0,0.5)
    L[[c]]$mu_b1v   = rnorm(1,0,0.5)
    L[[c]]$mu_b2v   = rnorm(1,0,0.5)
    L[[c]]$mu_g     = rnorm(1,0,0.5)
    
    L[[c]]$sigma_a     = 0.001
    L[[c]]$sigma_ndt   = 0.001
    L[[c]]$sigma_z0    = 0.001
    L[[c]]$sigma_bz    = 0.001
    L[[c]]$sigma_v0    = 0.001
    L[[c]]$sigma_b1v   = 0.001
    L[[c]]$sigma_b2v   = 0.001
    L[[c]]$sigma_g     = 0.001
    
    S = length(unique(df$sub))
    L[[c]]$z_a     = rnorm(S,0,0.5)
    L[[c]]$z_ndt   = rnorm(S,0,0.5)
    L[[c]]$z_z0    = rnorm(S,0,0.1)
    L[[c]]$z_bz    = rnorm(S,0,0.5)
    L[[c]]$z_v0    = rnorm(S,0,0.5)
    L[[c]]$z_b1v   = rnorm(S,0,0.5)
    L[[c]]$z_b2v   = rnorm(S,0,0.5)
    L[[c]]$z_g     = rnorm(S,0,0.5)
    
    # L[[c]]$z0      = rep(0.5,S)
    # L[[c]]$bz      = rep(0.0,S)
    
  }
  return (L)
}

setwd("/users/lukas/documents/github/timing_discrimination/models")
m9 <-  stan("m9.stan",
            init=init(4),
            data=stan_data,
            chains=4,
            iter = 500,
            cores=parallel::detectCores(),
            control = list(adapt_delta=0.9))

