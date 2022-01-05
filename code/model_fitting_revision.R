library(tidyverse)
library(rstan)
library(magrittr)


df <- read_csv("/users/lukas/documents/UniHeidel/Project_Discrimination/data_exp2_condD_revision.csv")

# get number of participants
S = length(unique(df$sub))
`T`  = nrow(df)

# arrange rows of df after sub and then trial
df <- df %>% 
  arrange(sub, trl)

# create stan data
stan_data = list(
  `T`    = `T`,
  S      = S,
  sub    = df$sub,
  min_rt = min(df$rt / 1000),
  d1     = (df$d1 - 500) / 100,
  d2     = (df$d2 - 500) / 100,
  resp   = df$resp,
  rt     = df$rt / 1000
)


# set initial values
init = function(chains=4) {
  L = list()
  for (c in 1:chains) {
    L[[c]] = list()
    
    L[[c]]$mu_a    = 1.0
    L[[c]]$mu_ndt  = 3
    L[[c]]$mu_z0   = 0.5
    L[[c]]$mu_bz   = 0.0
    L[[c]]$mu_v0   = 1.0
    L[[c]]$mu_bv   = 0.0
    L[[c]]$mu_b1v  = 0.0
    L[[c]]$mu_b2v  = 0.0
    
    L[[c]]$sd_a    = 0.001
    L[[c]]$sd_ndt  = 0.001
    L[[c]]$sd_z0   = 0.001
    L[[c]]$sd_bz   = 0.001
    L[[c]]$sd_v0   = 0.001
    L[[c]]$sd_bv   = 0.001
    L[[c]]$sd_b1v  = 0.001
    L[[c]]$sd_b2v  = 0.001
    
    L[[c]]$ndt_var = 15
    
    L[[c]]$a       = rep(1.0, S)
    L[[c]]$z0      = rep(0.5, S)
    L[[c]]$bz      = rep(0.0, S)
    L[[c]]$v0      = rep(1.0, S)
    L[[c]]$bv      = rep(0.0, S)
    L[[c]]$b1v     = rep(0.0, S)
    L[[c]]$b2v     = rep(0.0, S)
    L[[c]]$ndt     = rep(3.0, S)
    L[[c]]$ndt_t   = rep(0.09, `T`)
  }
  return (L)
}

setwd("/Users/lukas/Documents/UniHeidel/Project_Discrimination/code")
fit <-  stan("model_6_revision.stan",
             init = init(4),
             data = stan_data,
             chains = 4,
             iter = 2000,
             cores = parallel::detectCores()
             # control = list(adapt_delta=0.99))
            )
