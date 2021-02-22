library(tidyverse)
library(brms)
library(rstan)
library(magrittr)
library(bayesplot)

# read data
setwd("/users/lukas/documents/github/timing_discrimination/data")
df <- read_csv("exp1_cond_C.csv")

##------------------------------------------------
## DATA SUBSET
##------------------------------------------------
df <- df %>% 
  filter(sub<5)

##------------------------------------------------
## FITTING PREPARATION
##------------------------------------------------

# create stan data
stan_data = list(
  `T`  = nrow(df),
  sub  = df$sub,
  d1   = df$d1,
  d2   = df$d2,
  resp = df$resp,
  rt   = df$RT / 1000
)

# set initial values
init = function(chains=4) {
  L = list()
  for (c in 1:chains) {
    L[[c]]=list()
    
    L[[c]]$a   = 1.0
    L[[c]]$z0  = 0.5
    L[[c]]$v0  = 1.0
    L[[c]]$b1v = 0.0
    L[[c]]$ndt = 0.1
    
    L[[c]]$noise_sd = 5
    L[[c]]$kalman_q = 9

  }
  return (L)
}

##------------------------------------------------
## FIT MODEL
##------------------------------------------------
setwd("/users/lukas/documents/github/timing_discrimination/models")
fit_kalman_diff <-  stan("kalman_diffusion_test_2.stan",
                         init=init(4),
                         data=stan_data,
                         chains=4,
                         iter = 500,
                         cores=parallel::detectCores(),
                         control = list(adapt_delta=0.95))

setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/fits")
saveRDS(fit_kalman_diff,"fit_kalman_diff_actual.rds")

pars_fit_kalman_diff <- c("a","ndt","z0","v0","b1v","kalman_q","noise_sd")

print(fit_kalman_diff, pars=pars_fit_kalman_diff,digits=4)

mcmc_pairs(fit_kalman_diff,
           pars = pars_fit_kalman_diff)


pars_fit_kalman_diff <- c("delta[1]","delta[2]","delta[3]","delta[4]","delta[5]","delta[100]")
print(fit_kalman_diff,pars=pars_fit_kalman_diff)

pars_fit_kalman_diff <- c("kalman_k[2]","kalman_k[3]","kalman_k[4]","kalman_k[5]","kalman_k[6]","kalman_k[100]")
print(fit_kalman_diff,pars=pars_fit_kalman_diff)

pars_fit_kalman_diff <- c("I_1[1]","I_1[2]","I_1[3]","I_1[4]","I_1[5]","I_1[100]")
print(fit_kalman_diff,pars=pars_fit_kalman_diff)

pars_fit_kalman_diff <- c("I_2[1]","I_2[2]","I_2[3]","I_2[4]","I_2[5]","I_2[100]")
print(fit_kalman_diff,pars=pars_fit_kalman_diff)

mat <- rstan::extract(fit_kalman_diff, pars=pars_fit_kalman_diff)
mat <- as.data.frame(mat)


##----------------------------------------------##
##     NON-CENTERED HIERARCHICAL KALMAN DDM     ##
##----------------------------------------------##
# create stan data
stan_data = list(
  S    = length(unique(df$sub)),
  `T`  = nrow(df),
  sub  = df$sub,
  d1   = df$d1,
  d2   = df$d2,
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
    L[[c]]$mu_z0    = rnorm(1,0,1)
    L[[c]]$mu_v0    = rnorm(1,0,1)
    L[[c]]$mu_b1v   = rnorm(1,0,1)
    L[[c]]$mu_q     = rnorm(1,0,1)
    L[[c]]$mu_noise = rnorm(1,0,1)
    
    L[[c]]$sd_a     = 0.001
    L[[c]]$sd_ndt   = 0.001
    L[[c]]$sd_z0    = 0.001
    L[[c]]$sd_v0    = 0.001
    L[[c]]$sd_b1v   = 0.001
    L[[c]]$sd_q     = 0.001
    L[[c]]$sd_noise = 0.001
    
    S = length(unique(df$sub))
    L[[c]]$z_a     = rnorm(S,0,1)
    L[[c]]$z_ndt   = rnorm(S,0,1)
    L[[c]]$z_z0    = rnorm(S,0,1)
    L[[c]]$z_v0    = rnorm(S,0,1)
    L[[c]]$z_b1v   = rnorm(Ss,0,1)
    L[[c]]$z_q     = rnorm(S,0,1)
    L[[c]]$z_noise = rnorm(S,0,0.5)
    
  }
  return (L)
}

setwd("/users/lukas/documents/github/timing_discrimination/models")
fit_kalman_diff <-  stan("hierarchical_kalman_DDM.stan",
                         init=init(4),
                         data=stan_data,
                         chains=4,
                         iter = 500,
                         cores=parallel::detectCores(),
                         control = list(adapt_delta=0.95))

setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/fits")
fit_kalman_diff <- readRDS("fit_kalman_diff_actual.rds")

params <- c("mu_a","mu_ndt","mu_z0","mu_v0","mu_b1v","mu_q","mu_noise")
mcmc_pairs(fit_kalman_diff,
           pars=params)

params <- c("mu_a","sd_a","mu_ndt","sd_ndt","mu_z0","sd_z0","mu_v0","sd_v0","mu_b1v","sd_b1v","mu_q","sd_q","mu_noise","sd_noise")
traceplot(fit_kalman_diff,pars=params)

