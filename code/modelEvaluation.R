library(tidyverse)
library(brms)
library(rstan)
library(magrittr)
library(lemon)
library(ggthemes)
library(latex2exp)

# read data
setwd("/users/lukas/documents/github/timing_discrimination/data")
df <- read_csv("exp1_cond_C.csv")

# read model fits
setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/fits")
m6 <- readRDS("fit_cond_C_m6.rds")
m7 <- readRDS("fit_cond_C_m7.rds")
m8 <- readRDS("fit_cond_C_m8.rds")
m9 <- readRDS("fit_cond_C_H_m9[actual].rds")
m12 <- readRDS("fit_cond_C_m12.rds")

shinystan::launch_shinystan(m6)

check_divergences(m7)
check_divergences(m8)
check_divergences(m9)
check_divergences(m12)

# model comparison
loo_m7 <- loo(m7)
loo_m8 <- loo(m8)
loo_m9 <- loo(m9)
loo_m12 <- loo(m12)

comparison <- loo_compare(loo_m7,
                          loo_m8,
                          loo_m9,
                          loo_m12)



setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/fits")
m1 <- readRDS("fit_cond_C_m1.rds")
m2 <- readRDS("fit_cond_C_m2.rds")
m3 <- readRDS("fit_cond_C_m3.rds")
m4 <- readRDS("fit_cond_C_m4.rds")
m5 <- readRDS("fit_cond_C_m5.rds")
m6 <- readRDS("fit_cond_C_m6.rds")

loo_m1 <- loo(m1)
loo_m2 <- loo(m2)
loo_m3 <- loo(m3)
loo_m6 <- loo(m6)

compare <- loo_compare(loo_m1,loo_m1_transformed)

comparison <- loo_compare(loo_m1,
                          loo_m2,
                          loo_m3,
                          loo_m6,
                          loo_m6_new)

pars_m1 <- c("mu_a","mu_ndt","mu_z0","mu_v0","mu_b1v")
mcmc_pairs(m1,
           pars = pars_m1)

check_divergences(m1)

print(m4, pars="mu_bz",digits=8)


setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/fits")
m1_new <- readRDS("fit_m1_new.rds")
m2_new <- readRDS("fit_m3_new.rds")
m3_new <- readRDS("fit_m3_new.rds")
m6_new <- readRDS("fit_m6_new.rds")


loo_m1_new <- loo(m1_new)
loo_m2_new <- loo(m2_new)
loo_m3_new <- loo(m3_new)
loo_m6_new <- loo(m6_new)

comparison <- loo_compare(loo_m1_new,
                          loo_m2_new,
                          loo_m3_new,
                          loo_m6_new)
