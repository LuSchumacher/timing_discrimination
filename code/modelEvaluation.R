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
m7 <- readRDS("fit_cond_C_m7.rds")
m8 <- readRDS("fit_cond_C_m8.rds")
m9 <- readRDS("fit_cond_C_H_m9[actual].rds")
m12 <- readRDS("fit_cond_C_m12.rds")

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

