library(tidyverse)
library(rmatio)
library(brms)
library(latex2exp)
library(ggthemes)
library(rstan)

#------------------------------------------------------------------------------------#
# DATA PREPARATION
#------------------------------------------------------------------------------------#
setwd("/Users/lukas/Documents/UniHeidel/Project_Discrimination/data_deJong2020")
# read all files into one data frame
df = NULL
for (f in list.files()) {
  
  tmp = read.mat(f)
  colnames <- unlist(tmp$VariableNames)
  data <- as.tibble(tmp$ResultMatrix)
  colnames(data) = colnames
  
  data$id <- tmp$SubjDetails$Name
  data$sex <- tmp$SubjDetails$Sex
  data$age <- unlist(tmp$SubjDetails$Age)
  data$trial <- 1:nrow(data)
  
  df=rbind(df, data)
}

# tidy data
df <- df %>% 
  rename(rt=RT,
         correct=Acc,
         resp=button,
         d1=StimDur,
         d2=Probe_dur,
         comparison_diff=Probe_Dur_diff,
         comparison_dir=Probe_Dur_dir) %>% 
  select(-trlOnset,
         -trlEnd) %>% 
  relocate(id,sex,age,trial,
           comparison_dir,
           comparison_diff,
           Dur_step,
           d1,d2,resp,
           correct,rt) %>% 
  mutate(id=as.numeric(id))

# exclude some participants
df <- df %>% 
  group_by(id) %>% 
  mutate(included=ifelse(max(comparison_diff) < 1 & max(trial)==200,T,F)) %>% 
  filter(included==T) %>% 
  select(-included)

# resp: 0 = c is shorter; 1 = c is longer
# center Standard and Comparison around the geometric mean
df <- df %>% 
  mutate(resp=resp-1,
         d1_centered=d1-exp(mean(log(d1))),
         d2_centered=d2-exp(mean(log(d2))))

# reassign id nr
df %<>%
  mutate(id=as.numeric(id)) %>% 
  arrange(id)

df$id[df$id==2] <- 1

nr <- 2
for (i in 2:length(df$id)){
  if(df$id[i]!=df$id[i-1]){
    df$id[df$id==df$id[i]] <- nr
    nr <- nr+1
  }
}
# signed difference between S and C
df$delta_d <- df$comparison_diff*df$comparison_dir


#------------------------------------------------------------------------------------#
# MODEL FITTING
#------------------------------------------------------------------------------------#
# cut RT's above 5 sec
df <- df %>% 
  filter(rt<=5)

S = length(unique(df$id))

# arrange rows of df after sub and then trial
df <- df %>% 
  arrange(id,trial)

# create stan data
stan_data = list(
  `T`  = nrow(df),
  S    = S,
  sub  = df$id,
  d1   = df$d1_centered,
  d2   = df$d2_centered,
  resp = df$resp,
  rt   = df$rt
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
    
    L[[c]]$sd_a   = 0.001
    L[[c]]$sd_ndt = 0.001
    L[[c]]$sd_z0  = 0.001
    L[[c]]$sd_bz  = 0.001
    L[[c]]$sd_v0  = 0.001
    L[[c]]$sd_bv  = 0.001
    L[[c]]$sd_b1v = 0.001
    L[[c]]$sd_b2v = 0.001
    L[[c]]$sd_g   = 0.001
    
    L[[c]]$a = rep(1,S)
    L[[c]]$z0 = rep(0.5,S)
    L[[c]]$bz = rep(0.0,S)
    L[[c]]$v0 = rep(1.0,S)
    L[[c]]$bv = rep(0.0,S)
    L[[c]]$b1v = rep(0.0,S)
    L[[c]]$b2v = rep(0.0,S)
    L[[c]]$ndt = rep(0.1,S)
    L[[c]]$g = rep(0.5,S)
  }
  return (L)
}

setwd("/users/lukas/documents/github/timing_discrimination/models")
m7 <-  stan("m7_deJong.stan",
            init=init(4),
            data=stan_data,
            chains=4,
            iter = 500,
            cores=parallel::detectCores())

