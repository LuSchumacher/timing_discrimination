library(tidyverse)
library(rmatio)
library(brms)
library(latex2exp)
library(ggthemes)
library(rstan)
library(bayesplot)
library(rmatio)
#------------------------------------------------------------------------------------#
# LOAD EMPIRICAL DATA
#------------------------------------------------------------------------------------#
setwd("/Users/lukas/Documents/UniHeidel/Project_Discrimination/code")
df <- read.csv("deJong_final_data.csv")

# cut RT's above 5 sec
# arrange rows of df after sub and then trial
df <- df %>% 
  filter(rt<=5) %>% 
  arrange(id,trial)

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

# setwd("/Users/lukas/Documents/UniHeidel/Project_Discrimination/code")
# write_csv(df,"deJong_final_data.csv")
#------------------------------------------------------------------------------------#
# MODEL FITTING: M7
#------------------------------------------------------------------------------------#
S = length(unique(df$id))

# create stan data
stan_data = list(
  `T`  = nrow(df),
  S    = S,
  sub  = df$id,
  d1   = df$d1_centered,
  d2   = df$d2_centered,
  resp = df$resp+1,
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
m7_2 <-  stan("m7_deJong_2.stan",
            init=init(4),
            data=stan_data,
            chains=4,
            iter = 500,
            cores=parallel::detectCores())

#------------------------------------------------------------------------------------#
# MODEL EVALUATION: M7
#------------------------------------------------------------------------------------#
# read fitted model
setwd("/Users/lukas/Documents/UniHeidel/Project_Discrimination/fits")
m7 <- readRDS("fit_dejong_m7_new.rds")


# plot marginal posteriors
pars <- c("mu_a","mu_ndt","mu_z0","mu_bz","mu_v0","mu_b1v","mu_b2v")
mcmc_pairs(m7, pars=pars)

# extract posterior samples
mat <- as.data.frame(rstan::extract(m7, pars=pars))

# sample parameter sets
mat$index <- 1:nrow(mat)
mat <- mat[sample(max(mat$index),10),]

# experiment settings
nTrials <- nrow(df)
d1 = df$d1_centered
d2 = df$d2_centered
d1_nC = df$d1
d2_nC = df$d2
delta_d <- df$delta_d

# simulate new data
pp_sim_m7 <- function(nTrials, mat){
  
  setwd("/Users/lukas/documents/UniHeidel/code")
  source("wienerProcess2.R")
  library(svMisc)
  
  data <- data.frame(matrix(ncol=8,nrow=0, dimnames=list(NULL, c("resp","rt","d1","d2","d1_nC","d2_nC","delta_d","dataset"))))
  
  a <- mat[,1]
  ndt <- mat[,2]
  z <- mat[,3]
  bz <- mat[,4]
  v0 <- mat[,5]
  b1v <- mat[,6]
  b2v <- mat[,7]
  
  for (t in 1:nrow(mat)) {
    
    beta = z[t] + bz[t] * d1
    delta <-  v0[t] + b1v[t] * d1 + b2v[t] * d2
    
    start_time <- Sys.time()
    for (i in 1:nTrials) {
      data <- data %>% add_row(wienerProcess2(v=delta[i], a=a[t], z=beta[i], ndt=ndt[t]),
                               d1=d1[i],d2=d2[i],d1_nC=d1_nC[i],d2_nC=d2_nC[i],delta_d=delta_d[i],dataset=t)
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    
  }
  
  return(data)
}

pred_m7 <- pp_sim_m7(nTrials, mat)

#------------------------------------------------------------------------------------#
# SUMMARIZE SIMULATED DATA
#------------------------------------------------------------------------------------#
# sumsum <- df %>%
#   group_by(d1,
#            delta_d) %>%
#   summarise(n=length(resp))
summary_pred_m7 <- pred_m7 %>% 
  group_by(dataset,
           delta_d,
           d1_nC) %>% 
  mutate(d2_new=mean(d2_nC)) %>% 
  ungroup() %>% 
  mutate(resp=ifelse(resp==0,1,0)) %>% 
  group_by(dataset,
           d1,
           d2_new) %>% 
  summarise(acc=mean(resp),
            n=length(resp)) %>% 
  ungroup() %>% 
  group_by(d1,
           d2_new) %>% 
  summarise(acc=mean(acc),
            n=n[1])
  mutate(d1=as.factor(d1)) %>% 
  filter(n>10)


# summarize within each dataset
summary_pred_m7 <- pred_m7 %>% 
  mutate(resp=ifelse(resp==0,1,0)) %>% 
  group_by(dataset,
           d1,
           delta_d) %>% 
  summarise(acc=mean(resp),
            rt=median(rt),
            n=length(resp))

# summarize over datasets
summary_pred_m7 <- summary_pred_m7 %>% 
  group_by(d1,
           delta_d) %>%
  summarise(acc_quantile_low=quantile(acc,probs=0.025),
            acc_quantile_high=quantile(acc,probs=0.975),
            acc=mean(acc),
            rt_quantile_low=quantile(rt,probs=0.025),
            rt_quantile_high=quantile(rt,probs=0.975),
            rt=mean(rt),
            n=n[1]) %>% 
  mutate(d1=as.factor(d1))
  # filter(n>10)


#------------------------------------------------------------------------------------#
# SUMMARIZE EMPIRICAL DATA
#------------------------------------------------------------------------------------#
summary_df <- df %>% 
  #mutate(resp=ifelse(resp==0,1,0)) %>% 
  group_by(d1_centered,
           delta_d) %>% 
  summarise(acc=mean(resp),
            rt=median(rt),
            n=length(resp)) %>% 
  mutate(d1_centered=as.factor(d1_centered)) %>% 
  rename(d1=d1_centered)
  # filter(n>10)

summary_kalman_sim <- sim_dt_kalman %>% 
  #mutate(resp=ifelse(resp==0,1,0)) %>% 
  group_by(s_centered,
           delta) %>% 
  summarise(acc=mean(longer)) %>% 
  mutate(s_centered=as.factor(s_centered)) %>% 
  rename(d1=s_centered)


#------------------------------------------------------------------------------------#
# PLOTTING
#------------------------------------------------------------------------------------#
color_palette <- c("#D8E0BB","#B6CEC7","#B6A3C3","#7268A6")

# plot different d1 in individual facet
summary_pred_m7 %>% 
  ggplot(aes(x=delta_d,
             y=acc))+
  # geom_point()+
  # geom_line()+
  # geom_smooth(se=F)+
  geom_ribbon(aes(ymin = acc_quantile_low,
                  ymax = acc_quantile_high),
              alpha=0.5,
              color=NA)+
  geom_point(data=summary_df,
             mapping=aes(x=delta_d,
                         y=acc))+
  geom_line(data=summary_df,
             mapping=aes(x=delta_d,
                         y=acc))+
  ggthemes::theme_tufte()+
  scale_fill_manual(values = color_palette)+
  scale_color_manual(values = color_palette)+
  facet_wrap(~d1)

# plot different d1 in individual color
summary_pred_m7 %>% 
  ggplot(aes(x=delta_d,
             y=acc,
             color=d1))+
  # geom_point()+
  # geom_line()+
  # geom_smooth(se=F)+
  geom_ribbon(aes(ymin = acc_quantile_low,
                  ymax = acc_quantile_high,
                  fill=d1),
              alpha=0.5,
              color=NA)+
  geom_point(data=summary_df,
             mapping=aes(x=delta_d,
                         y=acc,
                         color=d1))+
  geom_line(data=summary_df,
            mapping=aes(x=delta_d,
                        y=acc,
                        color=d1))+
  ggthemes::theme_tufte()+
  scale_fill_manual(values = color_palette)+
  scale_color_manual(values = color_palette)

# emp data only
summary_df %>%
  filter(n>10) %>% 
  ggplot()+
  geom_point(aes(x=delta_d,
                 y=acc,
                 color=d1))+
  geom_line(aes(x=delta_d,
                y=acc,
                color=d1))

# emp data only
summary_kalman_sim %>%
  # filter(n>10) %>% 
  ggplot()+
  geom_point(aes(x=delta,
                 y=acc,
                 color=d1))+
  geom_line(aes(x=delta,
                y=acc,
                color=d1))

summary_df_2 <- df %>% 
  group_by(delta_d,
           d1) %>% 
  mutate(d2_new=mean(d2)) %>% 
  ungroup() %>% 
  group_by(d1,
           d2_new) %>% 
  summarise(acc=mean(resp),
            n=length(resp)) %>% 
  mutate(d1=as.factor(d1)) %>% 
  filter(n>5)

summary_sim_dt_kalman <- sim_dt_kalman %>% 
  group_by(c,s) %>% 
  summarise(acc=mean(longer),
            n=length(longer)) %>% 
  mutate(s=as.factor(s)) %>% 
  filter(n>10)

summary_sim_dt_IRM2 <- sim_dt_IRM1 %>% 
  group_by(c,s) %>% 
  summarise(acc=mean(longer),
            n=length(longer)) %>% 
  mutate(s=as.factor(s)) %>% 
  filter(n>10)

summary_pred_m3 <- pred_m3_3 %>% 
  group_by(delta_d,
           d1) %>% 
  mutate(d2_new=mean(d2)) %>% 
  ungroup() %>% 
  group_by(d1,
           d2_new) %>% 
  summarise(acc=mean(resp),
            n=length(resp)) %>% 
  mutate(d1=as.factor(d1)) %>% 
  filter(n>10)


summary_df_2 %>%
  ggplot(aes(x=d2_new,
             y=acc,
             color=d1))+
  geom_point(shape=3)+
  geom_line(linetype="dashed")+
  # scale_x_continuous(limits = c(0,2.4),
  #                    breaks = c(0.3,0.6,1.2,2.4))+
  # geom_point(data = summary_sim_dt_IRM2,
  #            mapping = aes(x=c,
  #                          y=acc,
  #                          color=s))+
  # geom_line(data = summary_sim_dt_IRM2,
  #           mapping = aes(x=c,
  #                         y=acc,
  #                         color=s))+
  # geom_point(data = summary_sim_dt_kalman,
  #            mapping = aes(x=c,
  #                          y=acc,
  #                          color=s))+
  # geom_line(data = summary_sim_dt_kalman,
  #            mapping = aes(x=c,
  #                          y=acc,
  #                          color=s))+
  # geom_point(data = summary_pred_m7,
  #            mapping = aes(x=d2_new,
  #                          y=acc,
  #                          group=d1),
  # color="black")+
  # geom_line(data = summary_pred_m7,
  #           mapping = aes(x=d2_new,
  #                         y=acc,
  #                         group=d1),
  #           color="black")+
  # geom_point(data = summary_pred_m3,
  #            mapping = aes(x=d2_new,
  #                          y=acc,
  #                          group=d1),
  #            color="maroon")+
  geom_line(data = summary_pred_m3,
            mapping = aes(x=d2_new,
                          y=acc,
                          group=d1),
            color="maroon")

summary_df <- df %>% 
  #mutate(resp=ifelse(resp==0,1,0)) %>% 
  group_by(d1,
           delta_d) %>% 
  summarise(acc=mean(resp),
            rt=median(rt),
            n=length(resp)) %>% 
  mutate(d1=as.factor(d1))

  
ggplot(data = dat, aes(ymin = 0)) +
  scale_color_brewer(palette = "RdBu", labels = c("0.3", "0.6", "1.2", "2.4")) +
  scale_fill_brewer(palette = "RdBu", labels = c("0.3", "0.6", "1.2", "2.4")) +
  geom_line(data=dat_newdata, aes(x=Comparison, y=fit, colour=factor(Standard)), size = 1.5) +
  #geom_segment(aes(x = lower, y = 0.5, xend = upper, yend = 0.5), colour = "black", data = HDIs_PSE, show.legend = FALSE, size = 2) +
  # geom_segment(data=PSE_dt, aes(x=s, xend=PSE, y=0.5, yend=0.5, color=factor(s)), size=2) +
  # geom_segment(data=PSE_dt, aes(x=s, xend=s, y=0.5, yend=0, color=factor(s)), size=1, linetype='dashed') +
  # geom_segment(aes(x=geometric_mean, xend=geometric_mean, y=0.5, yend=0), size=1, linetype='dashed') +
  # geom_density(data=CE_dt, 
  #              aes(x=PSE, y=0.04*..density.., 
  #                  fill=as.factor(Standard), color=as.factor(Standard)), 
  #              alpha=0.15) +
  # geom_point(data=PSE_dt, aes(x=s, y=0.5, color=factor(s)), size=3) +
  # geom_point(aes(x=geometric_mean, y=0.5), size=3) +
  guides(fill = 'none', colour = guide_legend(reverse = FALSE, title = "Standard"), size = FALSE) +
  scale_x_continuous(name = "", breaks = c(0.3, 0.6, 1.2, 2.4), limits = c(0,2.4)) +
  scale_y_continuous(name = "P('C Longer')", breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  #ggtitle('Empirical') +
  theme(legend.position = 'none', 
        legend.title = element_text(size = 16),
        #plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
        axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14))+
  geom_point(data=summary_df_2,
             aes(x=d2_new,
                 y=acc,
                 color=d1),
             size=2)+
  geom_line(data=summary_df_2,
            aes(x=d2_new,
                y=acc,
                color=d1),
            linetype="dashed",
            size=1)+
  theme_tufte()
  
  
  
  
#------------------------------------------------------------------------------------#
# FIT LOGISTIC REGRESSION
#------------------------------------------------------------------------------------#
df$d1_centered <- as.factor(df$d1_centered)
emp_logReg <- brm(resp~delta_d*d1_centered,
                  data=df,
                  family = bernoulli(),
                  chains = 2,
                  iter = 500,
                  cores = parallel::detectCores())

conditional_effects(emp_logReg)


pred_m7_swtich <- pred_m7 %>% 
  mutate(resp=ifelse(resp==0,1,0),
         d1=as.factor(d1))

pred_logReg <- brm(resp~delta_d*d1,
                   data=pred_m7_swtich,
                   family = bernoulli(),
                   chains = 2,
                   iter = 500,
                   cores = parallel::detectCores())

conditional_effects(pred_logReg)




#------------------------------------------------------------------------------------#
# NON-CENTERED
#------------------------------------------------------------------------------------#
# summarize within each dataset
summary_pred_m7 <- pred_m7 %>% 
  mutate(resp=ifelse(resp==0,1,0)) %>% 
  group_by(dataset,
           d1_nC,
           d2_nC) %>% 
  summarise(acc=mean(resp),
            rt=median(rt),
            n=length(resp))

# summarize over datasets
summary_pred_m7 <- summary_pred_m7 %>% 
  group_by(d1_nC,
           d2_nC) %>%
  summarise(acc_quantile_low=quantile(acc,probs=0.025),
            acc_quantile_high=quantile(acc,probs=0.975),
            acc=mean(acc),
            rt_quantile_low=quantile(rt,probs=0.025),
            rt_quantile_high=quantile(rt,probs=0.975),
            rt=mean(rt),
            n=n[1]) %>% 
  mutate(d1_nC=as.factor(d1_nC)) %>% 
  filter(n>10)


#------------------------------------------------------------------------------------#
# SUMMARIZE EMPIRICAL DATA
#------------------------------------------------------------------------------------#
summary_df <- df %>% 
  #mutate(resp=ifelse(resp==0,1,0)) %>% 
  group_by(d1,
           d2) %>% 
  summarise(acc=mean(resp),
            rt=median(rt),
            n=length(resp)) %>% 
  mutate(d1=as.factor(d1)) %>% 
  filter(n>10)


#------------------------------------------------------------------------------------#
# PLOTTING
#------------------------------------------------------------------------------------#
color_palette <- c("#D8E0BB","#B6CEC7","#B6A3C3","#7268A6")

# plot different d1 in individual color
summary_pred_m7 %>% 
  ggplot(aes(x=d2_nC,
             y=acc,
             color=d1_nC))+
  # geom_point()+
  # geom_line()+
  # geom_smooth(se=F)+
  geom_ribbon(aes(ymin = acc_quantile_low,
                  ymax = acc_quantile_high,
                  fill=d1_nC),
              alpha=0.5,
              color=NA)+
  geom_point(data=summary_df,
             mapping=aes(x=d2,
                         y=acc,
                         color=d1))+
  geom_line(data=summary_df,
            mapping=aes(x=d2,
                        y=acc,
                        color=d1))+
  ggthemes::theme_tufte()+
  scale_fill_manual(values = color_palette)+
  scale_color_manual(values = color_palette)

# emp data only
summary_df %>%
  filter(n>10) %>% 
  ggplot()+
  geom_point(aes(x=d2,
                 y=acc,
                 color=d1))+
  geom_line(aes(x=d2,
                y=acc,
                color=d1))

df$d1 <- as.factor(df$d1)
emp_logReg <- brm(resp~d2*d1,
                  data=df,
                  family = bernoulli(),
                  chains = 2,
                  iter = 500,
                  cores = parallel::detectCores())

conditional_effects(emp_logReg,
                    spaghetti = T)
color_palette <- c("#D8E0BB","#B6CEC7","#B6A3C3","#7268A6")
plot(p1,plot = FALSE)[[3]]+
  geom_point(data=summary_df,
             aes(x=d2,
                 y=acc,
                 color=d1))+
  geom_line(data=summary_df,
            aes(x=d2,
                y=acc,
                color=d1))

pred_m7_swtich <- pred_m7 %>% 
  mutate(resp=ifelse(resp==0,1,0),
         d1_nC=as.factor(d1_nC))

pred_logReg <- brm(resp~d2_nC*d1_nC,
                   data=pred_m7_swtich,
                   family = bernoulli(),
                   chains = 2,
                   iter = 500,
                   cores = parallel::detectCores())

conditional_effects(pred_logReg,
                    spaghetti = T,
                    points=T)

#------------------------------------------------------------------------------------#
# MODEL EVALUATION: M3
#------------------------------------------------------------------------------------#
# read fitted model
setwd("/Users/lukas/Documents/UniHeidel/Project_Discrimination/fits")
m3 <- readRDS("deJong_m3_new.RDS")

# plot marginal posteriors
pars <- c("mu_a","mu_ndt","mu_z0","mu_bz","mu_v0","mu_b1v","mu_g")
mcmc_pairs(m3, pars=pars)

# extract posterior samples
mat <- as.data.frame(rstan::extract(m3, pars=pars))

# sample parameter sets
mat$index <- 1:nrow(mat)
mat <- mat[sample(max(mat$index),10),]

# mean posteriors
mat <- summary(m3,pars=pars)$summary[,"mean"]


# experiment settings
nTrials <- nrow(df)
d1 = df$d1_centered
d2 = df$d2_centered
d1_nC = df$d1
d2_nC = df$d2
delta_d <- df$delta_d

# simulate new data
pp_sim_m3 <- function(nTrials, mat){
  
  setwd("/Users/lukas/documents/UniHeidel/code")
  source("wienerProcess2.R")
  library(svMisc)
  
  data <- data.frame(matrix(ncol=9,nrow=0, dimnames=list(NULL, c("resp","rt","d1","d2","d1_nC","d2_nC","delta_d","X_1","X_2"))))
  
  a <- mat[1]
  ndt <- mat[2]
  z <- mat[3]
  bz <- mat[4]
  v0 <- mat[5]
  b1v <- mat[6]
  g <- mat[7]
  
  X_1 <- c()
  X_2 <- c()
  
  start_time <- Sys.time()
  for (i in 1:nTrials) {

    if (i == 1) {
      I <-  d1[1]
      X_1[1] <- d1[1]
      X_2[1] <- d2[1]
    } else if (df$id[i] != df$id[i-1]) {
      I <- d1[i]
      X_1[i] <- d1[i]
      X_2[i] <- d2[i]
    } else{
      I <- g*I + (1-g)*d1[i]
      X_1[i] <- I
      I <- g*I + (1-g)*d2[i]
      X_2[i] <- I
    }

    beta <- z + bz * X_1
    delta <-  v0 + b1v * (X_1 - X_2)
    
    
    data <- data %>% add_row(wienerProcess2(v=delta[i], a=a, z=beta[i], ndt=ndt),
                             d1=d1[i],d2=d2[i],d1_nC=d1_nC[i],d2_nC=d2_nC[i],delta_d=delta_d[i],X_1=X_1[i],X_2=X_2[i])
  }
  end_time <- Sys.time()
  print(end_time - start_time)
    
  
  return(data)
}

pred_m3 <- pp_sim_m3(nTrials, mat)


pred_m3$d1 <- as.factor(pred_m3$d1)
pred_logReg <- brm(resp~d2_nC*d1,
                   data=pred_m3,
                   family = bernoulli(),
                   chains = 2,
                   iter = 500,
                   cores = parallel::detectCores())

conditional_effects(pred_logReg)

#------------------------------------------------------------------------------------#
# SUMMARIZE SIMULATED DATA
#------------------------------------------------------------------------------------#
# sumsum <- df %>%
#   group_by(d1,
#            delta_d) %>%
#   summarise(n=length(resp))

# summarize over datasets
summary_pred_m3 <- pred_m3 %>% 
  # mutate(resp=ifelse(resp==0,1,0)) %>% 
  group_by(d1,
           delta_d) %>%
  summarise(acc=mean(resp),
            rt=median(rt)) %>% 
  mutate(d1=as.factor(d1))
# filter(n>10)


#------------------------------------------------------------------------------------#
# SUMMARIZE EMPIRICAL DATA
#------------------------------------------------------------------------------------#
summary_df <- df %>% 
  #mutate(resp=ifelse(resp==0,1,0)) %>% 
  group_by(d1_centered,
           delta_d) %>% 
  summarise(acc=mean(resp),
            rt=median(rt),
            n=length(resp)) %>% 
  mutate(d1_centered=as.factor(d1_centered)) %>% 
  rename(d1=d1_centered)
# filter(n>10)


#------------------------------------------------------------------------------------#
# PLOTTING
#------------------------------------------------------------------------------------#
color_palette <- c("#D8E0BB","#B6CEC7","#B6A3C3","#7268A6")

# plot different d1 in individual facet
summary_pred_m3 %>% 
  ggplot(aes(x=delta_d,
             y=acc))+
  geom_point()+
  geom_line()+
  # geom_smooth(se=F)+
  geom_point(data=summary_df,
             mapping=aes(x=delta_d,
                         y=acc),
             shape=3)+
  geom_line(data=summary_df,
            mapping=aes(x=delta_d,
                        y=acc),
            linetype="dashed")+
  ggthemes::theme_tufte()+
  scale_fill_manual(values = color_palette)+
  scale_color_manual(values = color_palette)+
  facet_wrap(~d1)

# plot different d1 in individual color
summary_pred_m3 %>% 
  ggplot(aes(x=delta_d,
             y=acc,
             color=d1))+
  # geom_point()+
  # geom_line()+
  # geom_smooth(se=F)+
  geom_ribbon(aes(ymin = acc_quantile_low,
                  ymax = acc_quantile_high,
                  fill=d1),
              alpha=0.5,
              color=NA)+
  geom_point(data=summary_df,
             mapping=aes(x=delta_d,
                         y=acc,
                         color=d1))+
  geom_line(data=summary_df,
            mapping=aes(x=delta_d,
                        y=acc,
                        color=d1))+
  ggthemes::theme_tufte()+
  scale_fill_manual(values = color_palette)+
  scale_color_manual(values = color_palette)

#------------------------------------------------------------------------------------#
# SUBJECT PREDICTION
#------------------------------------------------------------------------------------#
individual_g <- summary(m3,pars=c("g[1]","g[2]","g[3]","g[4]","g[5]","g[6]","g[7]","g[8]","g[9]","g[10]",
                                  "g[11]","g[12]","g[13]","g[14]","g[15]","g[16]","g[17]","g[18]","g[19]","g[20]",
                                  "g[21]","g[22]","g[23]","g[24]","g[25]","g[26]","g[27]","g[28]","g[29]","g[30]",
                                  "g[31]","g[32]","g[33]","g[34]","g[35]","g[36]","g[37]","g[38]","g[39]","g[40]",
                                  "g[41]","g[42]","g[43]","g[44]","g[45]","g[46]"))$summary[,"mean"]
  
mat <- as.data.frame(m3) %>%
  select(starts_with(c("a[","ndt[","z0[","bz[","v0[","b1v[","g[")))

numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
}

a <- mat %>% 
  select(`a[1]`:`a[46]`) %>% 
  pivot_longer(cols = `a[1]`:`a[46]`,
               names_to = "sub",
               values_to = "a") %>% 
  mutate(sample=rep(seq(1,nrow(mat)),each=nrow(.)/nrow(mat)),
         sub=numextract(sub)) %>% 
  select(sample,sub,a)

ndt <- mat %>% 
  select(`ndt[1]`:`ndt[46]`) %>% 
  pivot_longer(cols = `ndt[1]`:`ndt[46]`,
               names_to = "sub",
               values_to = "ndt") %>% 
  select(ndt)

z0 <- mat %>% 
  select(`z0[1]`:`z0[46]`) %>% 
  pivot_longer(cols = `z0[1]`:`z0[46]`,
               names_to = "sub",
               values_to = "z0") %>% 
  select(z0)

bz <- mat %>% 
  select(`bz[1]`:`bz[46]`) %>% 
  pivot_longer(cols = `bz[1]`:`bz[46]`,
               names_to = "sub",
               values_to = "bz") %>% 
  select(bz)

v0 <- mat %>% 
  select(`v0[1]`:`v0[46]`) %>% 
  pivot_longer(cols = `v0[1]`:`v0[46]`,
               names_to = "sub",
               values_to = "v0") %>% 
  select(v0)

b1v <- mat %>% 
  select(`bz[1]`:`b1v[46]`) %>% 
  pivot_longer(cols = `b1v[1]`:`b1v[46]`,
               names_to = "sub",
               values_to = "b1v") %>% 
  select(b1v)

g <- mat %>% 
  select(`g[1]`:`g[46]`) %>% 
  pivot_longer(cols = `g[1]`:`g[46]`,
               names_to = "sub",
               values_to = "g") %>% 
  select(g)


mat <- cbind(a,ndt,z0,bz,v0,b1v,g)

# sample posteriors
posterior_samples <- sample(max(mat$sample),1)
mat <- mat %>% 
  filter(sample %in% posterior_samples)




sub_sim <- function(mat){
  
  setwd("/Users/lukas/documents/UniHeidel/code")
  source("wienerProcess.R")
  library(svMisc)
  
  data <- data.frame(matrix(ncol=11,nrow=0,
                            dimnames=list(NULL, c("dataset","sub","d1","d2","d1_nC","d2_nC","delta_d","X_1","X_2","resp", "rt"))))
  
  # iterate over posterior samples
  for (t in unique(mat$sample)){
    
    
    # iterate over subjects
    for (s in 1:length(unique(mat$sub))){
      a   <- mat$a[mat$sample==t][s]
      z   <- mat$z0[mat$sample==t][s]
      bz  <- mat$bz[mat$sample==t][s]
      v0  <- mat$v0[mat$sample==t][s]
      b1v <- mat$b1v[mat$sample==t][s]
      g <- mat$g[mat$sample==t][s]
      ndt <- mat$ndt[mat$sample==t][s]
      
      tmp_df <- df %>% 
        filter(id==s)
      
      nTrials <- nrow(tmp_df)
      d1 = tmp_df$d1_centered
      d2 = tmp_df$d2_centered
      d1_nC = tmp_df$d1
      d2_nC = tmp_df$d2
      delta_d <- tmp_df$delta_d
      
      X_1 <- c()
      X_2 <- c()
      start_time <- Sys.time() 
      for (i in 1:nTrials) {

        if (i == 1) {
          I <-  d1[1]
          X_1[1] <- d1[1]
          X_2[1] <- d2[1]
        } else if (df$id[i] != df$id[i-1]) {
          I <- d1[i]
          X_1[i] <- d1[i]
          X_2[i] <- d2[i]
        } else{
          I <- g*I + (1-g)*d1[i]
          X_1[i] <- I
          I <- g*I + (1-g)*d2[i]
          X_2[i] <- I
        }
        
        beta <- z + bz * X_1[i]
        delta <-  v0 + b1v * (X_1[i] - X_2[i])
        
        data <- data %>%
          add_row(dataset=t, sub=s,
                  d1=d1[i],d2=d2[i],d2_nC=d2_nC[i],d1_nC=d1_nC[i],
                  delta_d=delta_d[i],X_1=X_1[i],X_2=X_2[i],
                  wienerProcess(v=delta, a=a, z=beta, ndt=ndt))
      }
      end_time <- Sys.time()
      print(end_time - start_time)
    }

  }
  
  return(data)
}

pred_sub_m3 <- sub_sim(mat = mat)


# summarize within dataset
summary_pred_sub_m3 <- pred_sub_m3 %>% 
  group_by(sub,
           d1,
           delta_d) %>% 
  summarise(acc=mean(resp),
            rt=median(rt))  %>% 
  rename(id=sub,
         d1_centered=d1) %>% 
  mutate(d1_centered=as.factor(d1_centered))

# summarize within dataset
summary_pred_sub_m3 <- pred_sub_m3 %>% 
  group_by(dataset,
           sub,
           d1,
           delta_d) %>% 
  summarise(acc=mean(resp),
            rt=median(rt)) %>% 
  mutate(d1=as.factor(d1))

# summarize over datasets
summary_pred_sub_m7 <- summary_pred_sub_m7 %>% 
  group_by(sub,
           cpos,
           cdur) %>%
  summarise(acc_quantile_low=quantile(acc,probs=0.025),
            acc_quantile_high=quantile(acc,probs=0.975),
            acc=mean(acc),
            rt_quantile_low=quantile(rt,probs=0.025),
            rt_quantile_high=quantile(rt,probs=0.975),
            rt=median(rt))

summary <- df %>%
  group_by(id,
           d1_centered,
           delta_d) %>% 
  mutate(d1_centered=as.factor(d1_centered)) %>% 
  summarise(prob_c = mean(resp),
            rt=median(rt),
            n=length(resp))

summary$prob_c[summary$cpos==2] <- 1 - summary$prob_c[summary$cpos==2]
summary$cpos <- as.factor(summary$cpos)



summary_pred_sub_m3 %>% 
  ggplot(aes(x=delta_d,
             y=acc,
             color=d1_centered,
             fill=d1_centered)) +
  geom_point()+
  geom_line()+
  # geom_ribbon(aes(ymin = acc_quantile_low,
  #                 ymax = acc_quantile_high,
  #                 fill=cpos),
  #             alpha=0.2,
  #             color=NA)+
  geom_point(data=summary,
             mapping=aes(x=delta_d,
                         y=prob_c,
                         color=d1_centered),
             shape=3)+
  geom_line(data=summary,
            mapping=aes(x=delta_d,
                        y=prob_c,
                        color=d1_centered),
            alpha=0.5,
            linetype="dashed")+
  facet_wrap(~id,nrow=4)+
  ggthemes::theme_tufte()+
  scale_color_manual(values = color_palette)+
  scale_fill_manual(values = color_palette,
                    guide=F)+
  ylab("Probability for c > s response")+
  xlab("\nDuration of c")+
  labs(color = "Position of c")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(size = 16),
        plot.title = element_text(size =18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))

  
#------------------------------------------------------------------------------------#
# FIT LOGISTIC REGRESSION
#------------------------------------------------------------------------------------#
df$d1_centered <- as.factor(df$d1_centered)
emp_logReg <- brm(resp~d2*d1_centered,
                  data=df,
                  family = bernoulli(),
                  chains = 2,
                  iter = 500,
                  cores = parallel::detectCores())

conditional_effects(emp_logReg)


pred_m7_swtich <- pred_sub_m3 %>% 
  mutate(resp=ifelse(resp==0,1,0),
         d1=as.factor(d1))


pred_sub_m3$d1 <- as.factor(pred_sub_m3$d1)

pred_logReg <- brm(resp~d2*d1,
                   data=pred_sub_m3,
                   family = bernoulli(),
                   chains = 2,
                   iter = 500,
                   cores = parallel::detectCores())

conditional_effects(pred_logReg)


#------------------------------------------------------------------------------------#
# MODEL FITTING: KALMAN MODEL
#------------------------------------------------------------------------------------#
# create stan data
stan_data = list(
  `T`  = nrow(df),
  sub  = df$id,
  d1   = df$d1_centered,
  d2   = df$d2_centered,
  resp = df$resp+1,
  rt   = df$rt
)


# set initial values
init = function(chains=4) {
  L = list()
  for (c in 1:chains) {
    L[[c]]=list()
    
    L[[c]]$a   = 1.0
    L[[c]]$z0  = 0.5
    L[[c]]$v0  = 1.0
    L[[c]]$bv  = 0.0
    L[[c]]$b1v = 0.0
    L[[c]]$b2v = 0.0
    L[[c]]$ndt = 0.1
    
    L[[c]]$kalman_q = 9
    
  }
  return (L)
}

setwd("/users/lukas/documents/github/timing_discrimination/models")
fit_kalman_diff <-  stan("kalman_diffusion_test.stan",
                         init=init(4),
                         data=stan_data,
                         chains=4,
                         iter = 500,
                         cores=parallel::detectCores(),
                         control = list(adapt_delta=0.95))

setwd("/Users/lukas/Documents/UniHeidel/Project_Discrimination/fits")
saveRDS(fit_kalman_diff,"fit_kalman_diff_new.rds")
# fit_kalman_diff <- readRDS("fit_kalman_ddm.RDS")


params <- c("a","ndt","z0","v0","b1v","kalman_q")
mcmc_pairs(fit_kalman_diff,params)

mat <- summary(fit_kalman_diff,pars=params)$summary[,"mean"]

# experiment settings
nTrials <- nrow(df)
d1 = df$d1_centered
d2 = df$d2_centered
delta_d <- df$delta_d


kalman_ddm_sim <- function(nTrials, mat){
  
  setwd("/Users/lukas/documents/UniHeidel/code")
  source("wienerProcess.R")
  library(svMisc)
  
  data <- data.frame(matrix(ncol=9,nrow=0, dimnames=list(NULL, c("resp","rt","delta","d1","I_1","d2","I_2","delta_d","kalman_k"))))
  
  # Fitted posterior parameters
  a <- mat[1]
  ndt <- mat[2]
  z <- mat[3]
  v0 <- mat[4]
  b1v <- mat[5]
  kalman_q <- mat[6]
  noise_sd <- 2
  
  # variable initialization
  xd1 <- rep(NA,nTrials)
  xd2 <- rep(NA,nTrials)
  I_1 <- rep(NA,nTrials)
  I_2 <- rep(NA,nTrials)
  I <- NA
  
  kalman_k <- rep(NA,nTrials)
  kalman_prior_var <- rep(NA,nTrials)
  
  # sampling perceptive noise
  percept_noise_1 <- rnorm(nTrials,0,noise_sd)
  percept_noise_2 <- rnorm(nTrials,0,noise_sd)
  
  # iterate over trials and simulate
  for (i in 1:nTrials) {
    
    # initial values on trial 1
    if(i==1){
      I_1[1] <-  d1[1] + percept_noise_1[1]
      I_2[1] <-  d2[1] + percept_noise_2[1]
      I = I_1[1]
      kalman_prior_var[1] <-  5
    }else{
      # noisy stimulus perception
      xd1[i] = d1[i] + percept_noise_1[i];
      xd2[i] = d2[i] + percept_noise_2[i];
      
      # final percept stimulus 1
      kalman_k[i] <- (kalman_prior_var[i-1] + kalman_q) / (kalman_prior_var[i-1] + kalman_q + noise_sd)
      kalman_prior_var[i] <- kalman_k[i] * noise_sd
      I <- (1-kalman_k[i]) * I + kalman_k[i] * xd1[i]
      I_1[i] <- I
      
      # final percept stimulus 2
      kalman_k[i] <- (kalman_prior_var[i-1] + kalman_q) / (kalman_prior_var[i-1] + kalman_q + noise_sd)
      kalman_prior_var[i] <- kalman_k[i] * noise_sd
      I <- (1-kalman_k[i]) * I + kalman_k[i] * xd2[i]
      I_2[i] <- I
    }
    
    # trial-by-trial drift
    delta <- v0 + b1v * (I_1[i] - I_2[i])
    
    # diffusion process
    data <- data %>%
      add_row(wienerProcess(v=delta, a=a, z=z, ndt=ndt),
              delta=delta,d1=d1[i],I_1=I_1[i],d2=d2[i],I_2=I_2[i],delta_d=delta_d[i],kalman_k=kalman_k[i])
    
    # update progress bar
    # setTxtProgressBar(pb,i)
  }
  
  return(data)
}

prediction <- kalman_ddm_sim(nTrials,mat)
prediction$d1 <- as.factor(prediction$d1)


# create a summary of the empirical data
df_summary <- df %>% 
  group_by(d1_centered,
           delta_d) %>% 
  mutate(d2_new=mean(d2_centered)) %>% 
  ungroup() %>% 
  group_by(d1_centered,
           d2_new) %>% 
  summarise(acc = mean(resp))

df_summary$d1_centered <- as.factor(df_summary$d1_centered)


# Summarize accuracy of predicted data
summary_prediction <- prediction %>% 
  group_by(d1,
           delta_d) %>% 
  mutate(d2_new=mean(d2)) %>% 
  ungroup() %>% 
  group_by(d1,
           d2_new) %>% 
  summarise(acc=mean(resp)) %>% 
  mutate(acc=1-acc)

summary_prediction$d1 <- as.factor(summary_prediction$d1)

summary_prediction %>%
  ggplot(aes(x=d2_new,
             y = acc,
             color=d1))+
  geom_line(size=0.8)+
  geom_point(size=1)+
  geom_point(data=df_summary,
             mapping=aes(x=d2_new,
                         y=acc,
                         color=d1_centered),
             shape=3,
             size=4)+
  geom_line(data=df_summary,
            mapping=aes(x=d2_new,
                        y=acc,
                        color=d1_centered),
            linetype="dashed")+
  ggthemes::theme_tufte()
  scale_color_manual(values = color_palette)



loo_kalman <- loo(fit_kalman_diff)

#------------------------------------------------------------------------------------#
# MODEL COMPARISON
#------------------------------------------------------------------------------------#
loo_m7 <- loo(m7) 
loo_m3 <- loo(m3)
loo_m3_2 <- loo(m3_2)

loo::loo_compare(loo_m3,loo_m7,loo_m3_2)


# read fitted model
setwd("/Users/lukas/Documents/UniHeidel/Project_Discrimination/fits")
m3_2 <- readRDS("deJong_m3_new_2.RDS")

# plot marginal posteriors
pars <- c("mu_a","mu_ndt","mu_z0","mu_bz","mu_v0","mu_b1v","mu_g")
mcmc_pairs(m3_2, pars=pars)

# extract posterior samples
mat <- as.data.frame(rstan::extract(m3_2, pars=pars))

# sample parameter sets
mat$index <- 1:nrow(mat)
mat <- mat[sample(max(mat$index),10),]

# mean posteriors
mat <- summary(m3_3,pars=pars)$summary[,"mean"]


# experiment settings
nTrials <- nrow(df)
d1 = df$d1
d2 = df$d2
delta_d <- df$delta_d

# simulate new data
pp_sim_m3_2 <- function(nTrials, mat){
  
  setwd("/Users/lukas/documents/UniHeidel/code")
  source("wienerProcess2.R")
  library(svMisc)
  
  data <- data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c("resp","rt","d1","d2","delta_d","X_1","X_2"))))
  
  a <- mat[1]
  ndt <- mat[2]
  z <- mat[3]
  bz <- mat[4]
  v0 <- mat[5]
  b1v <- mat[6]
  g <- mat[7]
  
  X_1 <- c()
  X_2 <- c()
  
  start_time <- Sys.time()
  for (i in 1:nTrials) {
    
    if (i == 1) {
      I <-  d1[1]
      X_1[1] <- d1[1]
      X_2[1] <- d2[1]
    } else if (df$id[i] != df$id[i-1]) {
      I <- d1[i]
      X_1[i] <- d1[i]
      X_2[i] <- d2[i]
    } else{
      I <- g*I + (1-g)*d1[i]
      X_1[i] <- I
      I <- g*I + (1-g)*d2[i]
      X_2[i] <- I
    }
    
    beta <- z + bz * X_1
    delta <-  v0 + b1v * ((X_1 - X_2)/(X_1 + X_2))
    
    
    data <- data %>% add_row(wienerProcess2(v=delta[i], a=a, z=beta[i], ndt=ndt),
                             d1=d1[i],d2=d2[i],delta_d=delta_d[i],X_1=X_1[i],X_2=X_2[i])
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  
  
  return(data)
}

pred_m3_3 <- pp_sim_m3_2(nTrials, mat)

pred_m3_3$d1 <- as.factor(pred_m3_3$d1)
pred_logReg <- brm(resp~d2*d1,
                   data=pred_m3_3,
                   family = bernoulli(),
                   chains = 2,
                   iter = 500,
                   cores = parallel::detectCores())

conditional_effects(pred_logReg)

#------------------------------------------------------------------------------------#
# SUMMARIZE SIMULATED DATA
#------------------------------------------------------------------------------------#
# sumsum <- df %>%
#   group_by(d1,
#            delta_d) %>%
#   summarise(n=length(resp))

# summarize over datasets
summary_pred_m3_3 <- pred_m3_3 %>% 
  # mutate(resp=ifelse(resp==0,1,0)) %>% 
  group_by(d1,
           delta_d) %>%
  summarise(acc=mean(resp),
            rt=median(rt),
            n=length(resp)) %>% 
  mutate(d1=as.factor(d1)) %>% 
  filter(n>5)

summary_kalman_sim <- sim_dt_kalman %>% 
  #mutate(resp=ifelse(resp==0,1,0)) %>% 
  group_by(s_centered,
           delta) %>% 
  summarise(acc=mean(longer)) %>% 
  mutate(s_centered=as.factor(s_centered)) %>% 
  rename(d1=s_centered)
#------------------------------------------------------------------------------------#
# SUMMARIZE EMPIRICAL DATA
#------------------------------------------------------------------------------------#
summary_df <- df %>% 
  #mutate(resp=ifelse(resp==0,1,0)) %>% 
  group_by(d1,
           delta_d) %>% 
  summarise(acc=mean(resp),
            rt=median(rt),
            n=length(resp)) %>% 
  mutate(d1=as.factor(d1)) %>% 
  filter(n>5)


#------------------------------------------------------------------------------------#
# PLOTTING
#------------------------------------------------------------------------------------#
color_palette <- c("#D8E0BB","#B6CEC7","#B6A3C3","#7268A6")

# plot different d1 in individual facet
summary_pred_m3_3 %>% 
  ggplot(aes(x=delta_d,
             y=acc))+
  geom_point()+
  geom_line()+
  # geom_smooth(se=F)+
  geom_point(data=summary_df,
             mapping=aes(x=delta_d,
                         y=acc),
             shape=3)+
  geom_line(data=summary_df,
            mapping=aes(x=delta_d,
                        y=acc),
            linetype="dashed")+
  ggthemes::theme_tufte()+
  scale_fill_manual(values = color_palette)+
  scale_color_manual(values = color_palette)+
  facet_wrap(~d1)

# plot different d1 in individual color
summary_pred_m3_3 %>% 
  ggplot(aes(x=delta_d,
             y=acc,
             color=d1))+
  geom_point()+
  geom_line()+
  geom_point(data=summary_df,
             mapping=aes(x=delta_d,
                         y=acc,
                         color=d1),
             shape=3)+
  geom_line(data=summary_df,
            mapping=aes(x=delta_d,
                        y=acc,
                        color=d1),
            linetype="dashed")+
  ggthemes::theme_tufte()+
  scale_fill_manual(values = color_palette)+
  scale_color_manual(values = color_palette)

summary_kalman_sim %>% 
  ggplot(aes(x=delta,
             y=acc,
             color=d1))+
  geom_point()+
  geom_line()+
  geom_point(data=summary_df,
             mapping=aes(x=delta_d,
                         y=acc,
                         color=d1),
             shape=3)+
  geom_line(data=summary_df,
            mapping=aes(x=delta_d,
                        y=acc,
                        color=d1),
            linetype="dashed")+
  ggthemes::theme_tufte()+
  scale_fill_manual(values = color_palette)+
  scale_color_manual(values = color_palette)

summary_kalman_sim %>% 
  ggplot(aes(x=delta,
             y=acc))+
  geom_point()+
  geom_line()+
  # geom_smooth(se=F)+
  geom_point(data=summary_df,
             mapping=aes(x=delta_d,
                         y=acc),
             shape=3)+
  geom_line(data=summary_df,
            mapping=aes(x=delta_d,
                        y=acc),
            linetype="dashed")+
  ggthemes::theme_tufte()+
  scale_fill_manual(values = color_palette)+
  scale_color_manual(values = color_palette)+
  facet_wrap(~d1)
#------------------------------------------------------------------------------------#
# MODEL FITTING: M3_2
#------------------------------------------------------------------------------------#
S = length(unique(df$id))

# create stan data
stan_data = list(
  `T`  = nrow(df),
  S    = S,
  sub  = df$id,
  d1   = df$d1,
  d2   = df$d2,
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
m3_3 <-  stan("m3_deJong_3.stan",
            init=init(4),
            data=stan_data,
            chains=4,
            iter = 500,
            cores=parallel::detectCores())

# setwd("/Users/lukas/Documents/UniHeidel/Project_Discrimination/fits")
# saveRDS(m3_3,"fit_dejong_m3_3.rds")
setwd("/Users/lukas/Documents/UniHeidel/Project_Discrimination/fits")
m3_3 <- readRDS("fit_dejong_m3_3.rds")

pars <- c("mu_a","mu_ndt","mu_z0","mu_bz","mu_v0","mu_b1v","mu_g")
mcmc_pairs(m3_3, pars=pars)


# summarize predicted data with d2 averaged over delta_d
summary_pred_m3 <- pred_m3_3 %>% 
  group_by(delta_d,
           d1) %>% 
  mutate(d2_new=mean(d2)) %>% 
  ungroup() %>% 
  group_by(d1,
           d2_new) %>% 
  summarise(acc=mean(resp),
            n=length(resp)) %>% 
  mutate(d1=as.factor(d1)) %>% 
  filter(n>10)

# summarize emp data with d2 averaged over delta_d
summary_df_2 <- df %>% 
  group_by(delta_d,
           d1) %>% 
  mutate(d2_new=mean(d2)) %>% 
  ungroup() %>% 
  group_by(d1,
           d2_new) %>% 
  summarise(acc=mean(resp),
            n=length(resp)) %>% 
  mutate(d1=as.factor(d1)) %>% 
  filter(n>5)


# plot emp and pred data
summary_df_2 %>%
  ggplot(aes(x=d2_new,
             y=acc,
             color=d1))+
  geom_point(shape=3)+
  geom_line(linetype="dashed")+
geom_line(data = summary_pred_m3,
          mapping = aes(x=d2_new,
                        y=acc,
                        group=d1),
          color="maroon")

# emp log fit
df <- df %>%
  group_by(delta_d,
           d1) %>% 
  mutate(d2_new=mean(d2),
         d1=as.factor(d1))
  
log_emp <- brm(resp~d1*d2_new,
               data = df,
               family = bernoulli(),
               chains = 4,
               iter = 1000,
               cores = parallel::detectCores())

conditional_effects(log_emp)

# emp log fit
pred_m3_3 <- pred_m3_3 %>%
  group_by(delta_d,
           d1) %>% 
  mutate(d2_new=mean(d2),
         d1=as.factor(d1))

log_pred <- brm(resp~d1*d2_new,
               data = pred_m3_3,
               family = bernoulli(),
               chains = 4,
               iter = 1000,
               cores = parallel::detectCores())

conditional_effects(log_pred)


summary_pred_m3 <- pred_m3_3 %>% 
  group_by(delta_d,
           d1) %>% 
  summarise(acc=mean(resp),
            n=length(resp)) %>% 
  mutate(d1=as.factor(d1)) %>% 
  filter(n>10)

# summarize emp data with d2 averaged over delta_d
summary_df_2 <- df %>% 
  group_by(delta_d,
           d1) %>% 
  summarise(acc=mean(resp),
            n=length(resp)) %>% 
  mutate(d1=as.factor(d1)) %>% 
  filter(n>5)


# plot emp and pred data
summary_df_2 %>%
  ggplot(aes(x=delta_d,
             y=acc,
             color=d1))+
  geom_point(shape=3)+
  geom_line(linetype="dashed")+
  geom_line(data = summary_pred_m3,
            mapping = aes(x=delta_d,
                          y=acc,
                          group=d1,
                          color=d1))

# emp log fit
df <- df %>%
  group_by(delta_d,
           d1) %>% 
  mutate(d2_new=mean(d2),
         d1=as.factor(d1))