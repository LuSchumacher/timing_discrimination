library(tidyverse)
library(brms)
library(rstan)
library(magrittr)
library(lemon)
library(ggthemes)
library(latex2exp)
library(bayesplot)
library(lemon)
library(viridis)
library(bridgesampling)
library(HDInterval)
library(cowplot)
library(bayestestR)
library(tidybayes)
library(rcartocolor)
color_palette <- c("#2E5868","#B06988")

# read empirical data
setwd("/users/lukas/documents/github/timing_discrimination/data")
df <- read_csv("exp1_cond_C.csv")

# summarize empirical data
df_summary <- df %>% 
  group_by(cpos,
           cdur) %>% 
  summarise(Acc=mean(resp),
            Acc_sd=sd(resp),
            rt_median=median(RT)) %>% 
  ungroup()

df_summary$Acc <- df_summary$Acc - 1
df_summary$Acc[df_summary$cpos==1] <- 1 - df_summary$Acc[df_summary$cpos==1]
df_summary$cpos <- as.factor(df_summary$cpos)


# model comparison
setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/fits")
# Difference model
m1 <- readRDS("fit_m1_new.rds")
loo_m1 <- loo(m1)
# IRM
m2 <- readRDS("fit_m2_new.rds")
loo_m2 <- loo(m2)
# SWM
m3 <- readRDS("fit_m6_new.rds")
loo_m3 <- loo(m3)
# IRM 2x IR
m4 <- readRDS("fit_m4_exp1.rds")
loo_m4 <- loo(m4)
# IRM + SWM
m5 <- readRDS("fit_m5_new.rds")
loo_m5 <- loo(m5)
# IRM + bias
m6 <- readRDS("fit_m4_new.rds")
loo_m6 <- loo(m6)
# SWM + bias
m7 <- readRDS("fit_m7_new.rds")
loo_m7 <- loo(m7)
# IRM 2xIR
m8 <- readRDS("fit_m9_exp1_actual.rds")
loo_m8 <- loo(m8)
# IRM + SWM + bias
m9 <- readRDS("fit_m9_new.rds")
loo_m9 <- loo(m9)

# compare goodness-of-fit
comparison <- loo_compare(loo_m1,
                          loo_m2,
                          loo_m3,
                          loo_m4,
                          loo_m5,
                          loo_m6,
                          loo_m7,
                          loo_m8)

## plot model comparison
#------------------------------------------------------------------------#
# df with cols: model, elpd_diff, SE
comparison %<>% 
  as_tibble(rownames = "model") %>% 
  mutate(model=str_replace_all(model,"(.{5})", "\\1 "),
         Experiment=1)

# write.csv(comparison,"/users/lukas/documents/UniHeidel/Project_Discrimination/fits/loo_comparison_exp1.csv")
# 
# # read comparison
# comparison <- read_csv("/users/lukas/documents/UniHeidel/Project_Discrimination/fits/loo_comparison_exp1.csv")

comparison %>% 
  ggplot(aes(x = elpd_diff,
         y = model))+
  geom_point(color=carto_pal(7, "BurgYl")[7],
             size=3)+
  geom_linerange(aes(xmin=elpd_diff-se_diff*5,
                      xmax=elpd_diff+se_diff*5),
                 color="#B06988",
                 size=1)+
  geom_linerange(aes(xmin=elpd_diff-se_diff,
                      xmax=elpd_diff+se_diff),
                  color=carto_pal(7, "BurgYl")[7],
                  size=2)+
  geom_vline(xintercept = 0,
             linetype="dashed",
             color="#969696")+
  ggthemes::theme_tufte()+
  ylab("Model")+
  xlab("\nDifference in Elpd")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size =18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))

#------------------------------------------------------------------------#
# POSTERIOR PREDICTIVE CHECK: ACCURACY
#------------------------------------------------------------------------#
setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/fits")
m7 <- readRDS("fit_m7_new.rds")
pars <- c("mu_a","mu_ndt","mu_z0","mu_bz","mu_v0","mu_b1v","mu_b2v")

# extract posterior samples
mat <- as.data.frame(rstan::extract(m7, pars=pars))

# autocorrelation lag function
# mcmc_acf(m7,pars=pars)

# sample 500 parameter sets
mat$index <- 1:nrow(mat)
mat <- mat[sample(max(mat$index),100),]

# experiment settings
nTrials <- nrow(df)
d1 = df$d1-500
d2 = df$d2-500
cdur <- df$cdur
cpos <- as.numeric(df$cpos)

# simulate new data
discrimination_ddm_sim <- function(nTrials, mat){
  
  setwd("/Users/lukas/documents/UniHeidel/code")
  source("wienerProcess2.R")
  library(svMisc)
  
  data <- data.frame(matrix(ncol=7,nrow=0, dimnames=list(NULL, c("resp","rt","d1","d2","cdur","cpos","dataset"))))
  
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
                               d1=d1[i],d2=d2[i],cdur=cdur[i],cpos=cpos[i],dataset=t)
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    
  }
  
  return(data)
}

pred_m7 <- discrimination_ddm_sim(nTrials,mat)

# setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/fits")
# write.csv(pred_m7,"pred_m7.csv")
# 
# setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/fits")
# pred_m7 <- read_csv("pred_m7.csv")

# summarize within dataset
summary_pred_m7 <- pred_m7 %>% 
  group_by(dataset,
           cpos,
           cdur) %>% 
  summarise(acc=mean(resp),
            rt=median(rt))

# switcheroo
summary_pred_m7$acc[summary_pred_m7$cpos==2] <- 1 - summary_pred_m7$acc[summary_pred_m7$cpos==2]

# summarize over datasets
summary_pred_m7 <- summary_pred_m7 %>% 
  group_by(cpos,
           cdur) %>%
  summarise(acc_quantile_low=quantile(acc,probs=0.025),
            acc_quantile_high=quantile(acc,probs=0.975),
            acc_hdi_low=HDInterval::hdi(acc,credMass = 0.89)['lower'],
            acc_hdi_high=HDInterval::hdi(acc,credMass = 0.89)['upper'],
            acc=mean(acc),
            rt_quantile_low=quantile(rt,probs=0.025),
            rt_quantile_high=quantile(rt,probs=0.975),
            rt_hdi_low=HDInterval::hdi(rt,credMass = 0.89)['lower'],
            rt_hdi_high=HDInterval::hdi(rt,credMass = 0.89)['upper'],
            rt=median(rt))

summary_pred_m7$cpos <- as.factor(summary_pred_m7$cpos)

# plot posterior predictive check
tiff('/users/lukas/desktop/pp_check_acc.tiff', units="in", width=8, height=5, res=400)
summary_pred_m7 %>% ggplot(aes(x=cdur,
                                   y = acc,
                                   color=cpos))+
  geom_ribbon(aes(ymin = acc_hdi_low,
                  ymax = acc_hdi_high,
                  fill=cpos),
              alpha=0.2,
              color=NA)+
  geom_point(data=df_summary,
             mapping=aes(x=cdur,
                         y=Acc,
                         color=cpos),
             size=2)+
  geom_line(data=df_summary,
            mapping=aes(x=cdur,
                        y=Acc,
                        color=cpos),
            size=0.5)+
  scale_x_continuous(breaks=unique(summary_pred_m7$cdur),
                     labels=unique(summary_pred_m7$cdur),
                     expand = c(0.01,0.01))+
  ggthemes::theme_tufte()+
  scale_color_manual(values = color_palette)+
  scale_fill_manual(values = color_palette,
                    guide=F)+
  scale_linetype_manual(values=c("solid","dashed"))+
  ylab("Probability for c > s response")+
  xlab("\nDuration of c (ms)")+
  labs(color = "Position of c")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size =18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))

dev.off()
#------------------------------------------------------------------------#
# POSTERIOR PREDICTIVE CHECK: SUBJECT LEVEL
#------------------------------------------------------------------------#
mat <- as.data.frame(m7) %>%
  select(starts_with(c("a[","ndt[","z0[","bz[","v0[","b1v[","b2v[")))

numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
}

a <- mat %>% 
  select(`a[1]`:`a[21]`) %>% 
  pivot_longer(cols = `a[1]`:`a[21]`,
               names_to = "sub",
               values_to = "a") %>% 
  mutate(sample=rep(seq(1,nrow(mat)),each=nrow(.)/nrow(mat)),
         sub=numextract(sub)) %>% 
  select(sample,sub,a)

ndt <- mat %>% 
  select(`ndt[1]`:`ndt[21]`) %>% 
  pivot_longer(cols = `ndt[1]`:`ndt[21]`,
               names_to = "sub",
               values_to = "ndt") %>% 
  select(ndt)

z0 <- mat %>% 
  select(`z0[1]`:`z0[21]`) %>% 
  pivot_longer(cols = `z0[1]`:`z0[21]`,
               names_to = "sub",
               values_to = "z0") %>% 
  select(z0)

bz <- mat %>% 
  select(`bz[1]`:`bz[21]`) %>% 
  pivot_longer(cols = `bz[1]`:`bz[21]`,
               names_to = "sub",
               values_to = "bz") %>% 
  select(bz)

v0 <- mat %>% 
  select(`v0[1]`:`v0[21]`) %>% 
  pivot_longer(cols = `v0[1]`:`v0[21]`,
               names_to = "sub",
               values_to = "v0") %>% 
  select(v0)

b1v <- mat %>% 
  select(`bz[1]`:`b1v[21]`) %>% 
  pivot_longer(cols = `b1v[1]`:`b1v[21]`,
               names_to = "sub",
               values_to = "b1v") %>% 
  select(b1v)

b2v <- mat %>% 
  select(`b2v[1]`:`b2v[21]`) %>% 
  pivot_longer(cols = `b2v[1]`:`b2v[21]`,
               names_to = "sub",
               values_to = "b2v") %>% 
  select(b2v)


mat <- cbind(a,ndt,z0,bz,v0,b1v,b2v)

# sample posteriors
n_samples <- 100
posterior_samples <- sample(max(mat$sample),n_samples)
posteriors <- mat %>% 
  filter(sample %in% posterior_samples)

sub_sim <- function(mat){
  setwd("/Users/lukas/documents/UniHeidel/code")
  source("wienerProcess2.R")
  library(svMisc)
  
  col_names <- c("dataset","sub","d1","d2","cdur","cpos","resp", "rt")
  data <- matrix(ncol=8,
                 nrow=nrow(df)*n_samples,
                 dimnames=list(NULL, col_names))
  df_nrow <- nrow(df)
  
  x <- 1
  # iterate over posterior samples
  for (t in unique(mat$sample)){
    
    start_time <- Sys.time()
    cum_nrow <- 0
    # iterate over subjects
    for (s in 1:length(unique(mat$sub))){
      a   <- mat$a[mat$sample==t][s]
      z   <- mat$z[mat$sample==t][s]
      bz  <- mat$bz[mat$sample==t][s]
      v0  <- mat$v0[mat$sample==t][s]
      b1v <- mat$b1v[mat$sample==t][s]
      b2v <- mat$b2v[mat$sample==t][s]
      ndt <- mat$ndt[mat$sample==t][s]
      
      tmp_df <- df %>% 
        filter(sub==s)
      
      nTrials <- nrow(tmp_df)
      d1 = tmp_df$d1-500
      d2 = tmp_df$d2-500
      cdur <- tmp_df$cdur
      cpos <- as.numeric(tmp_df$cpos)

      # beta <- z + bz * d1
      # delta <- v0 + b1v * d1 + b2v * d2
      
      for (i in 1:nTrials) {
        
        delta <-  v0 + b1v * d1[i] + b2v * d2[i]
        beta <- z + bz * d1[i]
        
        diff_process <- wienerProcess2(v=delta, a=a, z=beta, ndt=ndt)
        data[((df_nrow*(x-1))+cum_nrow+i),] <- c(x, s, d1[i], d2[i],cdur[i],cpos[i],
                                                 diff_process[1],diff_process[2])
      }
      cum_nrow <- cum_nrow+nTrials
    }
    end_time <- Sys.time()
    print(end_time - start_time)
    x <- x+1
  }
  
  return(as.data.frame(data))
}

pred_sub_m7 <- sub_sim(mat = posteriors)

# setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/fits")
# write.csv(pred_sub_m7,"pred_sub_m7_exp1.csv")
# setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/fits")
# pred_sub_m7_exp1 <- read_csv("pred_sub_m7_exp1.csv")

# summarize within dataset
summary_pred_sub_m7 <- pred_sub_m7 %>% 
  group_by(dataset,
           sub,
           cpos,
           cdur) %>% 
  summarise(acc=mean(resp),
            rt=median(rt))

# switcheroo
summary_pred_sub_m7$acc[summary_pred_sub_m7$cpos==2] <- 1 - summary_pred_sub_m7$acc[summary_pred_sub_m7$cpos==2]

# summarize over datasets
summary_pred_sub_m7 <- summary_pred_sub_m7 %>% 
  group_by(sub,
           cpos,
           cdur) %>%
  summarise(acc_hdi_low=HDInterval::hdi(acc,credMass = 0.89)['lower'],
            acc_hdi_high=HDInterval::hdi(acc,credMass = 0.89)['upper'],
            acc=mean(acc),
            rt_hdi_low=HDInterval::hdi(rt,credMass = 0.89)['lower'],
            rt_hdi_high=HDInterval::hdi(rt,credMass = 0.89)['upper'],
            rt=median(rt))

summary_pred_sub_m7$cpos <- as.factor(summary_pred_sub_m7$cpos)

summary <- df %>%
  group_by(sub,
           cdur,
           cpos) %>% 
  mutate(resp=2-resp,
         cpos=as.factor(cpos)) %>% 
  summarise(prob_c = mean(resp),
            rt=median(RT))

summary$prob_c[summary$cpos==2] <- 1 - summary$prob_c[summary$cpos==2]
summary$cpos <- as.factor(summary$cpos)

tiff('/users/lukas/desktop/pp_check_sub_acc_exp1.tiff', units="in", width=9, height=10, res=400)
summary_pred_sub_m7 %>% 
  ggplot(aes(x=cdur,
             y=acc,
             color=cpos,
             fill=cpos)) +
  geom_ribbon(aes(ymin = acc_hdi_low,
                  ymax = acc_hdi_high,
                  fill=cpos),
              alpha=0.2,
              color=NA)+
  geom_point(data=summary,
             mapping=aes(x=cdur,
                         y=prob_c,
                         color=cpos))+
  geom_line(data=summary,
            mapping=aes(x=cdur,
                        y=prob_c,
                        color=cpos),
            alpha=0.5)+
  facet_wrap(~sub,nrow=7)+
  ggthemes::theme_tufte()+
  scale_color_manual(values = color_palette)+
  scale_fill_manual(values = color_palette,
                    guide=F)+
  ylab("Probability for c > s response")+
  xlab("\nDuration of c")+
  labs(color = "Position of c")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 12),
        text=element_text(size = 16),
        plot.title = element_text(size =18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))

dev.off()
#------------------------------------------------------------------------#
# POSTERIOR PREDICTIVE CHECK: DIFFERENCE LIMEN
#------------------------------------------------------------------------#
# function for transform probability to logit
prob2logit <- function(prob){
  logit <- log(prob/(1-prob))
  return(logit)
}

pred_sub_m7$DL_1 <- NA
pred_sub_m7$DL_2 <- NA
pred_sub_m7$PSE_1 <- NA
pred_sub_m7$PSE_2 <- NA

for (s in unique(pred_sub_m7$sub)) {
  tmp <- pred_sub_m7 %>% 
    filter(sub==s) %>% 
    mutate(cpos=as.factor(cpos),
           cdur=(cdur-500))
  
  # tmp$resp[tmp$cpos==2] <- 1 - tmp$resp[tmp$cpos==2]
  
  # logstic regression for psychometric curve
  logReg <- tmp %>%
    brm(formula = resp~cdur*cpos,
        family = bernoulli(),
        chains = 4,
        iter = 1000,
        cores=parallel::detectCores(),
        control = list(adapt_delta=0.95))
  
  params <- as.numeric(as.data.frame(logReg) %>% colMeans())
  
  # DL and PSE for cpos = 1
  c25_1 <- (prob2logit(0.25)-params[1])/params[2]
  pred_sub_m7$PSE_1[pred_sub_m7$sub==s] <- (prob2logit(0.50)-params[1])/params[2]
  c75_1 <- (prob2logit(0.75)-params[1])/params[2]
  pred_sub_m7$DL_1[pred_sub_m7$sub==s] <- (c75_1 - c25_1)/2
  # DL and PSE for cpos = 2
  c25_2 <- (prob2logit(0.25)-params[1]-params[3])/(params[2]+params[4])
  pred_sub_m7$PSE_2[pred_sub_m7$sub==s] <- (prob2logit(0.50)-params[1]-params[3])/(params[2]+params[4])
  c75_2 <- (prob2logit(0.75)-params[1]-params[3])/(params[2]+params[4])
  pred_sub_m7$DL_2[pred_sub_m7$sub==s] <- (c75_2 - c25_2)/2
}

summary_DL_PSE <- pred_sub_m7 %>% 
  group_by(sub) %>% 
  summarise(DL_1=DL_1[1],
            DL_2=abs(DL_2[1]),
            PSE_1=PSE_1[1],
            PSE_2=PSE_2[1]) %>% 
  mutate(TypeB=DL_2-DL_1,
         TOE=abs(PSE_2)-abs(PSE_1))

summary_emp_DL_PSE <- df %>% 
  group_by(sub) %>% 
  summarise(DL_1=DL_1[1],
            DL_2=DL_2[1],
            PSE_1=PSE_1[1],
            PSE_2=PSE_2[1]) %>% 
  mutate(DL_diff=DL_2-DL_1,
         PSE_dff=abs(PSE_2)-abs(PSE_1))


# Type B error
sumsum <- data_frame(emp=summary_emp_DL_PSE$DL_diff,
                     pred=summary_DL_PSE$TypeB)

tiff('/users/lukas/desktop/pp_DL.tiff', units="in", width=6, height=5, res=400)
sumsum %>%
  ggplot(aes(x=emp,
             y=pred))+
  geom_point()+
  geom_abline(intercept = 0,
              slope = 1,
              linetype="dashed",
              color="#969696")+
  scale_x_continuous(limits = c(-300,100))+
  scale_y_continuous(limits = c(-300,100))+
  geom_text(aes(x=0,
                y=-200),
            label=paste("r =",as.character(round(cor(sumsum$emp,sumsum$pred),digits=2))),
            family = "serif",
            size=5)+
  ggthemes::theme_tufte()+
  ylab("Predicted type B error")+
  xlab("\nEmpirical type B error")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size =18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))

dev.off()

# Type A error
sumsum <- data_frame(emp=summary_emp_DL_PSE$PSE_dff,
                     pred=summary_DL_PSE$TOE)

cor(sumsum$emp,sumsum$pred)

tiff('/users/lukas/desktop/pp_TOE.tiff', units="in", width=6, height=5, res=400)
sumsum %>%
  ggplot(aes(x=emp,
             y=pred))+
  geom_point()+
  geom_abline(intercept = 0,
              slope = 1,
              linetype="dashed",
              color="#969696")+
  scale_x_continuous(limits = c(-110,110))+
  scale_y_continuous(limits = c(-110,110))+
  geom_text(aes(x=50,
                y=-50),
            label=paste("r =",as.character(round(cor(sumsum$emp,sumsum$pred),digits=2))),
            family = "serif",
            size=5)+
  ggthemes::theme_tufte()+
  ylab("Predicted type A error")+
  xlab("\nEmpirical type A error")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size =18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))
dev.off()
#------------------------------------------------------------------------#
# POSTERIOR PREDICTIVE CHECK: WEIGHT DIFFERENCE
#------------------------------------------------------------------------#
# get individual point estimates of b1v and b2v parameter
individual_b1v <- summary(m7,pars=c("b1v[1]","b1v[2]","b1v[3]","b1v[4]","b1v[5]","b1v[6]","b1v[7]","b1v[8]","b1v[9]","b1v[10]","b1v[11]","b1v[12]","b1v[13]","b1v[14]","b1v[15]","b1v[16]","b1v[17]","b1v[18]","b1v[19]","b1v[20]","b1v[21]"))$summary[,"mean"]
individual_b2v <- summary(m7,pars=c("b2v[1]","b2v[2]","b2v[3]","b2v[4]","b2v[5]","b2v[6]","b2v[7]","b2v[8]","b2v[9]","b2v[10]","b2v[11]","b2v[12]","b2v[13]","b2v[14]","b2v[15]","b2v[16]","b2v[17]","b2v[18]","b2v[19]","b2v[20]","b2v[21]"))$summary[,"mean"]

summary_DL_PSE$individual_b1v <- individual_b1v
summary_DL_PSE$individual_b2v <- individual_b2v

# summary_DL_PSE$weightdiff <- abs(summary_DL_PSE$individual_b1v) - abs(summary_DL_PSE$individual_b2v)

summary_DL_PSE$weightdiff_percent <- (abs(summary_DL_PSE$individual_b1v) - abs(summary_DL_PSE$individual_b2v))/(abs(summary_DL_PSE$individual_b1v) + abs(summary_DL_PSE$individual_b2v))

# weights difference on DL
sumsum <- data_frame(emp=summary_emp_DL_PSE$DL_diff,
                     weight_diff=summary_DL_PSE$weightdiff_percent)

p1 <- sumsum %>%
  ggplot(aes(x=emp,
             y=weight_diff))+
  geom_point(size=2,
             alpha=0.9)+
  scale_x_continuous(expand = c(0.01,1))+
  # scale_y_continuous(limits = c(-110,110))+
  geom_smooth(method='lm',se=T,
              color="black",
              fill="gray",
              linetype="dashed",
              size=0.5)+
  geom_text(aes(x=-50,
                y=-1),
            label=paste("r =",as.character(round(cor(sumsum$emp,sumsum$weight_diff),digits=2))),
            family = "serif",
            size=5)+
  ggthemes::theme_tufte()+
  ylab("Weights difference")+
  xlab("\nEmpirical type B error")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size =18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  scale_x_continuous(breaks = c(-250,-200,-150,-100,-50,0),
                     expand = c(0.01,0.01))+
  scale_y_continuous(breaks = seq(-1.2,0,0.2),
                     limits = c(-1.3,0.1))
p1

# weights difference on PSE
sumsum <- data_frame(emp=summary_emp_DL_PSE$PSE_dff,
                     weight_diff=summary_DL_PSE$weightdiff_percent)

p2 <- sumsum %>%
  ggplot(aes(x=emp,
             y=weight_diff))+
  geom_point(size=2,
             alpha=0.9)+
  scale_x_continuous(expand = c(0.01,1))+
  # scale_y_continuous(limits = c(-110,110))+
  geom_smooth(method='lm',se=T,
              color="black",
              fill="gray",
              linetype="dashed",
              size=0.5)+
  geom_text(aes(x=0,
                y=-1),
            label=paste("r =",as.character(round(cor(sumsum$emp,sumsum$weight_diff),digits=2))),
            family = "serif",
            size=5)+
  ggthemes::theme_tufte()+
  ylab("Weights difference")+
  xlab("\nEmpirical type A error")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size =18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  scale_x_continuous(expand = c(0.01,0.01))+
  scale_y_continuous(breaks = seq(-1.2,0,0.2),
                     limits = c(-1.3,0.1))

p2

tiff('/users/lukas/desktop/pp_check_DL.tiff', units="in", width=9, height=5, res=400)
plot_grid(p2,p1,
          labels = "AUTO",
          label_fontfamily= "serif",
          label_size = 18)
dev.off()
#------------------------------------------------------------------------#
# POSTERIOR PREDICTIVE CHECK: RESPONSE TIMES
#------------------------------------------------------------------------#
setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/fits")
pred_sub_m7 <- read_csv("pred_sub_m7.csv")
pred_m7 <- read_csv("pred_m7.csv")


tiff('/users/lukas/desktop/pp_check_rt.tiff', units="in", width=8, height=5, res=400)
summary_pred_m7 %>% ggplot(aes(x=cdur,
                               y = rt,
                               color=cpos))+
  geom_ribbon(aes(ymin = rt_hdi_low,
                  ymax = rt_hdi_high,
                  fill=cpos),
              alpha=0.2,
              color=NA)+
  # geom_line(size=0.8)+
  # geom_point(size=1)
  geom_point(data=df_summary,
             mapping=aes(x=cdur,
                         y=rt_median/1000,
                         color=cpos),
             size=2)+
  geom_line(data=df_summary,
            mapping=aes(x=cdur,
                        y=rt_median/1000,
                        color=cpos),
            size=0.5,
            linetype="dashed")+
  scale_x_continuous(breaks=unique(summary_pred_m7$cdur),
                     labels=unique(summary_pred_m7$cdur),
                     expand = c(0.01,0.01))+
  scale_y_continuous(limits = c(0.1,0.9))+
  ggthemes::theme_tufte()+
  scale_color_manual(values = color_palette)+
  scale_fill_manual(values = color_palette,
                    guide=F)+
  ylab("Median response time (s)")+
  xlab("\nDuration of c")+
  labs(color = "Position of c")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size =18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))

dev.off()


tiff('/users/lukas/desktop/pp_check_sub_rt.tiff', units="in", width=16, height=8, res=400)
summary_pred_sub_m7 %>% 
  ggplot(aes(x=cdur,
             y=rt,
             color=cpos,
             fill=cpos)) +
  # geom_point(shape=4)+
  # geom_line()+
  geom_ribbon(aes(ymin = rt_quantile_low,
                  ymax = rt_quantile_high,
                  fill=cpos),
              alpha=0.2,
              color=NA)+
  geom_point(data=summary,
             mapping=aes(x=cdur,
                         y=rt/1000,
                         color=cpos))+
  geom_line(data=summary,
            mapping=aes(x=cdur,
                        y=rt/1000,
                        color=cpos),
            alpha=0.5,
            linetype="dashed")+
  facet_wrap(~sub,nrow=4)+
  ggthemes::theme_tufte()+
  scale_color_manual(values = color_palette)+
  scale_fill_manual(values = color_palette,
                    guide=F)+
  ylab("Median response time")+
  xlab("\nDuration of c")+
  labs(color = "Position of c")+
  theme(axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        text=element_text(size = 16),
        plot.title = element_text(size =18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))

dev.off()
#------------------------------------------------------------------------#
# POSTERIOR PREDICTIVE CHECK: RT QUANTILE PLOT
#------------------------------------------------------------------------#
df_rt_quantiles <- df %>% 
  filter(cdur != 500) %>% 
  rename(rt=RT)

qfun <- function(x, quantiles = c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  enframe(quantile(x$rt, probs = quantiles, na.rm = TRUE),
          name = "quantile",
          value = "rt")
}

df_rt_quantiles %<>% 
  drop_na() %>% 
  droplevels() %>%
  group_by(cpos, correct, cdur) %>%
  nest() %>%
  mutate(data = purrr::map(data, qfun)) %>%
  unnest() %>% 
  mutate(quantile = as_factor(quantile)) %>% 
  mutate(rt=rt/1000)

sumsum <- df %>% 
  drop_na() %>% 
  droplevels() %>%
  group_by(cpos, correct, cdur) %>%
  summarise(n=length(RT))
  
# have to flip the responses
pred_m7_rt_quantiles <- pred_m7 %>% 
  mutate(resp=ifelse(resp==1,0,1))

# calculate errors
pred_m7_rt_quantiles$correct <- NA
for (i in 1:nrow(pred_m7_rt_quantiles)) {
  
  if(pred_m7_rt_quantiles$resp[i]==0 && pred_m7_rt_quantiles$cdur[i]<500){
    pred_m7_rt_quantiles$correct[i] <- 1
  }else if(pred_m7_rt_quantiles$resp[i]==1 && pred_m7_rt_quantiles$cdur[i]>500){
    pred_m7_rt_quantiles$correct[i] <- 1
  }else{
    pred_m7_rt_quantiles$correct[i] <- 0
  }
}

# remove trials cdur=500, because the response is always wrong
pred_m7_rt_quantiles %<>% 
  filter(cdur != 500)

# calculate RT quantiles
pred_m7_rt_quantiles %<>% 
  drop_na() %>% 
  droplevels() %>%
  group_by(dataset, cpos, correct, cdur) %>%
  nest() %>%
  mutate(data = purrr::map(data, qfun)) %>%
  unnest() %>% 
  mutate(quantile = as_factor(quantile)) %>% 
  ungroup() %>% 
  group_by(cpos,
           cdur,
           correct,
           quantile) %>% 
  summarise(rt_median=median(rt),
            quantile_low=HDInterval::hdi(rt,credMass=0.89)["lower"],
            quantile_high=HDInterval::hdi(rt,credMass=0.89)["upper"])

pred_m7_rt_quantiles %>%
  ggplot(aes(x = as_factor(cdur), y = rt_median,
             group = quantile,
             color = quantile)) +
  # geom_line(alpha = 1, size = 1) +
  # geom_point(alpha = 1, size = 4, shape=1) +
  geom_ribbon(aes(ymin = quantile_low,
                  ymax = quantile_high,
                  group = quantile,
                  fill = quantile),
              alpha=0.2,
              color=NA)+
  geom_point(data = df_rt_quantiles,
             mapping = aes(x = as_factor(cdur),
                           y = rt,
                           group = quantile,
                           color = quantile))+
  geom_line(data = df_rt_quantiles,
            mapping = aes(x = as_factor(cdur),
                          y = rt,
                          group = quantile,
                          color = quantile))+
  facet_grid(correct~cpos,
             labeller = label_both) +
  scale_color_viridis(direction=-1,discrete = TRUE, option = "E") +
  scale_fill_viridis(direction=-1,discrete = TRUE, option = "E",guide=F) +
  xlab("") +
  ylab("Response time (s)") +
  ggtitle("") +
  ggthemes::theme_tufte()+
  ylab("Median response time")+
  xlab("\nDuration of c")+
  labs(color = "Quantile")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size =18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))

##------------##
## Second Attempt
##------------##
df_rt_quantiles <- df %>% 
  rename(rt=RT) %>% 
  droplevels() %>%
  group_by(cpos, cdur) %>%
  nest() %>%
  mutate(data = purrr::map(data, qfun)) %>%
  unnest() %>% 
  mutate(quantile = as_factor(quantile)) %>% 
  mutate(rt=rt/1000) %>% 
  rename(`Position of c`=cpos)

sumsum_pred_m7 <- pred_m7 %>% 
  mutate(resp=ifelse(resp==1,0,1)) %>% 
  droplevels() %>%
  group_by(dataset, cpos, cdur) %>%
  nest() %>%
  mutate(data = purrr::map(data, qfun)) %>%
  unnest() %>% 
  mutate(quantile = as_factor(quantile)) %>% 
  ungroup() %>% 
  group_by(cpos,
           cdur,
           quantile) %>% 
  summarise(rt_median=median(rt),
            quantile_low=HDInterval::hdi(rt,credMass=0.89)["lower"],
            quantile_high=HDInterval::hdi(rt,credMass=0.89)["upper"]) %>% 
  rename(`Position of c`=cpos)

tiff('/users/lukas/desktop/pp_check_rt.tiff', units="in", width=9, height=5, res=400)
sumsum_pred_m7 %>%
  ggplot(aes(x = as_factor(cdur), y = rt_median,
             group = quantile,
             color = quantile)) +
  # geom_line(alpha = 1, size = 1) +
  # geom_point(alpha = 1, size = 4, shape=1) +
  geom_ribbon(aes(ymin = quantile_low,
                  ymax = quantile_high,
                  group = quantile,
                  fill = quantile),
              alpha=0.2,
              color=NA)+
  geom_point(data = df_rt_quantiles,
             mapping = aes(x = as_factor(cdur),
                           y = rt,
                           group = quantile,
                           color = quantile))+
  geom_line(data = df_rt_quantiles,
            mapping = aes(x = as_factor(cdur),
                          y = rt,
                          group = quantile,
                          color = quantile))+
  facet_grid(~`Position of c`,
             labeller = label_both) +
  scale_color_viridis(direction=-1,discrete = TRUE, option = "E") +
  scale_fill_viridis(direction=-1,discrete = TRUE, option = "E",guide=F) +
  ggtitle("") +
  ggthemes::theme_tufte()+
  ylab("Response time (s)")+
  xlab("\nDuration of c (ms)")+
  labs(color = "Quantile")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14,
                                   angle = 45, vjust = 0.5, hjust=1),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size =18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))
dev.off()