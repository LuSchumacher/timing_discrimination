library(tidyverse)
library(magrittr)
library(rstan)
source("C:/Users/Lukas.Schumacher/Documents/GitHub/timing_discrimination/code/wienerProcess2.R")
`%notin%` = Negate(`%in%`)

qfun <- function(x, quantiles = c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  enframe(quantile(x$rt, probs = quantiles, na.rm = TRUE),
          name = "quantile",
          value = "rt")}

qfun_new <- function(x, quantiles = c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  enframe(quantile(x$rt_new, probs = quantiles, na.rm = TRUE),
          name = "quantile",
          value = "rt")}

color_palette <- c("#2E5868","#B06988")
#---------------------------------------------------------------------------#
# DATA PREPARATION
#---------------------------------------------------------------------------#
setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/Dyjas Bausenhart Ulrich 2012/Experiment 2")
files <- list.files(pattern = "*D.txt")
emp_data <- NULL
for (f in files) {
  sub <- as.numeric(substr(f, 10, nchar(f) - 5))
  if (sub == 99) next
  data <- read.table(f, skip = 30)
  colnames(data) <- c("trl", "cpos", "cdur", "resp", "RT", "type")
  data$d1 <- ifelse(data$cpos == 2, 500, data$cdur)
  data$d2 <- ifelse(data$cpos == 2, data$cdur, 500)
  emp_data <- rbind(emp_data, cbind(sub = sub, data)
  )
}

# data wrangling
emp_data %<>%
  rename(rt = RT) %>% 
  mutate(sub = as.numeric(sub),
         rt = rt / 1000) %>% 
  filter(rt != 0.0,
         rt <= 5.0,
         resp > 0) %>%
  arrange(sub)

# exclude non-cooperative subjects
emp_data %<>%
  filter(sub %notin% c(10, 19))

# reassign subject id
nr <- 1
for (i in 2:length(emp_data$sub)){
  if (emp_data$sub[i] != emp_data$sub[i - 1]){
    nr <- nr + 1
    emp_data$sub[emp_data$sub == emp_data$sub[i]] <- nr
  }
}

# evaluate correct response
emp_data$correct <- NA
for (i in 1:nrow(emp_data)) {
  if(emp_data$cpos[i] == 1 & emp_data$cdur[i] > 500){
    if(emp_data$resp[i] == 1){
      emp_data$correct[i] <- 1
    }else{
      emp_data$correct[i] <- 0
    }
  }else if(emp_data$cpos[i] == 1 & emp_data$cdur[i] < 500){
    if(emp_data$resp[i] == 2){
      emp_data$correct[i] <- 1
    }else{
      emp_data$correct[i] <- 0
    }
  }else if(emp_data$cpos[i] == 2 & emp_data$cdur[i] > 500){
    if(emp_data$resp[i] == 2){
      emp_data$correct[i] <- 1
    }else{
      emp_data$correct[i] <- 0
    }
  }else if(emp_data$cdur[i] == 500){
    emp_data$correct[i] <- 0
  }
  else{
    if(emp_data$resp[i] == 1){
      emp_data$correct[i] <- 1
    }else{
      emp_data$correct[i] <- 0
    }
  }
}

# rt manipulation
# whenever d2 was shorter than d1, then the rt decreased by the difference d1 - d2
# emp_data$rt_new <- NA
# for (i in 1:nrow(emp_data)){
#   if (emp_data$d2[i] < emp_data$d1[i]){
#     emp_data$rt_new[i] <- emp_data$rt[i] + ((emp_data$d2[i] - emp_data$d1[i]) / 1000)
#   }
#   else {
#     emp_data$rt_new[i] <- emp_data$rt[i]
#   }
# }

# # rt manipulation
# # whenever d2 was longer than d1, then the rt increased by the difference d2 - d1
# emp_data$rt_new <- NA
# for (i in 1:nrow(emp_data)){
#   if (emp_data$d2[i] > emp_data$d1[i]){
#     emp_data$rt_new[i] <- emp_data$rt[i] + ((emp_data$d2[i] - emp_data$d1[i]) / 1000)
#   }
#   else {
#     emp_data$rt_new[i] <- emp_data$rt[i]
#   }
# }

# emp_data$rt_new <- NA
# for (i in 1:nrow(emp_data)){
#   if (emp_data$d2[i] > emp_data$d1[i]){
#     emp_data$rt_new[i] <- emp_data$rt[i] + ((emp_data$d2[i] - emp_data$d1[i]) / 1000)
#   }
#   # else if (emp_data$d2[i] < emp_data$d1[i]){
#   #   emp_data$RT[i] <- emp_data$rt[i] + ((emp_data$d2[i] - emp_data$d1[i]) / 1000)
#   # }
#   else {
#     emp_data$rt_new[i] <- emp_data$rt[i]
#   }
# }

performance_summary <- emp_data %>% 
  group_by(sub) %>% 
  summarise(accuracy = mean(correct),
            median_rt = median(rt),
            # median_rt_new = median(rt_new),
            n_trials = length(rt))

# simulate random monkey
accuracy_data <- rbinom(1e4, 655, 0.5)
hist(accuracy_data/655, breaks=100)
binom.test(280, 521 , 0.5, alternative = "greater")

# check accuracy of fast guesses
lower_rt_cutoff <- 0.15
fast_guesses_summary <- emp_data %>% 
  filter(cdur != 500,
         rt <= lower_rt_cutoff) %>% 
  group_by(sub) %>% 
  summarise(accuracy = round(mean(correct), 3),
            mean_rt = round(mean(rt), 3),
            n = length(correct))

# check accuracy of delayed start up's
upper_rt_cutoff <- 2
delayed_start_ups_summary <- emp_data %>% 
  filter(cdur != 500,
         rt >= upper_rt_cutoff) %>% 
  group_by(sub) %>% 
  summarise(accuracy = round(mean(correct), 3),
            mean_rt = round(mean(rt), 3),
            max_rt = max(rt),
            n = length(correct))

# EWMA method
# Define constants
lambda   <- 0.01  # weight param (how many previous data points)
c_0      <- 0.5   # in-control mean (expected average performance for fast guess)
sigma_0  <- 0.5   # in-control sd
L        <- 1.5   # sensitivity (1.5 is a relative low value)

# Compute EMWA
data <- data.frame()
for (s in unique(emp_data$sub)) {
  
  tmp <- emp_data %>% 
    filter(sub==s) %>% 
    arrange(rt)
  
  tmp$c_s <- NA
  tmp$inCm <- "empty"
  tmp$UCL <- NA
  
  c_before <- c_0
  
  # Compute EWMA for each data point
  for (i in 1:length(tmp$rt)) {
    if(tmp$cdur[i]==500) next
    tmp$c_s[i] <- (lambda*tmp$correct[i]) + ((1-lambda)*c_before)
    tmp$UCL[i] <- c_0 + (L*sigma_0)*sqrt((lambda/(2-lambda))*(1-(1-lambda)^(2*i)))
    tmp$inCm[i] <- tmp$c_s[i] < tmp$UCL[i]
    c_before <- tmp$c_s[i]
  }
  
  data <- rbind(data,tmp)
}

# plot EWMA control chart
data %>%
  filter(inCm != "empty") %>% 
  mutate(inCm=as.logical(inCm)) %>% 
  ggplot(aes(x=rt,
             y=c_s,
             color=(inCm>0)))+
  geom_line(aes(group=2),
            size=0.3)+
  geom_hline(yintercept = 0.5,
             linetype="dotted")+
  geom_line(aes(x=rt,
                y=UCL,
                color="black"))+
  scale_color_manual(values=c("#453b44","#39C8C6","#630C3A"))+
  xlab("\nReaction time (s)")+
  ylab(latex2exp::TeX("State of the system ($c_s$)"))+
  ggtitle("EWMA control chart")+
  ggthemes::theme_tufte(base_family = "GillSans")+
  theme(legend.position="none",
        axis.line = element_line(size = .5, color = "#444444"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size =18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  # scale_x_continuous(limits = c(0,2.5))+
  facet_wrap(~sub)

data_new <- NULL
trash <- NULL
data$inCm[data$inCm == "empty"] <- FALSE
window <- 5
for (sub in unique(data$sub)){
  tmp <- data[data$sub == sub, ]
  for (i in 1:nrow(tmp)){
    if (!any(as.logical(tmp$inCm[i:(i + window)]))){
      data_new <- rbind(data_new, tmp[i:nrow(tmp), ])
      trash <- rbind(trash, tmp[1:i, ])
      print(i)
      break
    }
  }
}

data <- data_new

# summarize emp data
emp_rt_quantiles <- data %>% 
  filter(rt <= 2.5) %>%
  droplevels() %>%
  group_by(cpos,
           cdur) %>%
  nest() %>%
  mutate(data = purrr::map(data, qfun)) %>%
  unnest(cols = c(data)) %>% 
  mutate(quantile = as_factor(quantile),
         cdur = as_factor(cdur)) %>% 
  rename(`Position of c` = cpos)

# plot rt quantiles
emp_rt_quantiles %>%
  ggplot(aes(x = cdur,
             y = rt,
             group = quantile,
             color = quantile)) +
  geom_point()+
  geom_line()+
  facet_wrap(~`Position of c`,
             labeller = label_both) +
  scale_y_continuous(breaks = seq(0, 2 , by = 0.2))+
  viridis::scale_color_viridis(direction = -1,discrete = TRUE, option = "E") +
  viridis::scale_fill_viridis(direction = -1,discrete = TRUE, option = "E", guide="none") +
  ggtitle("") +
  ggthemes::theme_tufte()+
  ylab("Response time (s)")+
  xlab("\nDuration of c (ms)")+
  labs(color = "Quantile")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14,
                                   angle = 45,
                                   vjust = 0.5,
                                   hjust = 1),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size = 18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))

performance_summary <- data_new %>% 
  group_by(sub) %>% 
  summarise(accuracy = mean(correct),
            median_rt = median(rt),
            # median_rt_new = median(rt_new),
            n_trials = length(rt))

# check accuracy of fast guesses
lower_rt_cutoff <- 0.15
fast_guesses_summary <- data %>% 
  filter(cdur != 500,
         rt <= lower_rt_cutoff) %>% 
  group_by(sub) %>% 
  summarise(accuracy = round(mean(correct), 3),
            mean_rt = round(mean(rt), 3),
            n = length(correct))

# check accuracy of delayed start up's
upper_rt_cutoff <- 2.5
delayed_start_ups_summary <- data %>% 
  filter(cdur != 500,
         rt >= upper_rt_cutoff) %>% 
  group_by(sub) %>% 
  summarise(accuracy = round(mean(correct), 3),
            mean_rt = round(mean(rt), 3),
            max_rt = max(rt),
            n = length(correct))

data %<>%
  filter(rt >= 0.15,
         rt <= 2.0)


write_csv(data, "/users/lukas/documents/github/timing_discrimination/data/data_exp2_cond_D_revision_lastTry.csv")


#---------------------------------------------------------------------------#
# MODEL FITTING
#---------------------------------------------------------------------------#
df <- read_csv("C:/Users/Lukas.Schumacher/Documents/GitHub/timing_discrimination/data/data_exp2_cond_D_revision_experimental.csv")

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
  d1     = (df$d1 - 500) / 100,
  d2     = (df$d2 - 500) / 100,
  resp   = df$resp,
  rt     = df$rt_new
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
    L[[c]]$mu_b1v = 0.0
    L[[c]]$mu_b2v = 0.0
    
    L[[c]]$sd_a   = 0.001
    L[[c]]$sd_ndt = 0.001
    L[[c]]$sd_z0  = 0.001
    L[[c]]$sd_bz  = 0.001
    L[[c]]$sd_v0  = 0.001
    L[[c]]$sd_b1v = 0.001
    L[[c]]$sd_b2v = 0.001   
    
    L[[c]]$a   = runif(S,0.5,2)
    L[[c]]$ndt = runif(S,0.08,0.1)
    L[[c]]$z0  = rnorm(S,0.5,0.1)
    L[[c]]$bz  = rep(0.0,S)
    L[[c]]$v0  = runif(S,0.5,2)
    L[[c]]$b1v = rep(0.0,S)
    L[[c]]$b2v = rep(0.0,S)

  }
  return (L)
}

fit_m6 <-  stan("C:/Users/Lukas.Schumacher/Documents/GitHub/timing_discrimination/models/model_6_new_prior.stan",
                init=init(4),
                data=stan_data,
                chains=4,
                iter = 2000,
                cores=parallel::detectCores())

fit_m6

# saveRDS(fit_m6, "C:/Users/Lukas.Schumacher/Documents/GitHub/timing_discrimination/code/fit_exp2_cond_D_revision_experimental.rds")

#---------------------------------------------------------------------------#
# PP CHECK
#---------------------------------------------------------------------------#
df <- read_csv("C:/Users/Lukas.Schumacher/Documents/GitHub/timing_discrimination/data/data_exp2_cond_D_revision_experimental.csv")
fit <- readRDS("C:/Users/Lukas.Schumacher/Documents/GitHub/timing_discrimination/code/fit_exp2_cond_D_revision_experimental.rds")

# extract posterior distributions
mat <- as.data.frame(fit) %>%
  select(starts_with(c("a[", "ndt[", "z0[", "bz[", "v0[", "b1v[", "b2v[")))

numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
}

a <- mat %>% 
  select(`a[1]`:`a[22]`) %>% 
  pivot_longer(cols = `a[1]`:`a[22]`,
               names_to = "sub",
               values_to = "a") %>% 
  mutate(sample=rep(seq(1,nrow(mat)), each=nrow(.) / nrow(mat)),
         sub=numextract(sub)) %>% 
  select(sample,sub,a)

ndt <- mat %>% 
  select(`ndt[1]`:`ndt[22]`) %>% 
  pivot_longer(cols = `ndt[1]`:`ndt[22]`,
               names_to = "sub",
               values_to = "ndt") %>% 
  select(ndt)

z0 <- mat %>% 
  select(`z0[1]`:`z0[22]`) %>% 
  pivot_longer(cols = `z0[1]`:`z0[22]`,
               names_to = "sub",
               values_to = "z0") %>% 
  select(z0)

bz <- mat %>% 
  select(`bz[1]`:`bz[22]`) %>% 
  pivot_longer(cols = `bz[1]`:`bz[22]`,
               names_to = "sub",
               values_to = "bz") %>% 
  select(bz)

v0 <- mat %>% 
  select(`v0[1]`:`v0[22]`) %>% 
  pivot_longer(cols = `v0[1]`:`v0[22]`,
               names_to = "sub",
               values_to = "v0") %>% 
  select(v0)

b1v <- mat %>% 
  select(`b1v[1]`:`b1v[22]`) %>% 
  pivot_longer(cols = `b1v[1]`:`b1v[22]`,
               names_to = "sub",
               values_to = "b1v") %>% 
  select(b1v)

b2v <- mat %>% 
  select(`b2v[1]`:`b2v[22]`) %>% 
  pivot_longer(cols = `b2v[1]`:`b2v[22]`,
               names_to = "sub",
               values_to = "b2v") %>% 
  select(b2v)

mat <- cbind(a, ndt, z0, bz, v0, b1v, b2v)
mat %<>%
  mutate(sub=as.numeric(sub))


# SIMULATION
# sample posteriors
n_samples <- 500
posterior_samples <- sample(max(mat$sample), n_samples)
posteriors <- mat %>% 
  filter(sample %in% posterior_samples)

n_sub <- max(mat$sub)

col_names <- c("sim" ,"sub", "d1", "d2", "cdur", "cpos", "resp", "rt")
pred_data <- matrix(ncol=8,
                    nrow=nrow(df) * n_samples,
                    dimnames=list(NULL, col_names))
df_nrow <- nrow(df)


x <- 0
for (sim in 1:n_samples){
  start_time <- Sys.time()
  cum_nrow <- 0
  
  for (s in unique(mat$sub)){
    a   <- mat$a[mat$sample == sim][s]
    z0   <- mat$z0[mat$sample == sim][s]
    bz   <- mat$bz[mat$sample == sim][s]
    v0  <- mat$v0[mat$sample == sim][s]
    b1v <- mat$b1v[mat$sample == sim][s]
    b2v <- mat$b2v[mat$sample == sim][s]
    ndt <- mat$ndt[mat$sample == sim][s]
    
    tmp_df <- df %>% 
      filter(sub == s)
    
    n_trials <- nrow(tmp_df)
    d1 <- (tmp_df$d1 - 500) / 100  
    d2 <- (tmp_df$d2 - 500) / 100
    cdur <- tmp_df$cdur
    cpos <- as.numeric(tmp_df$cpos)
    
    for (i in 1:n_trials) {
      delta <- v0 + b1v * d1[i] + b2v * d2[i]
      beta <-  z0 + bz * d1[i]
      diff_process <- wienerProcess2(v=delta, a=a, z=beta, ndt=ndt)
      pred_data[(df_nrow*x + cum_nrow + i), ] <- c(x, s, d1[i], d2[i], cdur[i], cpos[i],
                                                   diff_process[1], diff_process[2])
      
    }
    cum_nrow <- cum_nrow + n_trials
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  x <- x + 1
}

# pred_data_df <- as.data.frame(pred_data)
# write_csv(pred_data_df, 'C:/Users/Lukas.Schumacher/Documents/GitHub/timing_discrimination/data/pp_check_exp2_cond_D_revision_experimental.csv')


#---------------------------------------------------------------------------#
# SUMMARIZING AND PLOTTING
#---------------------------------------------------------------------------#
emp_data <- read_csv("C:/Users/Lukas.Schumacher/Documents/GitHub/timing_discrimination/data/data_exp2_cond_D_revision_experimental.csv")
pred_data <- read_csv('C:/Users/Lukas.Schumacher/Documents/GitHub/timing_discrimination/data/pp_check_exp2_cond_D_revision_experimental.csv')


# individual rt quantile plot
#---------------------------------------------------------------------------#
# summarize emp data
emp_rt_quantiles <- emp_data %>% 
  droplevels() %>%
  group_by(sub,
           cpos,
           cdur) %>%
  nest() %>%
  mutate(data = purrr::map(data, qfun_new)) %>%
  unnest(cols = c(data)) %>% 
  mutate(quantile = as_factor(quantile)) %>% 
  rename(`Position of c` = cpos)

pred_rt_quantiles <- pred_data %>% 
  mutate(resp = ifelse(resp == 1, 0, 1)) %>% 
  droplevels() %>%
  group_by(sim,
           sub,
           cpos,
           cdur) %>%
  nest() %>%
  mutate(data = purrr::map(data, qfun)) %>%
  unnest(cols = c(data)) %>% 
  mutate(quantile = as_factor(quantile)) %>% 
  ungroup() %>% 
  group_by(sub,
           cpos,
           cdur,
           quantile) %>% 
  summarise(rt_median = median(rt),
            quantile_low = HDInterval::hdi(rt, credMass = 0.89)["lower"],
            quantile_high = HDInterval::hdi(rt, credMass = 0.89)["upper"]) %>% 
  rename(`Position of c` = cpos)

# plotting
pred_rt_quantiles %>%
  filter(sub <= 4) %>%
  ggplot(aes(x = as_factor(cdur), y = rt_median,
             group = quantile,
             color = quantile)) +
  geom_ribbon(aes(ymin = quantile_low,
                  ymax = quantile_high,
                  group = quantile,
                  fill = quantile),
              alpha = 0.2,
              color = NA)+
  geom_point(data = emp_rt_quantiles[emp_rt_quantiles$sub <= 4, ],
             # geom_point(data = emp_rt_quantiles,
             mapping = aes(x = as_factor(cdur),
                           y = rt,
                           group = quantile,
                           color = quantile))+
  geom_line(data = emp_rt_quantiles[emp_rt_quantiles$sub <= 4, ],
            # geom_line(data = emp_rt_quantiles,
            mapping = aes(x = as_factor(cdur),
                          y = rt,
                          group = quantile,
                          color = quantile))+
  facet_wrap(sub ~ `Position of c`,
             labeller = label_both) +
  viridis::scale_color_viridis(direction = -1,discrete = TRUE, option = "E") +
  viridis::scale_fill_viridis(direction = -1,discrete = TRUE, option = "E", guide=F) +
  ggtitle("") +
  ggthemes::theme_tufte()+
  ylab("Response time (s)")+
  xlab("\nDuration of c (ms)")+
  labs(color = "Quantile")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size = 18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  scale_y_continuous(breaks = seq(0, 2 , by = 0.2))+
  # scale_x_continuous(guide = guide_axis(angle = 90))
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))



# individual emp and pred rt dist
#---------------------------------------------------------------------------#
emp_rt_quantiles <- emp_data %>% 
  droplevels() %>%
  group_by(sub) %>%
  summarise(median_rt = mean(rt_new))

pred_data %>% 
  ggplot()+
  geom_density(aes(x = rt,
                   fill = "Empirical",
                   color = "Empirical"),
               data = emp_data,
               alpha = 0.4)+
  geom_density(aes(x = rt,
                   color = "Predicted",
                   fill = "Predicted"),
               alpha = 0.6)+
  scale_fill_manual(name="", values=c("#2E5868", "#B06988")) +
  scale_color_manual(name="", values=c("#2E5868", "#B06988"))+ 
  geom_vline(aes(xintercept = median(rt)),
             color = "#2E5868",
             alpha = 0.7)+
  geom_vline(aes(xintercept = median_rt),
             data = emp_rt_quantiles,
             color = "#B06988",
             alpha = 0.7)+
  facet_wrap(~ sub)+
  ggthemes::theme_tufte()+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size = 18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))



# average rt quantile plot
#---------------------------------------------------------------------------#
# average over sub
emp_rt_quantiles <- emp_data %>% 
  droplevels() %>%
  group_by(cpos,
           cdur) %>%
  nest() %>%
  mutate(data = purrr::map(data, qfun_new)) %>%
  unnest(cols = c(data)) %>% 
  mutate(quantile = as_factor(quantile)) %>% 
  rename(`Position of c` = cpos)

pred_rt_quantiles_2 <- pred_data %>% 
  mutate(resp = ifelse(resp == 1, 0, 1)) %>% 
  droplevels() %>%
  group_by(sim,
           cpos,
           cdur) %>%
  nest() %>%
  mutate(data = purrr::map(data, qfun)) %>%
  unnest(cols = c(data)) %>% 
  mutate(quantile = as_factor(quantile)) %>% 
  ungroup() %>% 
  group_by(cpos,
           cdur,
           quantile) %>% 
  summarise(rt_median = median(rt),
            quantile_low = HDInterval::hdi(rt, credMass = 0.89)["lower"],
            quantile_high = HDInterval::hdi(rt, credMass = 0.89)["upper"]) %>%
  rename(`Position of c` = cpos)

# plotting
pred_rt_quantiles_2 %>%
  ggplot(aes(x = as_factor(cdur), y = rt_median,
             group = quantile,
             color = quantile)) +
  geom_ribbon(aes(ymin = quantile_low,
                  ymax = quantile_high,
                  group = quantile,
                  fill = quantile),
              alpha = 0.2,
              color = NA)+
  geom_line(linetype="dashed")+
  geom_point(data = emp_rt_quantiles,
             mapping = aes(x = as_factor(cdur),
                           y = rt,
                           group = quantile,
                           color = quantile))+
  geom_line(data = emp_rt_quantiles,
            mapping = aes(x = as_factor(cdur),
                          y = rt,
                          group = quantile,
                          color = quantile))+
  facet_grid(~`Position of c`,
             labeller = label_both) +
  viridis::scale_color_viridis(direction = -1,discrete = TRUE, option = "E") +
  viridis::scale_fill_viridis(direction = -1,discrete = TRUE, option = "E", guide=F) +
  ggtitle("") +
  ggthemes::theme_tufte()+
  ylab("Response time (s)")+
  xlab("\nDuration of c (ms)")+
  labs(color = "Quantile")+
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size = 18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0)))+
  scale_y_continuous(breaks = seq(0, 2 , by = 0.2))+
  # scale_x_continuous(guide = guide_axis(angle = 90))
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))


pred_data_df %>% 
  ggplot()+
  geom_density(aes(x=rt,
                   fill = "Prediction",
                   color = "Prediction"),
               alpha = 0.3) +
  geom_density(aes(x=rt / 1000,
                   fill = "Empirical",
                   color = "Empirical"),
               data = df,
               alpha = 0.3) +
  facet_wrap(~cpos,
             labeller = "label_both") +
  ggtitle("Predicted vs. empirical response time distribution") +
  ggthemes::theme_tufte() +
  ylab("Density") +
  xlab("\nResponse time (ms)") +
  theme(axis.line = element_line(size = .5, color = "#969696"),
        axis.ticks = element_line(color = "#969696"),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        text=element_text(size = 16),
        plot.title = element_text(size = 18,
                                  hjust = 0.5),
        axis.title.y = element_text(vjust = 0.5,
                                    margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  # scale_y_continuous(breaks = seq(0, 2 , by = 0.2))+
  # scale_x_continuous(limits = c(0, 5)) +
  # theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 1))+
  scale_fill_manual(name="", values=c("#B06988", "#2E5868")) +
  scale_color_manual(name="", values=c("#B06988", "#2E5868"))


# emp_rt_quantiles %>% 
#   ggplot(aes(x=cdur,
#              y=rt,
#              group=quantile,
#              fill=quantile))+
#   geom_line(linetype="dashed")+
#   facet_wrap(~~`Position of c`,
#              labeller = "label_both")

