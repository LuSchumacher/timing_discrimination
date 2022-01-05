library(tidyverse)
library(magrittr)
`%notin%` = Negate(`%in%`)

qfun <- function(x, quantiles = c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  enframe(quantile(x$rt, probs = quantiles, na.rm = TRUE),
          name = "quantile",
          value = "rt")}

exclude_exp1_cond_A <- c(10, 24)
exclude_exp1_cond_B <- c(6, 7, 10, 17, 20, 22, 23, 25)
exclude_exp1_cond_D <- c(6, 10, 18, 20, 24)
exclude_exp2_cond_A <- c(6)
exclude_exp2_cond_B <- c(6, 7, 19)
exclude_exp2_cond_D <- c(10, 19)
##------------------------------------------------------------------------##
## READ ALL DATA
##------------------------------------------------------------------------##
setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/Dyjas Bausenhart Ulrich 2012/Experiment 1")
# Exp 1, Cond A
files <- list.files(pattern = "*A.txt")
df_exp1_cond_A <- NULL
for (f in files) {
  sub <- as.numeric(substr(f, 10, nchar(f) - 5))
  if (sub == 99) next
  data <- read.table(f, skip = 30)
  colnames(data) <- c("trl", "cpos", "cdur", "resp", "RT", "type")
  data$d1 <- ifelse(data$cpos == 2, 500, data$cdur)
  data$d2 <- ifelse(data$cpos == 2, data$cdur, 500)
  df_exp1_cond_A <- rbind(df_exp1_cond_A, cbind(sub = sub, data)
  )
}

# Exp 1, Cond B
files <- list.files(pattern = "*B.txt")
df_exp1_cond_B <- NULL
for (f in files) {
  sub <- as.numeric(substr(f, 10, nchar(f) - 5))
  if (sub == 99) next
  data <- read.table(f, skip = 30)
  colnames(data) <- c("trl", "cpos", "cdur", "resp", "RT", "type")
  data$d1 <- ifelse(data$cpos == 2, 500, data$cdur)
  data$d2 <- ifelse(data$cpos == 2, data$cdur, 500)
  df_exp1_cond_B <- rbind(df_exp1_cond_B, cbind(sub = sub, data)
  )
}

# Exp 1, Cond D 
files <- list.files(pattern = "*D.txt")
df_exp1_cond_D <- NULL
for (f in files) {
  sub <- as.numeric(substr(f, 10, nchar(f) - 5))
  if (sub == 99) next
  data <- read.table(f, skip = 30)
  colnames(data) <- c("trl", "cpos", "cdur", "resp", "RT", "type")
  data$d1 <- ifelse(data$cpos == 2, 500, data$cdur)
  data$d2 <- ifelse(data$cpos == 2, data$cdur, 500)
  df_exp1_cond_D <- rbind(df_exp1_cond_D, cbind(sub = sub, data)
  )
}

setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/Dyjas Bausenhart Ulrich 2012/Experiment 2")
# Exp 2, Cond A
files <- list.files(pattern = "*A.txt")
df_exp2_cond_A <- NULL
for (f in files) {
  sub <- as.numeric(substr(f, 10, nchar(f) - 5))
  if (sub == 99) next
  data <- read.table(f, skip = 30)
  colnames(data) <- c("trl", "cpos", "cdur", "resp", "RT", "type")
  data$d1 <- ifelse(data$cpos == 2, 500, data$cdur)
  data$d2 <- ifelse(data$cpos == 2, data$cdur, 500)
  df_exp2_cond_A <- rbind(df_exp2_cond_A, cbind(sub = sub, data)
  )
}

setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/Dyjas Bausenhart Ulrich 2012/Experiment 2")
# Exp 2, Cond B
files <- list.files(pattern = "*B.txt")
df_exp2_cond_B <- NULL
for (f in files) {
  sub <- as.numeric(substr(f, 10, nchar(f) - 5))
  if (sub == 99) next
  data <- read.table(f, skip = 30)
  colnames(data) <- c("trl", "cpos", "cdur", "resp", "RT", "type")
  data$d1 <- ifelse(data$cpos == 2, 500, data$cdur)
  data$d2 <- ifelse(data$cpos == 2, data$cdur, 500)
  df_exp2_cond_B <- rbind(df_exp2_cond_B, cbind(sub = sub, data)
  )
}

setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/Dyjas Bausenhart Ulrich 2012/Experiment 2")
# Exp 2, Cond D
files <- list.files(pattern = "*D.txt")
df_exp2_cond_D <- NULL
for (f in files) {
  sub <- as.numeric(substr(f, 10, nchar(f) - 5))
  if (sub == 99) next
  data <- read.table(f, skip = 30)
  colnames(data) <- c("trl", "cpos", "cdur", "resp", "RT", "type")
  data$d1 <- ifelse(data$cpos == 2, 500, data$cdur)
  data$d2 <- ifelse(data$cpos == 2, data$cdur, 500)
  df_exp2_cond_D <- rbind(df_exp2_cond_D, cbind(sub = sub, data)
  )
}

##------------------------------------------------------------------------##
## REASSGIN SUB NR & REMOVE NA RESPONSES & EVALUATE CORRECTNESS OF RESPONSE
##------------------------------------------------------------------------##
df_manipulation <- function(x, exclude) {
  x %<>%
    mutate(sub = as.numeric(sub)) %>% 
    rename(rt = RT) %>% 
    filter(rt != 0) %>% # sub %notin% exclude
    arrange(sub)
  
  nr <- 1
  for (i in 2:length(x$sub)){
    if (x$sub[i] != x$sub[i - 1]){
      nr <- nr + 1
      x$sub[x$sub == x$sub[i]] <- nr
    }
  }
  
  x %<>%
    filter(sub %notin% exclude) %>%
    arrange(sub)
  
  nr <- 1
  for (i in 2:length(x$sub)){
    if (x$sub[i] != x$sub[i - 1]){
      nr <- nr + 1
      x$sub[x$sub == x$sub[i]] <- nr
    }
  }
  
  x[x$resp < 1 | x$resp > 2, "resp"] <- NA
  x <- na.omit(x)
  
  # calculate correct response
  x$correct <- NA
  for (i in 1:nrow(x)) {
    if(x$cpos[i] == 1 & x$cdur[i] > 500){
      if(x$resp[i] == 1){
        x$correct[i] <- 1
      }else{
        x$correct[i] <- 0
      }
    }else if(x$cpos[i] == 1 & x$cdur[i] < 500){
      if(x$resp[i] == 2){
        x$correct[i] <- 1
      }else{
        x$correct[i] <- 0
      }
    }else if(x$cpos[i] == 2 & x$cdur[i] > 500){
      if(x$resp[i] == 2){
        x$correct[i] <- 1
      }else{
        x$correct[i] <- 0
      }
    }else if(x$cdur[i] == 500){
      x$correct[i] <- 0
    }
    else{
      if(x$resp[i] == 1){
        x$correct[i] <- 1
      }else{
        x$correct[i] <- 0
      }
    }
  }
  return(x)
}

df_exp1_cond_A <- df_manipulation(df_exp1_cond_A, exclude_exp1_cond_A)
df_exp1_cond_B <- df_manipulation(df_exp1_cond_B, exclude_exp1_cond_B)
df_exp1_cond_D <- df_manipulation(df_exp1_cond_D, exclude_exp1_cond_D)
df_exp2_cond_A <- df_manipulation(df_exp2_cond_A, exclude_exp2_cond_A)
df_exp2_cond_B <- df_manipulation(df_exp2_cond_B, exclude_exp2_cond_B)
df_exp2_cond_D <- df_manipulation(df_exp2_cond_D, exclude_exp2_cond_D)

##------------------------------------------------------------------------##
## SUMMARIZE PERFORMANCE
##------------------------------------------------------------------------##
get_performance_summary <- function(x) {
  x %<>% 
    group_by(sub) %>% 
    summarise(accuracy = mean(correct),
              median_rt = median(rt),
              n_trials = length(rt))
  return(x)
}

summary_exp1_cond_A <- get_performance_summary(df_exp1_cond_A)
summary_exp1_cond_B <- get_performance_summary(df_exp1_cond_B)
summary_exp1_cond_D <- get_performance_summary(df_exp1_cond_D)
summary_exp2_cond_A <- get_performance_summary(df_exp2_cond_A)
summary_exp2_cond_B <- get_performance_summary(df_exp2_cond_B)
summary_exp2_cond_D <- get_performance_summary(df_exp2_cond_D)


get_proportion_fast_guesses <- function(x, rt_lim) {
  x %<>%
    filter(cdur != 500,
           rt < rt_lim) %>% 
    group_by(sub) %>% 
    summarise(accuracy = round(mean(correct), 3),
              mean_rt = round(mean(rt), 3),
              n = length(correct))
  return(x)
}

rt_lim <- 200
fast_guesses_exp1_cond_A <- get_proportion_fast_guesses(df_exp1_cond_A, rt_lim)
fast_guesses_exp1_cond_B <- get_proportion_fast_guesses(df_exp1_cond_B, rt_lim)
fast_guesses_exp1_cond_D <- get_proportion_fast_guesses(df_exp1_cond_D, rt_lim)
fast_guesses_exp2_cond_A <- get_proportion_fast_guesses(df_exp2_cond_A, rt_lim)
fast_guesses_exp2_cond_B <- get_proportion_fast_guesses(df_exp2_cond_B, rt_lim)
fast_guesses_exp2_cond_D <- get_proportion_fast_guesses(df_exp2_cond_D, rt_lim)

##------------------------------------------------------------------------##
## PLOT RT QUANTILES
##------------------------------------------------------------------------##
plot_rt_quantiles <- function(x) {
  emp_rt_quantiles <- x %>% 
    droplevels() %>%
    group_by(cpos,
             cdur) %>%
    nest() %>%
    mutate(data = purrr::map(data, qfun)) %>%
    unnest(cols = c(data)) %>% 
    mutate(quantile = as_factor(quantile)) %>% 
    mutate(rt = rt / 1000) %>% 
    rename(`Position of c` = cpos)
  
  emp_rt_quantiles %>%
    ggplot(aes(x = as_factor(cdur),
               y = rt,
               group = quantile,
               color = quantile)) +
    geom_point()+
    geom_line()+
    facet_grid(~`Position of c`,
               labeller = label_both) +
    viridis::scale_color_viridis(direction = -1,discrete = TRUE, option = "E") +
    viridis::scale_fill_viridis(direction = -1,discrete = TRUE, option = "E", guide="none") +
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
}

plot_rt_quantiles(df_exp1_cond_A)
plot_rt_quantiles(df_exp1_cond_B)
plot_rt_quantiles(df_exp1_cond_D)
plot_rt_quantiles(df_exp2_cond_A)
plot_rt_quantiles(df_exp2_cond_B)
plot_rt_quantiles(df_exp2_cond_D)

check_fast_resp <- df_exp1_cond_B %>%
  filter(rt < 100) %>%
  mutate(acc = mean(correct))


##------------------------------------------------------------------------##
## PREPARE DATA FROM EXP 2 COND D
##------------------------------------------------------------------------##
# Define constants
lambda   <- 0.01  # weight param (how many previous data points)
c_0      <- 0.5   # in-control mean (expected average performance for fast guess)
sigma_0  <- 0.5   # in-control sd
L        <- 1.5   # sensitivity (1.5 is a relative low value)

# Compute EMWA
data <- data.frame()
for (s in unique(df_exp2_cond_D$sub)) {
  
  tmp <- df_exp2_cond_D %>% 
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
  scale_x_continuous(limits = c(0,2000))+
  facet_wrap(~sub)

data %<>%
  filter(rt >= 100,
         rt <= 2000)

sumsum <- data %>% 
  group_by(sub) %>% 
  summarise(n=length(rt))

write_csv(data, "/users/lukas/documents/UniHeidel/Project_Discrimination/data_exp2_condD_revision.csv")
