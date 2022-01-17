library(tidyverse)
library(magrittr)
`%notin%` = Negate(`%in%`)

qfun <- function(x, quantiles = c(0.1, 0.3, 0.5, 0.7, 0.9)) {
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
emp_data$rt_new <- NA
for (i in 1:nrow(emp_data)){
  if (emp_data$d2[i] < emp_data$d1[i]){
    emp_data$rt_new[i] <- emp_data$rt[i] + ((emp_data$d2[i] - emp_data$d1[i]) / 1000)
  }
  else {
    emp_data$rt_new[i] <- emp_data$rt[i]
  }
}

# rt manipulation
# whenever d2 was longer than d1, then the rt increased by the difference d2 - d1
emp_data$rt_new <- NA
for (i in 1:nrow(emp_data)){
  if (emp_data$d2[i] > emp_data$d1[i]){
    emp_data$rt_new[i] <- emp_data$rt[i] + ((emp_data$d2[i] - emp_data$d1[i]) / 1000)
  }
  else {
    emp_data$rt_new[i] <- emp_data$rt[i]
  }
}

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
lower_rt_cutoff <- 0.1
fast_guesses_summary <- emp_data %>% 
  filter(cdur != 500,
         rt <= lower_rt_cutoff) %>% 
  group_by(sub) %>% 
  summarise(accuracy = round(mean(correct), 3),
            mean_rt = round(mean(rt), 3),
            n = length(correct))

# check accuracy of delayed start up's
upper_rt_cutoff <- 2.5
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
  for (i in 1:length(tmp$rt_new)) {
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

# summarize emp data
emp_rt_quantiles <- data %>% 
  filter(inCm == FALSE,
         rt_new <= 2.5) %>%
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
  # filter(sub <= 6) %>%
  ggplot(aes(x = cdur,
             y = rt,
             group = quantile,
             color = quantile)) +
  geom_point()+
  geom_line()+
  facet_wrap(~`Position of c`,
             labeller = label_both) +
  scale_y_continuous(breaks = seq(0, 2 , by = 0.2))+
  # scale_x_continuous(guide = guide_axis(angle = 90))+
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


performance_summary <- data %>% 
  group_by(sub) %>% 
  summarise(accuracy = mean(correct),
            median_rt = median(rt),
            # median_rt_new = median(rt_new),
            n_trials = length(rt))

# check accuracy of fast guesses
lower_rt_cutoff <- 0.1
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
  filter(rt_new > 0.1)


write_csv(data, "/users/lukas/documents/github/timing_discrimination/data/data_exp2_cond_D_revision_experimental.csv")
