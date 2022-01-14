library(tidyverse)
library(magrittr)
`%notin%` = Negate(`%in%`)

exclude_exp2_cond_A <- c(6)

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

df_manipulation <- function(x, exclude) {
  x %<>%
    mutate(sub = as.numeric(sub)) %>% 
    rename(rt = RT) %>% 
    filter(rt != 0) %>%
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

df_exp2_cond_A <- df_manipulation(df_exp2_cond_A, exclude_exp2_cond_A)


# Define constants
lambda   <- 0.01  # weight param (how many previous data points)
c_0      <- 0.5   # in-control mean (expected average performance for fast guess)
sigma_0  <- 0.5   # in-control sd
L        <- 1.5   # sensitivity (1.5 is a relative low value)

# Compute EMWA
data <- data.frame()
for (s in unique(df_exp2_cond_A$sub)) {
  
  tmp <- df_exp2_cond_A %>% 
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
  scale_x_continuous(limits = c(0,4000))+
  facet_wrap(~sub)

data %<>%
  filter(rt >= 100,
         rt <= 2000)

sumsum <- data %>% 
  group_by(sub) %>% 
  summarise(n=length(rt))

write_csv(data, "/users/lukas/documents/UniHeidel/Project_Discrimination/data_exp1_cond_A_revision.csv")
