library(tidyverse)
library(magrittr)

##------------------------------------------------
## READ RAW DATA
##------------------------------------------------
setwd("/users/lukas/documents/UniHeidel/Project_Discrimination/Dyjas Bausenhart Ulrich 2012/Experiment 1")

# Read all files into one data frame
files = list.files(pattern="*D.txt")
df = NULL
for (f in files) {
  sub=as.numeric(substr(f,10,nchar(f)-5))
  if(sub==99) next
  data = read.table(f, skip=30)
  colnames(data) = c("trl", "cpos", "cdur", "resp", "RT", "type")
  data$d1 = ifelse(data$cpos==2, 500, data$cdur)
  data$d2 = ifelse(data$cpos==2, data$cdur, 500)
  
  df=rbind(
    df,
    cbind(
      sub=sub,
      data)
  )
}

# reassign sub nr
df %<>%
  mutate(sub=as.numeric(sub)) %>% 
  arrange(sub)
nr <- 1
for (i in 2:length(df$sub)){
  if(df$sub[i]!=df$sub[i-1]){
    nr <- nr+1
    df$sub[df$sub==df$sub[i]] <- nr
  }
}

# Remove NA response and fast responses
df[df$resp<1 | df$resp>2,"resp"] = NA
df = na.omit(df)

##------------------------------------------------
## REMOVE FAST GUESSES AND NON COOPERATIVE PARTICIPANTS
##------------------------------------------------
# calculate correct response
df$correct <- NA
for (i in 1:nrow(df)) {
  if(df$cpos[i]==1 & df$cdur[i]>500){
    if(df$resp[i]==1){
      df$correct[i] <- 1
    }else{
      df$correct[i] <- 0
    }
  }else if(df$cpos[i]==1 & df$cdur[i]<500){
    if(df$resp[i]==2){
      df$correct[i] <- 1
    }else{
      df$correct[i] <- 0
    }
  }else if(df$cpos[i]==2 & df$cdur[i]>500){
    if(df$resp[i]==2){
      df$correct[i] <- 1
    }else{
      df$correct[i] <- 0
    }
  }else if(df$cdur[i]==500){
    df$correct[i] <- 0
  }
  else{
    if(df$resp[i]==1){
      df$correct[i] <- 1
    }else{
      df$correct[i] <- 0
    }
  }
}

# Define constants
lambda   <- 0.01  # weight param (how many previous data points)
c_0      <- 0.5   # in-control mean (expected average performance for fast guess)
sigma_0  <- 0.5   # in-control sd
L        <- 1.5   # sensitivity (1.5 is a relative low value)

# Compute EMWA
data <- data.frame()
for (s in 1:length(unique(df$sub))) {
  
  tmp <- df %>% 
    filter(sub==s) %>% 
    arrange(RT)
  
  tmp$c_s <- NA
  tmp$inCm <- "empty"
  tmp$UCL <- NA
  
  c_before <- c_0
  
  # Compute EWMA for each data point
  for (i in 1:length(tmp$RT)) {
    if(tmp$cdur[i]==500) next
    tmp$c_s[i] <- (lambda*tmp$correct[i]) + ((1-lambda)*c_before)
    tmp$UCL[i] <- c_0 + (L*sigma_0)*sqrt((lambda/(2-lambda))*(1-(1-lambda)^(2*i)))
    tmp$inCm[i] <- tmp$c_s[i] < tmp$UCL[i]
    c_before <- tmp$c_s[i]
  }
  
  data <- rbind(data,tmp)
}

# remove non-cooperative participants and fast guesses
df <- data %>% 
  filter(sub!=6,
         sub!=10,
         sub!=18,
         sub!=20,
         sub!=24) %>% 
  mutate(inCm=as.character(inCm)) %>% 
  filter(inCm!="TRUE")

# remove RT's below 100 ms
df %<>%
  filter(RT>100)

# reassign sub nr
df %<>%
  mutate(sub=as.numeric(sub)) %>% 
  arrange(sub)
nr <- 1
for (i in 2:length(df$sub)){
  if(df$sub[i]!=df$sub[i-1]){
    nr <- nr+1
    df$sub[df$sub==df$sub[i]] <- nr
  }
}

# arrange rows of df after sub and then trial
df <- df %>% 
  arrange(sub,trl)

##------------------------------------------------
## SAVE DATASET
##------------------------------------------------
setwd("/users/lukas/documents/github/timing_discrimination/data")
write_csv(df,"exp1_cond_C.csv")


