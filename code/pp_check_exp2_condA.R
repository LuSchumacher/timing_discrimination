library(tidyverse)
library(rstan)
library(magrittr)
setwd("/Users/lukas/documents/UniHeidel/code")
source("wienerProcess2.R")

df <- read_csv("/users/lukas/documents/UniHeidel/Project_Discrimination/data_exp1_cond_A_revision.csv")
fit <- readRDS("/users/lukas/documents/UniHeidel/Project_Discrimination/fit_model_6_cond_A.rds")

# extract posterior distributions
mat <- as.data.frame(fit) %>%
  select(starts_with(c("a[", "ndt[", "z0[", "v0[", "b2v[")))

numextract <- function(string){ 
  str_extract(string, "\\-*\\d+\\.*\\d*")
}

a <- mat %>% 
  select(`a[1]`:`a[24]`) %>% 
  pivot_longer(cols = `a[1]`:`a[24]`,
               names_to = "sub",
               values_to = "a") %>% 
  mutate(sample=rep(seq(1,nrow(mat)), each=nrow(.) / nrow(mat)),
         sub=numextract(sub)) %>% 
  select(sample,sub,a)

ndt <- mat %>% 
  select(`ndt[1]`:`ndt[24]`) %>% 
  pivot_longer(cols = `ndt[1]`:`ndt[24]`,
               names_to = "sub",
               values_to = "ndt") %>% 
  select(ndt)

z0 <- mat %>% 
  select(`z0[1]`:`z0[24]`) %>% 
  pivot_longer(cols = `z0[1]`:`z0[24]`,
               names_to = "sub",
               values_to = "z0") %>% 
  select(z0)

v0 <- mat %>% 
  select(`v0[1]`:`v0[24]`) %>% 
  pivot_longer(cols = `v0[1]`:`v0[24]`,
               names_to = "sub",
               values_to = "v0") %>% 
  select(v0)

b2v <- mat %>% 
  select(`b2v[1]`:`b2v[24]`) %>% 
  pivot_longer(cols = `b2v[1]`:`b2v[24]`,
               names_to = "sub",
               values_to = "b2v") %>% 
  select(b2v)

mat <- cbind(a, ndt, z0, v0, b2v)
mat %<>%
  mutate(sub=as.numeric(sub))


#-------------------------------------------------------------------------------#
# SIMULATION
#-------------------------------------------------------------------------------#
# sample posteriors
n_samples <- 500
posterior_samples <- sample(max(mat$sample), n_samples)
posteriors <- mat %>% 
  filter(sample %in% posterior_samples)

n_sub <- max(mat$sub)

col_names <- c("sim","sub","d1","d2","cdur","cpos","resp", "rt")
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
    z   <- mat$z[mat$sample == sim][s]
    v0  <- mat$v0[mat$sample == sim][s]
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
      delta <- v0 + b2v * d2[i]
      diff_process <- wienerProcess(v=delta, a=a, z=z, ndt=ndt)
      pred_data[(df_nrow*x + cum_nrow + i), ] <- c(x, s, d1[i], d2[i], cdur[i], cpos[i],
                                                   diff_process[1], diff_process[2])
      
    }
    cum_nrow <- cum_nrow + n_trials
  }
  end_time <- Sys.time()
  print(end_time - start_time)
  x <- x + 1
}

# write_csv(pred_data_df, "/Users/lukas/Documents/UniHeidel/Project_Discrimination/code/pp_check_exp2_cond_A.csv")

#-------------------------------------------------------------------------------#
# SUMMARY AND PLOTTING
#-------------------------------------------------------------------------------#
pred_data <- read_csv("/Users/lukas/Documents/UniHeidel/Project_Discrimination/code/pp_check_exp2_cond_A.csv")
pred_data <- as.data.frame(pred_data)

qfun <- function(x, quantiles = c(0.1, 0.3, 0.5, 0.7, 0.9)) {
  enframe(quantile(x$rt, probs = quantiles, na.rm = TRUE),
          name = "quantile",
          value = "rt")}


# individual rt quantile plot
#-------------------------------------------------------------------------------#
# summarize emp data
emp_rt_quantiles <- df %>% 
  rename(rt = rt) %>% 
  droplevels() %>%
  group_by(sub,
           cpos,
           cdur) %>%
  nest() %>%
  mutate(data = purrr::map(data, qfun)) %>%
  unnest(cols = c(data)) %>% 
  mutate(quantile = as_factor(quantile)) %>% 
  mutate(rt = rt / 1000) %>% 
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
#-------------------------------------------------------------------------------#
emp_rt_quantiles <- df %>% 
  rename(rt = rt) %>% 
  droplevels() %>%
  group_by(sub) %>%
  summarise(median_rt = mean(rt)) %>% 
  mutate(median_rt = median_rt / 1000)

pred_data %>% 
  ggplot()+
  geom_density(aes(x = rt / 1000,
                   fill = "Empirical",
                   color = "Empirical"),
               data = df,
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
#-------------------------------------------------------------------------------#
# average over sub
emp_rt_quantiles <- df %>% 
  droplevels() %>%
  group_by(cpos,
           cdur) %>%
  nest() %>%
  mutate(data = purrr::map(data, qfun)) %>%
  unnest(cols = c(data)) %>% 
  mutate(quantile = as_factor(quantile)) %>% 
  mutate(rt = rt / 1000) %>% 
  rename(`Position of c` = cpos)

pred_rt_quantiles_2 <- pred_data_df %>% 
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
  # summarise(rt_median = median(rt),
  #           quantile_low = quantile(rt, 0.025),
  #           quantile_high = quantile(rt, 0.975)) %>% 
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

