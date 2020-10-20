
###'######################################################################
###'
###' Category: Publication Tables and Figures 
###' 
###' Task: Regression tables
###' 
###'  [Table 03]. Meta-model regression results for percentile bias
###'           
###' Data: Simulated data
###' 
###' Data: 2020-04-29
###' 
###' Author: JoonHo Lee (joonho@berkeley.edu)
###' 
###' 

###'######################################################################
###'
###' Basic settings
###'
###'

### Start with a clean slate
gc(); rm(list=ls())   


### Set working directory 
work_dir <- c("~/Bayes-deconvolution")
setwd(work_dir)


### Set a data directory
data_dir <- file.path(work_dir, "datasets", 
                      "12_All collected and replicated quantile estimates")


### Call libraries
library(tidyverse)
library(cowplot)
library(multcomp)
library(jtools)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "functions/03_Simulation Implementation_04_Plot helpers.R")
source(file = "functions/03_Simulation Implementation_05_Data management helpers.R")



###'######################################################################
###'
###' Load the processed collection of  quantile estimates
###'
###'

### Load the dataset
setwd(data_dir)

df <- readRDS(file = "All collected and replicated quantile est.rds") 


### Convert the cluster_ID into factor
df <- df %>%
  mutate(cluster_ID = factor(cluster_ID)) %>%
  mutate(tail_p = tail_p*100, 
         quantile = quantile*100, 
         diff = diff*100)

classmode(df, everything())



###'######################################################################
###'
###' Simplify the factor labels
###'
###'

### Rename variables
df <- df %>%
  rename(I = ICC, R = rsr, anl = assumed, est = estimator)


### Define new labels
lev_N <- unique(df$N) 
lab_N <- lev_N %>% str_remove(" = ")

lev_I <- unique(df$I)
lab_I <- lev_I %>% str_remove(" = ") 

lev_R <- unique(df$R)
lab_R <- lev_R %>% str_remove(" = ") 

lev_anl <- unique(df$anl)
lab_anl <- lev_anl %>% str_replace("-", "_")


### Mutate factors
df <- df %>%
  mutate(N = factor(N, levels = lev_N, labels = lab_N), 
         I = factor(I, levels = lev_I, labels = lab_I),
         R = factor(R, levels = lev_R, labels = lab_R),
         anl = factor(anl, levels = lev_anl, labels = lab_anl))



###'######################################################################
###'
###' Subset only conditional cases
###'
###'

### Tabulate existing variables
names(df)
tabdf(df, truth)
tabdf(df, N)
tabdf(df, anl)
tabdf(df, I)
tabdf(df, est)


### Filter only conditional cases
df_cond <- df %>%
  filter(truth %in% c("Gaussian", "ALD", "Mixed")) %>%
  filter(est != c("ML")) %>%
  filter(anl != c("T")) %>%
  filter(!(N %in% c("N100"))) 



###'######################################################################
###'
###' Fit meta-models & save the tables as .csv file
###' 
###' 

### Prepare empty lists to collect results
list_fit <- list()
list_df_fit <- list()
list_df_fit_pub <- list()


### Generate a table for looping
vec_truth <- unique(df_cond$truth)
vec_tail_p <- unique(df_cond$tail_p)
df_grid <- expand.grid(vec_tail_p, vec_truth) %>%
  data.frame() %>%
  dplyr::select(Var2, Var1) %>%
  set_names(c("truth", "tail_p")) 


for (i in seq(nrow(df_grid))){
  
  ### Extract elements
  elem_truth <- df_grid$truth[i]
  elem_tail_p <- df_grid$tail_p[i]
  
    
  ### Filter only conditional cases
  df_sub <- df_cond %>%
    filter(truth %in% elem_truth) %>%
    filter(tail_p %in% elem_tail_p)
  
  names(df_sub)
  
  
  ### Set reference category for each factor
  ref_N <- "N25"
  ref_I <- "I0.1"
  ref_R <- "R1"
  ref_anl <- "Gaussian"
  ref_est <- "PM"
  
  vec_ref <- paste(ref_N, ref_I, ref_R, ref_anl, ref_est, sep = ", ")
  
  df_sub_relevel <- df_sub %>%
    mutate(N = relevel(N, ref = ref_N), 
           I = relevel(I, ref = ref_I), 
           R = relevel(R, ref = ref_R), 
           anl = relevel(anl, ref = ref_anl), 
           est = relevel(est, ref = ref_est))
  
  
  ### Convert diff scores into absolute or raw values
  df_sub_relevel <- df_sub_relevel 
  
  
  ### Fit the OLS model with clustered SE
  fit <- lm(diff ~ N + I + R + anl + est + 
              N:I + N:R + N:anl + N:est + 
              I:R + I:anl + I:est + 
              R:anl + R:est + 
              anl:est, 
            data = df_sub_relevel)
  
  df_fit <- get_robust_se(fit, cluster = "cluster_ID")$coefs %>% 
    data.frame()
  
  df_fit_pub <- df_fit %>%
    dplyr::select(Est., t.value) %>%
    rownames_to_column() %>%
    mutate(Est. = round(Est., 1), 
           t.value = round(t.value, 1)) %>%
    set_names(paste0(c("predictor", "estimate", "t-value"), "_", elem_truth))
  
  
  ### Collect the resulting dataframes
  list_fit[[i]] <- fit
  list_df_fit[[i]] <- df_fit
  list_df_fit_pub[[i]] <- df_fit_pub
  
}


### Bind columns and save the results
setwd("~/Bayes-deconvolution/tables/01_Regression tables")

saveRDS(list_fit, file = "01_Regression object_04_Raw bias of quantile estimates.rds")

write.csv(bind_cols(list_df_fit), 
          file = "01_Regression table_04_Raw bias of quantile estimates.csv")

write.csv(bind_cols(list_df_fit_pub), 
          file = "01_Regression table_04_Raw bias of quantile estimates_Pub.csv")





