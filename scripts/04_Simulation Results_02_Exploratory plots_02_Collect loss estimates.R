
###'######################################################################
###'
###' Category: Simulation Results
###' 
###' Task: Generate TABLES and FIGURES from the simulation results  
###' 
###'       (2) Collect loss estimates
###'       
###' Data: Results from the simulation
###' 
###' Data: 2020-04-09
###' 
###' Author: JoonHo Lee (`joonho@berkeley.edu`)
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
data_dir <- c("D:/Data/DP_Simulation")


### Call libraries
library(tidyverse)
library(cowplot)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "functions/03_Simulation Implementation_04_Plot helpers.R")
source(file = "functions/03_Simulation Implementation_05_Data management helpers.R")



###'######################################################################
###'
###' Load the table containing the simulation parameters
###'
###'

setwd(work_dir)
df_sim <- read.csv(file = "tables/simulation parameters_2_All with dg_name mark.csv") %>%
  select(-X)



###'######################################################################
###'
###' Define a function to define factors
###'
###'

turn_to_factors <- function(df_collect){
  
  ### Define factor levels and labels
  lev_truth <- lab_truth <- c("Gaussian", "T", "ALD", "Bimodal", "Skew", "Mixed")
  lev_assumed <- c("Gaussian", "T", "DP-diffuse", "DP-EB")
  lab_assumed <- c("Gaussian", "T", "DP-diffuse", "DP-Inform")
  lev_est <- lab_est <- c("ML", "PM", "CB", "GR")
  
  lev_N <- c(25, 50, 100, 200)
  lev_ICC <- c(0.1, 0.5, 0.9)
  lev_rsr <- c(1, 5, 10)
  
  lab_N <- paste0("N = ", lev_N)
  lab_ICC <- paste0("I = ", lev_ICC)
  lab_rsr <- paste0("R = ", lev_rsr)
  
  
  ### Assign factor levels and labels
  df_collect %>%
    mutate(truth = factor(truth, levels = lev_truth, labels = lab_truth), 
           assumed = factor(assumed, levels = lev_assumed, labels = lab_assumed),
           N = factor(N, levels = lev_N, labels = lab_N), 
           ICC = factor(ICC, levels = lev_ICC, labels = lab_ICC), 
           rsr = factor(rsr, levels = lev_rsr, labels = lab_rsr), 
           estimator = factor(estimator, levels = lev_est, labels = lab_est))
}



###'######################################################################
###'
###' Conditionally select a subgroup of folders 
###'
###'

### Parameter selection => drop DP_EB for now
tabdf(df_sim, assumed)

# df_select <- df_sim %>%
#   filter(assumed != c("DP-EB"))

df_select <- df_sim

nrow(df_select)


### Return a vector for the selected subgroup of folders
folder_list <- folder_select(df_select, data_dir, leading_zeros = FALSE)



###'######################################################################
###'
###' (1) Collection of loss estimates
###' 
###' Collect the selected .csv files within the selelcted folders
###'
###'

df_collect <- collect_selected_files(folder_dir = data_dir, 
                                     folder_list = folder_list, 
                                     file_name = "09_Collection of loss estimates.csv", 
                                     df_select)

### Extract true variances
df_truevar <- df_collect %>%
  filter(estimator == "true")

df_collect <- df_collect %>% 
  filter(estimator != "true")


### Turn to factors & tidy up
df_loss <- turn_to_factors(df_collect) %>%
  select(-X) %>%
  select(quantity, truth, N, ICC, rsr, assumed, estimator, value) %>%
  arrange(quantity, truth, N, ICC, rsr, assumed, estimator)


### Save as .rds and .csv format
setwd(work_dir)
saveRDS(df_loss, file = "datasets/All collection of loss estimates.rds")
write.csv(df_loss, file = "datasets/All collection of loss estimates.csv")



###'######################################################################
###'
###' (2) Percentile estimates of G
###' 
###' Collect the selected .csv files within the selelcted folders
###'
###'

df_collect <- collect_selected_files(folder_dir = data_dir, 
                                     folder_list = folder_list, 
                                     file_name = "06_Percentile estimates of G.csv", 
                                     df_select)


### Turn to factors & tidy up
vec_tail_p <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95) 

df_pct <- turn_to_factors(df_collect) %>% 
  mutate(tail_p = factor(tail_p, levels = vec_tail_p)) %>%
  select(-X) %>%
  select(truth, N, ICC, rsr, assumed, estimator, 
         tail_p, quantile, diff) %>%
  arrange(truth, N, ICC, rsr, assumed, estimator, tail_p)


### Save as .rds and .csv format
setwd(work_dir)
saveRDS(df_pct, file = "datasets/All percentile estimates of G.rds")
write.csv(df_pct, file = "datasets/All percentile estimates of G.csv")
  


###'######################################################################
###'
###' (3) Histogram bars 
###'  
###' Collect the selected .csv files within the selelcted folders
###'
###'

df_collect <- collect_selected_files(folder_dir = data_dir, 
                                     folder_list = folder_list, 
                                     file_name = "07_Histogram bars.csv", 
                                     df_select)


### Rename estimator variables and convert to the long format
df_collect <- df_collect %>%
  rename(PM = theta_pm, CB = theta_cb, GR = theta_gr, ML = theta_ml) %>%
  gather(key = estimator, value = value, PM, CB, GR, ML) %>%
  rename(Bin = X)


### Turn to factors & tidy up
df_hist <- turn_to_factors(df_collect) %>% 
  select(truth, N, ICC, rsr, assumed, estimator, 
         Bin, value) %>%
  arrange(truth, N, ICC, rsr, assumed, estimator)


### Save as .rds and .csv format
setwd(work_dir)
saveRDS(df_hist, file = "datasets/All histogram bars.rds")
write.csv(df_hist, file = "datasets/All histogram bars.csv")

