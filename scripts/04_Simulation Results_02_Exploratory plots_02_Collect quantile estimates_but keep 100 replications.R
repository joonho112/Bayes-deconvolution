
###'######################################################################
###'
###' Category: Simulation Results
###' 
###' Task: Generate TABLES and FIGURES from the simulation results  
###' 
###'       (2) Collect quantile estimates, but keep 100 replications
###'       
###' Data: Results from the simulation
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
###' Collection of quantile estimates
###' 
###' - But keep 100 replications
###' 
###' Collect the selected .csv files within the selelcted folders
###'
###'

### Collect all .csv files for the folder_list
file_name <- c("10_Percentile estimates of G_with 100 reps.csv")

df_collect <- collect_selected_files(folder_dir = data_dir, 
                                     folder_list = folder_list, 
                                     file_name = file_name, 
                                     df_select)

nrow(df_collect)

df_temp <- df_collect %>%
  dplyr::select(-X)



###'######################################################################
###'
###' Convert variables into factors
###'
###'

df_factor <- df_temp %>%
  turn_to_factors() %>%
  dplyr::select(N, ICC, rsr, truth, assumed, estimator, 
                iter, tail_p, quantile, diff) %>%
  arrange(N, ICC, rsr, truth, assumed, estimator, tail_p, iter) 



###'######################################################################
###'
###' Generate a cluster variable 
###' We fitted the same model for the same data-generating mechanism
###'
###'

df_final <- df_factor %>%
  select(N, ICC, rsr, truth, iter) %>%
  distinct() %>%
  mutate(cluster_ID = row_number()) %>%
  full_join_track(df_factor, by = c("N", "ICC", "rsr", "truth", "iter")) %>%
  select(cluster_ID, everything()) %>%
  arrange(N, ICC, rsr, truth, assumed, estimator, tail_p, iter) 



###'######################################################################
###' 
###' Save as .rds and .csv format
###' 
###' 

setwd(work_dir)
saveRDS(df_final, file = "datasets/All collected and replicated quantile est.rds")
write.csv(df_final, file = "datasets/All collected and replicated quantile est.csv")


