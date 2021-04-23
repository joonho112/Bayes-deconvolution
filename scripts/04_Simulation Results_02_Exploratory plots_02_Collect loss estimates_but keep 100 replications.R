
###'######################################################################
###'
###' Category: Simulation Results
###' 
###' Task: Generate TABLES and FIGURES from the simulation results  
###' 
###'       (1) Collect loss estimates, but keep 100 replications
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
###' Collection of loss estimates
###' 
###' (1) SSEL for the individual effects
###' 
###' Collect the selected .csv files within the selelcted folders
###'
###'

### Collect all .csv files for the folder_list
df_collect <- collect_selected_files(folder_dir = data_dir, 
                                     folder_list = folder_list, 
                                     file_name = "02_SSEL for the individual effects.csv", 
                                     df_select)


### Tidy up the resulting dataframe
nrow(df_collect)

df_temp <- df_collect %>%
  rename(rep_ID = X) %>%
  select(-(WSSEL_PM:WSSEL_ML))
  

### Gather by estimator
vec_est <- c("ML", "PM", "CB", "GR")

df_long <- df_temp %>%
  gather(key = estimator, value = SSEL, paste0("SSEL_", vec_est)) %>%
  mutate(estimator = factor(estimator, 
                            levels = paste0("SSEL_", vec_est), 
                            labels = vec_est)) %>%
  turn_to_factors() %>%
  select(N, ICC, rsr, truth, assumed, estimator, rep_ID, SSEL) %>%
  arrange(N, ICC, rsr, truth, assumed, estimator, rep_ID) 


###' Generate a cluster variable 
###' We fitted the same model for the same data-generating mechanism
df_SSEL <- df_long %>%
  select(N, ICC, rsr, truth, rep_ID) %>%
  distinct() %>%
  mutate(cluster_ID = row_number()) %>%
  full_join_track(df_long, by = c("N", "ICC", "rsr", "truth", "rep_ID")) %>%
  select(cluster_ID, everything())


### Save as .rds and .csv format
setwd(work_dir)
saveRDS(df_SSEL, file = "datasets/All collected and replicated loss est_01_SSEL.rds")
write.csv(df_SSEL, file = "datasets/All collected and replicated loss est_01_SSEL.csv")



###'######################################################################
###'
###' Collection of loss estimates
###' 
###' (2) ISEL for EDF
###' 
###' Collect the selected .csv files within the selelcted folders
###'
###'

### Collect all .csv files for the folder_list
df_collect <- collect_selected_files(folder_dir = data_dir, 
                                     folder_list = folder_list, 
                                     file_name = "04_ISEL for the EDF.csv", 
                                     df_select)


### Tidy up the resulting dataframe
nrow(df_collect)

df_temp <- df_collect %>%
  rename(rep_ID = X) 


### Gather by estimator
vec_est <- c("ML", "PM", "CB", "GR")

df_long <- df_temp %>%
  gather(key = estimator, value = ISEL, paste0("ISEL_", vec_est)) %>%
  mutate(estimator = factor(estimator, 
                            levels = paste0("ISEL_", vec_est), 
                            labels = vec_est)) %>%
  turn_to_factors() %>%
  select(N, ICC, rsr, truth, assumed, estimator, rep_ID, ISEL) %>%
  arrange(N, ICC, rsr, truth, assumed, estimator, rep_ID) 


###' Generate a cluster variable 
###' We fitted the same model for the same data-generating mechanismdf
df_ISEL <- df_long %>%
  select(N, ICC, rsr, truth, rep_ID) %>%
  distinct() %>%
  mutate(cluster_ID = row_number()) %>%
  full_join_track(df_long, by = c("N", "ICC", "rsr", "truth", "rep_ID")) %>%
  select(cluster_ID, everything())


### Save as .rds and .csv format
setwd(work_dir)
saveRDS(df_ISEL, file = "datasets/All collected and replicated loss est_02_ISEL.rds")
write.csv(df_ISEL, file = "datasets/All collected and replicated loss est_02_ISEL.csv")



###'######################################################################
###'
###' Collection of loss estimates
###' 
###' (3) SSEL for rank estimates
###' 
###' Collect the selected .csv files within the selelcted folders
###'
###'

### Collect all .csv files for the folder_list
df_collect <- collect_selected_files(folder_dir = data_dir, 
                                     folder_list = folder_list, 
                                     file_name = "05_SEL for the ranks.csv", 
                                     df_select)


### Tidy up the resulting dataframe
nrow(df_collect)

df_temp <- df_collect %>%
  rename(rep_ID = X) %>%
  select(-SEL_rank)
  

### Gather by estimator
vec_est <- c("ML", "PM", "CB", "GR")

df_long <- df_temp %>%
  gather(key = estimator, value = SELrank, paste0("SEL_rank_", vec_est)) %>%
  mutate(estimator = factor(estimator, 
                            levels = paste0("SEL_rank_", vec_est), 
                            labels = vec_est)) %>%
  turn_to_factors() %>%
  select(N, ICC, rsr, truth, assumed, estimator, rep_ID, SELrank) %>%
  arrange(N, ICC, rsr, truth, assumed, estimator, rep_ID) 


###' Generate a cluster variable 
###' We fitted the same model for the same data-generating mechanismdf
df_SELrank <- df_long %>%
  select(N, ICC, rsr, truth, rep_ID) %>%
  distinct() %>%
  mutate(cluster_ID = row_number()) %>%
  full_join_track(df_long, by = c("N", "ICC", "rsr", "truth", "rep_ID")) %>%
  select(cluster_ID, everything())


### Save as .rds and .csv format
setwd(work_dir)
saveRDS(df_SELrank, file = "datasets/All collected and replicated loss est_03_SEL_rank.rds")
write.csv(df_SELrank, file = "datasets/All collected and replicated loss est_03_SEL_rank.csv")



###'######################################################################
###'
###' Merge the 3 loss estimates
###'
###'

list_collect <- list(df_SSEL, df_ISEL, df_SELrank)

df_combine <- reduce(list_collect, 
                     full_join_track, 
                     by = names(df_SSEL)[1:8])


### Save as .rds and .csv format
setwd(work_dir)
saveRDS(df_combine, file = "datasets/All collected and replicated loss estimates.rds")
write.csv(df_combine, file = "datasets/All collected and replicated loss estimates.csv")


