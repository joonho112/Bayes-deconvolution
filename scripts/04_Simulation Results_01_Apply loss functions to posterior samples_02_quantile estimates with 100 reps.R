
###'######################################################################
###'
###' Category: Simulation Results
###' 
###' Task: Generate TABLES and FIGURES from the simulation results  
###' 
###'       (1) Collect quantile estimates, but keep 100 replications
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
df_sim <- read.csv(file = "tables/simulation parameters_2_All with dg_name mark.csv")



###'######################################################################
###'
###' Conditionally select a subgroup of folders 
###'
###'

### Parameter selection => drop DP_EB for now
tabdf(df_sim, assumed)

# df_select <- df_sim %>%
#   filter(assumed == c("DP-EB"))

df_select <- df_sim

nrow(df_select)


### Return a vector for the selected subgroup of folders
folder_list <- folder_select(df_select, data_dir, leading_zeros = FALSE)



###'######################################################################
###'
###' Start loops over simulation conditions
###'
###'

### The number of simulation conditions 
N_cond <- length(folder_list)   


### The number of Monte Carlo iterations to perform
N_sim <- 100


### Mark start time
start_time <- Sys.time()


for (j in 420:N_cond){
  
  ###'######################################################################
  ###'
  ###' Set required objects 
  ###' 
  ###'
  
  ### Extract simulation parameters from the folder name 
  sim_parm <- df_select %>%
    filter(fname == folder_list[j]) %>%
    as.vector()
  
  
  ### Extract simulation condition elements
  gen_name <- paste0("True_", sim_parm$truth)
  anl_name <- paste0("Model_", sim_parm$assumed)
  ICC <- sim_parm$ICC  
  rsr <- sim_parm$rsr
  nu <- sim_parm$df
  K_max <- sim_parm$N
  
  
  ### Define the number of posterior summary methods
  post_sum_methods <- c("PM", "CB", "GR", "ML")
  J_max <- length(post_sum_methods)
  
  
  
  ###'######################################################################
  ###'
  ###' Load the generated data
  ###'
  ###'
  
  temp_path <- file.path(work_dir, 
                         "datasets", 
                         "05_Pre-generated simulated data", 
                         paste0(sim_parm$dg_name, ".rds"))
  
  list_pregen <- readRDS(file = temp_path)
  
  
  
  ###'######################################################################
  ###'
  ###' Load the calculated posterior samples
  ###'
  ###'
  
  folder_dir <- file.path(data_dir, folder_list[j])
  
  setwd(folder_dir)
  
  list_posterior <- readRDS(file = "100 posterior samples.rds")
  
  
  
  ###'######################################################################
  ###'
  ###' Generate empty arrays or vectors 
  ###' to save or accumulate resulting losses
  ###' 
  ###' 
  
  ### Tail quantiles
  tail_p <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95) 
  # sim_tail <- rep(0, length(tail_p)*J_max)
  sim_tail <- NULL
  
  
  ###'######################################################################
  ###'
  ###' Calculate loss quantities across Mote Carlo iterations (N_sim)
  ###'
  ###'
  
  for (i in seq(N_sim)){
    
    
    ###'######################################################################
    ###'
    ###' Extract generated data (per iteration)
    ###'
    ###'
    
    list_meta <- list_pregen[[i]]$list_data
    
    df_meta <- list_pregen[[i]]$df_data 
    
    priquan <- list_pregen[[i]]$CDF_theory
    
    
    
    ###'######################################################################
    ###'
    ###' Extract posterior samples (per iteration)
    ###'
    ###'
    
    df_posterior <- list_posterior[i][[1]]
    
    
    
    ###'######################################################################
    ###'
    ###' (3) Compute the site-specific estimates
    ###' 
    ###' - Posterior means (PM) and standard deviations (PSD)
    ###' - Constrained Bayes estimates (CB)
    ###' - Triple goal estimates (GR)
    ###' - Posterior means of ranks (rbar)
    ###' - Integer ranks (rhat)
    ###'
    ###' and add the TRUE site-specific values etc. 
    ###' 
    ###'
    
    ### Compute PM, PSD, CB, GR, rbar, and rhat
    df_theta <- df_posterior %>%
      dplyr::select(contains("tau_k["))  # Extract only site-specific estimates
    
    df_temp <- HETOP::triple_goal(as.matrix(df_theta))
    
    
    ### Combine true theta and ML theta (maximum likelihood) with true ranks
    df_est <- cbind.data.frame(df_meta, df_temp) %>%
      dplyr::select(-index) %>%
      rename(theta_true = theta, theta_ml = Y) %>%
      mutate(rtrue = rank(theta_true), 
             rpm = rank(theta_pm),
             rcb = rank(theta_cb), 
             rgr = rank(theta_gr), 
             rml = rank(theta_ml))
    
    
    ###'######################################################################
    ###'
    ###' Apply a loss function
    ###' 
    ###'  - Get tail quantile estimates
    ###'
    ###' 
    
    tail <- get_tail_quantiles(df_est, tail_p, 
                               priquan, est_vec = post_sum_methods)
    
    sim_tail <- rbind(sim_tail, tail)
  
  } # End of loops over i (Monte Carlo iteration)
  
  
  
  ###'######################################################################
  ###'
  ###' Tidy up and save the results from one condition (j)
  ###' 
  ###' and save the results in the defined folders
  ###'
  ###'
  
  setwd(folder_dir)
  
  
  
  ###'######################################################################
  ###'
  ###' Tidy up and save the 100 replications of tail quantiles
  ###'
  ###'
  
  ### Generate a column name vector
  vec_est <- rep(post_sum_methods, times = length(tail_p))
  vec_tail <- rep(tail_p, each = J_max)
  
  vec_colname <- data.frame(vec_est, vec_tail) %>%
    mutate(vec_tail = sprintf("%02d", vec_tail*100)) %>%
    unite(varname, vec_est, vec_tail, sep = "_") %>%
    .$varname 
  
  
  ### Convert sim_tail matrix into the dataframe (long format)
  df_tail <- sim_tail %>%
    data.frame() %>%
    set_names(vec_colname) %>%
    mutate(iter = seq(N_sim)) %>%
    dplyr::select(iter, everything()) %>%
    gather(key = key, value = value, -iter) %>%
    separate(key, into = c("estimator", "quantile"), sep = "_") %>%
    dplyr::select(estimator, quantile, iter, value) %>%
    mutate(estimator = factor(estimator, levels = post_sum_methods), 
           quantile = as.numeric(quantile)*0.01, 
           value = value) %>%
    arrange(estimator, quantile, iter) %>%
    mutate(diff = value - quantile) %>%
    set_names(c("estimator", "tail_p", "iter", "quantile", "diff"))
  
  
  ### Save the resulting dataframe
  write.csv(df_tail, file = "10_Percentile estimates of G_with 100 reps.csv")
  
  
  
  ###'######################################################################
  ###'
  ###' Print progress
  ###' 
  ###' 
  
  cat(paste0(sprintf("%03d", j)), 
      gen_name, anl_name, paste0("N = ", K_max), 
      paste0("ICC = ", ICC), paste0("rsr = ", rsr), "\n")
  
}  # End of loops over j (simulation condition)



###'######################################################################
###'
###' Measure computation time
###'
###'

### Mark end time
end_time <- Sys.time()

difftime(end_time, start_time)
