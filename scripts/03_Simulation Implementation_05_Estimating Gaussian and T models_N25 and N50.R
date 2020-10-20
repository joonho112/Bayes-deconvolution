
###'######################################################################
###'
###' Category: Simulation Implementation
###' 
###' Task: Estimating Gaussian and T models  
###' 
###' Data: Simulated data
###' 
###' Data: 2020-04-06
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
gc()            # force R to release memory it is no longer using
rm(list=ls())   # delete all the objects in the workspace


### Set working directory 
work_dir <- c("~/Bayes-deconvolution")
setwd(work_dir)


### Set a data directory
data_dir <- file.path(work_dir, "datasets")


### Call libraries
library(tidyverse)
library(hhsim)
library(DPpackage)
library(HETOP)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "functions/03_Simulation Implementation_04_Plot helpers.R")
source(file = "functions/03_Simulation Implementation_05_Data management helpers.R")



###'######################################################################
###'
###' Load and save the table of simulation variables
###'
###'

### Load the simulation parameters: ALL variables
setwd(work_dir)

sim_parm <- read.csv(file = "tables/simulation parameters_2.csv") %>%
  dplyr::select(-X)


### Load the simulation parameters: Only data-generating variables
data_parm <- read.csv(file = "tables/simulation parameter_2_data-generation only_dg_name.csv") %>%
  dplyr::select(-X)


### Merge data_parm with sim_parm to mark data-generating models
df_sim <- sim_parm %>%
  full_join(data_parm, by = c("truth", "N", "ICC", "rsr", "df")) %>%
  dplyr::select(-j_cond)

write.csv(df_sim, file = "tables/simulation parameters_2_All with dg_name mark.csv")



###'######################################################################
###'
###' Filter only Gaussian and T models
###' 
###' which are NOT computationally intensive
###'
###'

tabdf(df_sim, assumed)

df_sim <- df_sim %>%
  filter(assumed %in% c("Gaussian", "T")) %>%
  filter(N %in% c(25, 50)) %>%
  filter(truth != c("Gaussian"))



###'######################################################################
###'
###' Start loops over simulation conditions
###'
###'

### The number of simulation conditions 
N_cond <- length(df_sim$truth)   


### Mark start time
start_time <- Sys.time()


for (j in 1:N_cond){
  
  
  ###'######################################################################
  ###'
  ###' Set required objects for simulation
  ###' 
  ###'
  
  ### Extract simulation condition elements
  gen_name <- paste0("True_", df_sim$truth[j])
  anl_name <- paste0("Model_", df_sim$assumed[j])
  K_max <- df_sim$N[j]
  ICC <- df_sim$ICC[j]
  rsr <- df_sim$rsr[j]
  nu <- df_sim$df[j]
  fname <- df_sim$fname[j]
  dg_name <- df_sim$dg_name[j]
  
  
  ### Define necessary parameters to control Bayesian model estimation
  nmcmc <- 2000    # number of MCMC iterations per each model fit
  nburn <- 2000    # burn-in iterations to discard at each round
  nDraws <- nmcmc + nburn
  
  
  
  ###'######################################################################
  ###'
  ###' Load pre-generated 100 random samples
  ###'
  ###'
  
  setwd(file.path(data_dir, "05_Pre-generated simulated data"))
  
  list_data <- readRDS(file = paste0(dg_name, ".rds"))
  
  setwd(work_dir)
  
  
  
  ###'######################################################################
  ###'
  ###' Perform Mote Carlo iterations (N_sim)
  ###'
  ###'
  
  ### The number of Monte Carlo iterations to perform
  N_sim <- 100
  
  
  ### Prepare an empty list to collect posterior samples
  list_collect <- list()
  
  
  for (i in seq(N_sim)){
    
    ###'######################################################################
    ###'
    ###' Extract "i"th list data
    ###'
    ###'
    
    list_meta <- list_data[[i]]$list_data
    
    df_meta <- list_data[[i]]$df_data
    
    
    
    ###'######################################################################
    ###'
    ###' Fit the assumed model 
    ###' 
    ###' - To get posterior samples 
    ###'
    ###' 
    
    if (anl_name == "Model_Gaussian"){
      
      ### Fit the model
      outp <- GaussianGaussian(list_meta$Y,
                               list_meta$sigma2,
                               nDraws)
      
      ### Extract posterior samples
      df_posterior <- get_posterior_hhsim(output = outp, nburn, nDraws)
      
      
    } else if (anl_name == "Model_T"){
      
      ### Fit the model
      outp <- GaussianT(list_meta$Y,
                        list_meta$sigma2,
                        nDraws, nu)
      
      ### Extract posterior samples
      df_posterior <- get_posterior_hhsim(output = outp, nburn, nDraws)
      
    } 
    
    
    ###'######################################################################
    ###'
    ###' Embed into the prepared list
    ###'
    ###'
    
    list_collect[[i]] <- df_posterior
    
    
  } # End of loops over i (Monte Carlo iteration)
  
  
  
  ###'######################################################################
  ###'
  ###' Tidy up and save the results from one condition (j)
  ###'
  ###'
  
  ### Create a directory to save the results
  folder_dir <- file.path(work_dir, "datasets", 
                          "07_Simulation results-Final", 
                          fname)
  
  dir.create(folder_dir, showWarnings = FALSE)
  
  
  ### Save the resulting 100 posterior samples
  setwd(folder_dir)
  
  saveRDS(list_collect, file = "100 posterior samples.rds")
  
  setwd(work_dir)
  
  
  ### Print progress
  cat(paste0(sprintf("%03d", j)), 
      gen_name, paste0("N = ", K_max), 
      paste0("ICC = ", ICC), paste0("rsr = ", rsr), 
      anl_name, "\n")
  
  
}  # End of loops over j (simulation condition)



###'######################################################################
###'
###' Measure computation time
###'
###'

### Mark end time
end_time <- Sys.time()

comp_time <- end_time - start_time

comp_time


