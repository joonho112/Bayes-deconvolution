
###'######################################################################
###'
###' Category: Simulation Implementation
###' 
###' Task: Pre-generate and save 100 random samples 
###'       per each simulation condition
###'       
###' Data: Simulated data
###' 
###' Data: 2020-04-05
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
library(cowplot)
library(sn)
library(LaplacesDemon)
library(hhsim)
library(DPpackage)
library(HETOP)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "functions/03_Simulation Implementation_04_Plot helpers.R")
source(file = "functions/03_Simulation Implementation_05_Data management helpers.R")



###'######################################################################
###'
###' Call simulation variables - Only data specific ones
###'
###'

setwd(work_dir)

df_temp <- read.csv(file = "tables/simulation parameters_2_data-generation only.csv")

sim_parm <- df_temp %>%
  rename(j_cond = X) %>%
  unite(fname, truth, N, ICC, rsr, sep = "_") %>%
  dplyr::select(-df) %>% 
  cbind.data.frame(df_temp) %>%
  dplyr::select(j_cond, truth, N, ICC, rsr, df, fname)

write.csv(sim_parm, file = "tables/simulation parameter_2_data-generation only_fname.csv")



###'######################################################################
###'
###' Start loops over simulation conditions
###'
###'

### The number of simulation conditions 
N_cond <- nrow(sim_parm)   


### The number of Monte Carlo iterations to perform
N_sim <- 100


### Mark start time
start_time <- Sys.time()



for (j in 1:N_cond){
  
  ###'######################################################################
  ###'
  ###' Set required objects for simulation
  ###' 
  ###'
  
  ### Extract simulation condition elements
  gen_name <- paste0("True_", sim_parm$truth[j])
  
  N <- sim_parm$N[j]
  
  ICC <- sim_parm$ICC[j]
  
  rsr <- sim_parm$rsr[j]
  
  nu <- sim_parm$df[j]
  
  fname <- sim_parm$fname[j]
  
  
  ### Set quantiles (tail probabilities)
  tail_p <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95) 
  
  
  
  ###'######################################################################
  ###'
  ###' Perform Mote Carlo iterations (N_sim)
  ###'
  ###'
  
  ### Prepare an empty list to save the generated data
  list_collect <- list()
  
  
  for (i in seq(N_sim)){
    
    ### Print progress
    cat(paste0(sprintf("%03d", j)), 
        gen_name, paste0("N = ", N), 
        paste0("ICC = ", ICC), paste0("rsr = ", rsr), 
        paste0("i = ", i), "\n")
  
    
    
    ###'######################################################################
    ###'
    ###'  Generate data 
    ###' 
    ###' - Conforming to the "true" distribution
    ###'
    ###'
    
    if (gen_name == "True_Gaussian"){ 
      
      list_meta <- Gen_Gaussian(n = N, ICC = ICC, rsr = rsr)
      df_meta <- list_meta %>% list_to_df_G()
      priquan <- Q2CDF_Gaussian(tail_p)
      
      
    } else if (gen_name == "True_T"){   
      
      nu <- sim_parm$df[j]  # degrees of freedom
      
      list_meta <- Gen_T(n = N, ICC = ICC, rsr = rsr, nu = nu)
      df_meta <- list_meta %>% list_to_df_G()
      priquan <- Q2CDF_T(tail_p)
      
      
    } else if (gen_name == "True_ALD"){
      
      p <- 0.1   # control the skewness parameter
      
      list_meta <- Gen_ALD(n = N, ICC = ICC, rsr = rsr, 
                           mean = 0, var = 1, p = p)
      
      df_meta <- list_meta %>% list_to_df_G()
      priquan <- Q2CDF_Mixture(tail_p)
      
      
    } else if (gen_name == "True_Bimodal"){
      
      delta <- 4   # distance between two means
      eps <- 0.3   # proportion of the small component
      ups <- 1     # ratio between two variances
      
      list_meta <- Gen_Mixture(n = N, ICC = ICC, rsr = rsr, 
                               delta = delta, eps = eps, ups = ups)
      
      df_meta <- list_meta %>% list_to_df_G()
      priquan <- Q2CDF_Mixture(tail_p)
      
      
    } else if (gen_name == "True_Skew"){
      
      slant <- 10   # control the degree of skewness
      
      list_meta <- Gen_SkewN(n = N, ICC = ICC, rsr = rsr, 
                             mean = 0, var = 1, slant = slant)
      
      df_meta <- list_meta %>% list_to_df_G()
      priquan <- Q2CDF_Mixture(tail_p)
      
      
    } else if (gen_name == "True_Mixed"){
      
      delta <- 5   # distance between two means
      eps <- 0.3   # proportion of the small component
      ups <- 2     # ratio between two variances
      
      list_meta <- Gen_Mixture(n = N, ICC = ICC, rsr = rsr, 
                               delta = delta, eps = eps, ups = ups)
      
      df_meta <- list_meta %>% list_to_df_G()
      priquan <- Q2CDF_Mixture(tail_p)
      
    }
    
    
    
    ###'######################################################################
    ###'
    ###' Save the 100 random samples 
    ###'
    ###'
    
    list_temp <- list(list_meta, df_meta, priquan, c(var(list_meta$theta)))
    
    names(list_temp) <- c("list_data", "df_data", "CDF_theory", "var_theta")
    
    list_collect[[i]] <- list_temp
    
    
  } # End of loops over i (Monte Carlo iteration)
  
  
  ###'######################################################################
  ###'
  ###' Save the generated data for each data-generating conditions
  ###'
  ###'

  setwd(work_dir)
  
  saveRDS(list_collect, file = file.path("datasets",
                                         "05_Pre-generated simulated data", 
                                         paste0(fname, ".rds")))

  
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


