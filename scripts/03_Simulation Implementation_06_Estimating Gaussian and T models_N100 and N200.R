
###'######################################################################
###'
###' Category: Simulation Implementation
###' 
###' Task: Estimating Gaussian & T models with N = 100, 200  
###'       
###' Data: Simulated data
###' 
###' Data: 2020-04-07
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
library(parallel)
library(foreach)
library(doAzureParallel)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "functions/03_Simulation Implementation_04_Plot helpers.R")
source(file = "functions/03_Simulation Implementation_05_Data management helpers.R")



###'######################################################################
###'
###' Generate a configuration file called "credentials.json"
###' 
###'  in the working directory 
###'
###' - And populate this file with 
###'   my Batch and storage account names and keys
###'
###'

### Generate a configuration file
setwd(work_dir)
generateCredentialsConfig("credentials.json")


### Set the credentials for my current R session
setCredentials("credentials.json")



###'######################################################################
###'
###' Create a Batch pool    
###'
###'

### Generate a cluster configuration JSON file in the working directory
setwd(work_dir)
generateClusterConfig("cluster_Main_estimation.json")


### Create your cluster if it does not exist; this takes a few minutes
cluster <- makeCluster("cluster_Main_estimation.json") 


### Register your parallel backend 
registerDoAzureParallel(cluster) 


### Check that the nodes are running 
getDoParWorkers()



###'######################################################################
###'
###' Load the simulation variables
###'
###'

df_sim <- read.csv(file = "tables/simulation parameters_2_All with dg_name mark.csv") %>%
  select(-X)



###'######################################################################
###'
###' Filter only Gaussian and T data-analytic models
###' 
###' with N = 100 and 200
###' 
###' which is computationally intensive
###'
###'

tabdf(df_sim, assumed)

df_sim <- df_sim %>%
  filter(truth %in% c("Mixed")) %>%
  filter(assumed %in% c("Gaussian")) %>%
  filter(N %in% c(100, 200)) %>%



###'######################################################################
###'
###' Start loops over simulation conditions
###'
###'

### The number of simulation conditions 
N_cond <- length(df_sim$truth)   


### Number of replication
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
  
  ###' Optimize runtime. 
  ###' Chunking allows running multiple iterations on a single R instance.
  opt <- list(chunkSize = 4) 
  
  
  list_collect <- foreach (i = seq(N_sim), 
                           .options.azure = opt, 
                           .errorhandling = "pass") %dopar% {
                             
                             ### Call libraries within the active nodes
                             library(tidyverse)
                             library(DPpackage)
                             library(hhsim)
                             
                             
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
                               
                               
                             } else if (anl_name == "Model_DP-diffuse"){
                               
                               ### Prepare dataset as a matrix form (with y and sigma2)
                               mat_DPmeta <- df_meta %>% select(Y, sigma2) %>% as.matrix()
                               
                               
                               ### Set MCMC parameters
                               state <- NULL    # the current value of the parameters
                               mcmc <- list(nburn = 4000,    # the number of burn-in scans, 
                                            nsave = 4000,    # the total number of scans to be saved,
                                            nskip = 20,      # the thinning interval, 
                                            ndisplay = 100)  # the number of saved scans to be displayed on screen
                               
                               
                               ### Diffuse prior with respect to \alpha (the precision parameter)
                               b <- 0.1
                               alpha_mean <- K_max/2
                               a <- alpha_mean*b
                               alpha_var <- a/b^2  
                               
                               prior <- list(a0 = a,      # alpha0: shape param.
                                             b0 = b,      # alpha0: rate param.  
                                             tau1 = 1,    # G0 variance: shape param.
                                             tau2 = 1,    # G0 variance: rate param. 
                                             mub = 0, 
                                             Sb = 100)    # G0 mean: mean 
                               
                               
                               ### Fit the model with the DP prior
                               outp <- DPmeta(formula = mat_DPmeta ~ 1,
                                              prior = prior, mcmc = mcmc, 
                                              state = state, status = TRUE)
                               
                               
                               ### Extract posterior samples
                               df_posterior <- get_posterior_DPmeta(outp, nburn_DP = 4000)
                               
                             }
                             
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


