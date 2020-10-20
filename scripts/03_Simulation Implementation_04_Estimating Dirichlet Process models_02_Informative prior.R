
###'######################################################################
###'
###' Category: Simulation Implementation
###' 
###' Task: Estimating Dirichlet Process models  
###' 
###'      (2) Informative Priors
###' 
###'     - Flexible distributions for triple-goal estimates in 
###'       two-stage hierarchical models
###'       
###'     - Using Dirichlet Processes for Modeling Heterogenous 
###'       Treatment Effects Across Sites  
###'       
###' Data: Simulated data
###' 
###' Data: 2020-04-18
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
generateCredentialsConfig("cluster/credentials.json")


### Set the credentials for my current R session
setCredentials("cluster/credentials.json")



###'######################################################################
###'
###' Create a Batch pool    
###'
###'

### Generate a cluster configuration JSON file in the working directory
setwd(work_dir)
generateClusterConfig("cluster/cluster_Main_estimation.json")


### Create your cluster if it does not exist; this takes a few minutes
cluster <- makeCluster("cluster/cluster_Main_estimation.json") 


### Register your parallel backend 
registerDoAzureParallel(cluster) 


### Check that the nodes are running 
getDoParWorkers()



###'######################################################################
###'
###' Prepare simulation variables
###' 
###' (1) N: number of sites (K = 25, 50, 100, 200)
###'
###' (2) ICC: Intraclass correlation (ICC = 0.10, 0.50, 0.90)
###'    - ICC = 1/(1 + rgt)
###'    - rgt: the informativeness of the data (rgt = 9, 1.0, 0.111)
###'    - Larger ICC = More information within sites
###' 
###' (3) rsr: the heterogeneity of the uncertainty estimates (rsr = 1^2, 5^2, 10^2)
###' 
###' (4) True population distribution of G
###'     All DGMs have zero mean and variance of one
###'     - 1. Gaussian distribution N(0, 1)
###'     - 2. T distribution T_2(0, 1) with df = 5  => Fatter tail
###'     - 3. ALD distribution with mean 0 and variance 1 
###'     - 4. Gaussian mixture 0.70 N(0, 1) + 0.30 (4, 1) with same variance
###'     
###'     [Optional]
###'     - 5. Skewed normal (tilt) with slant parameter of 5
###'     - 6. Gaussian mixture with different variances 
###'
###'  

### Set seed 
set.seed(10101)


### Generate a table of the simulation parameters
tbl_parm <- expand.grid(truth = c("Gaussian", "T", "ALD", "Bimodal", "Skew", "Mixed"),
                        assumed = c("Gaussian", "T", "DP-diffuse", "DP-EB"),
                        ICC = c(0.10, 0.50, 0.90),
                        rsr = c(1, 5, 10),
                        N = c(25, 50, 100, 200), 
                        df = 5)

nrow(tbl_parm)


### Reorder the table & add a variable
tbl_parm <- tbl_parm %>% 
  dplyr::select(truth, N, ICC, rsr, assumed, df) %>%
  arrange(truth, N, ICC, rsr, assumed, df)


### Add a concatenated design indicator
sim_parm <- tbl_parm %>%
  unite(fname, truth, N, ICC, rsr, assumed, sep = "_") %>%
  dplyr::select(-df) %>%
  cbind.data.frame(tbl_parm) %>%
  select(truth:df, fname)

head(sim_parm, 10)


### Save as a table 
write.csv(sim_parm, file = "tables/simulation parameters_2.csv")


### The number of Monte Carlo iterations to perform
N_sim <- 100



###'######################################################################
###'
###' Prepare data simulation variable
###' 
###' - Which excludes only the data-analytic model variable
###'
###'

### Extract only simulation variables related to the data-generation
data_parm_sub <- sim_parm %>%
  dplyr::select(truth, N, ICC, rsr, df) %>%
  distinct()

nrow(data_parm_sub)

write.csv(data_parm_sub, 
          file = "tables/simulation parameters_2_data-generation only.csv")


### Generate ID for each data-generating condition
data_parm <- data_parm_sub %>%
  unite(dg_name, truth, N, ICC, rsr, sep = "_") %>%
  dplyr::select(-df) %>% 
  cbind.data.frame(data_parm_sub) %>%
  mutate(j_cond = 1:nrow(data_parm_sub)) %>%
  dplyr::select(j_cond, truth, N, ICC, rsr, df, dg_name)

write.csv(data_parm, 
          file = "tables/simulation parameter_2_data-generation only_dg_name.csv")



###'######################################################################
###'
###' Merge data_parm with sim_parm to mark data-generating models
###'
###'

df_sim <- sim_parm %>%
  full_join(data_parm, by = c("truth", "N", "ICC", "rsr", "df")) %>%
  dplyr::select(-j_cond)



###'######################################################################
###'
###' Filter only DP-EB model
###' 
###' which is (or could be) computationally intensive
###'
###'

tabdf(df_sim, assumed)

df_sim <- df_sim %>%
  filter(assumed == "DP-EB")



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
  
  ###' Optimize runtime. 
  ###' Chunking allows running multiple iterations on a single R instance.
  opt <- list(chunkSize = 4) 
  
  
  list_collect <- 
    foreach (i = seq(N_sim), .options.azure = opt, .errorhandling = "pass") %dopar% {
             
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
       
     } else if (anl_name == "Model_DP-EB"){
       
       ### Prepare dataset as a matrix form (with y and sigma2)
       mat_DPmeta <- df_meta %>% select(Y, sigma2) %>% as.matrix()
       
       
       ### Set MCMC parameters
       state <- NULL    # the current value of the parameters
       mcmc <- list(nburn = 4000,    # the number of burn-in scans, 
                    nsave = 4000,    # the total number of scans to be saved,
                    nskip = 20,      # the thinning interval, 
                    ndisplay = 100)  # the number of saved scans to be displayed on screen
       
       
       ### Diffuse prior with respect to \alpha (the precision parameter)
       info_priors <- list(c(1.24, 0.64), 
                           c(1.60, 1.22), 
                           c(1.84, 1.82),
                           c(1.96, 2.38))
       
       names(info_priors) <- c("25", "50", "100", "200")
       
       ab_info <- info_priors[[as.character(K_max)]]
       
       a <- ab_info[1]; b <- ab_info[2]
       
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
      paste0("ICC = ", ICC), paste0("rsr = ", rsr), "\n")
  
  
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


