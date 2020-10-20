
###'######################################################################
###'
###' Category: Code experiment
###' 
###' Task: Run a parallel R simulation with Azure Batch
###'       (2) with hhsim and DPpackage
###'
###' Date: 2020-04-03
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
library(parallel)
library(foreach)
library(doAzureParallel)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "scripts/03_Simulation Implementation_04_Plot helpers.R")
source(file = "scripts/03_Simulation Implementation_05_Data management helpers.R")



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
generateClusterConfig("cluster.json")


### Create your cluster if it does not exist; this takes a few minutes
cluster <- makeCluster("cluster.json") 

# resizeCluster(cluster, 
#               dedicatedMin = 0, 
#               dedicatedMax = 0, 
#               lowPriorityMin = 5, 
#               lowPriorityMax = 10)


### Register your parallel backend 
registerDoAzureParallel(cluster) 


### Check that the nodes are running 
getDoParWorkers()


###' Optimize runtime. 
###' Chunking allows running multiple iterations on a single R instance.
opt <- list(chunkSize = 10) 



###'######################################################################
###'
###' Prepare a parallel simulation
###'
###'

### Load functions for data generation & converting quantiles to CDF
setwd(work_dir)
source(file = "scripts/03_Simulation Implementation_02_Data generation.R")
source(file = "scripts/03_Simulation Implementation_03_Convert quantiles to theoretical CDFs.R")


### Set data-generating parameters outside the foreach loop
K_max <- 50
rgt <- 10
rsr <- 10
nu <- 5

delta <- 4   # distance between two means
eps <- 0.2   # proportion of the small component
ups <- 1     # ratio between two variances

tail_p <- c(0.05, 0.10, 0.25, 0.75, 0.90, 0.95) 


### Set MCMC parameters
state <- NULL    # the current value of the parameters
mcmc <- list(nburn = 4000,    # the number of burn-in scans, 
             nsave = 4000,    # the total number of scans to be saved,
             nskip = 20,      # the thinning interval, 
             ndisplay = 100)  # the number of saved scans to be displayed on screen


nmcmc <- 2000    # number of MCMC iterations per each model fit
nburn <- 2000    # burn-in iterations to discard at each round
nDraws <- nmcmc + nburn


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



###'######################################################################
###'
###' Implement the parallel simulation sequentially with %do%
###' 
###' 85 seconds per model * 100 = 8500 seconds = 141 minutes
###'
###'

### Mark start time
start_time <- Sys.time()


### Implement a loop
temp_list <- foreach (i = 1:100, .options.azure = opt) %do% {
  
  ### Call libraries within the active nodes
  library(tidyverse)
  library(DPpackage)
  
  
  ### Generate data
  list_meta <- Generate_Mixture(n = K_max, rgt = rgt, rsr = rsr, 
                                delta = delta, eps = eps, ups = ups)
  
  df_meta <- list_meta %>% list_to_df_G()
  priquan <- Q2CDF_Mixture(tail_p)
  
  
  ### Prepare dataset as a matrix form (with y and sigma2)
  mat_DPmeta <- df_meta %>% select(Y, sigma2) %>% as.matrix()
  
  
  ### Fit the model with the DP prior
  outp <- DPmeta(formula = mat_DPmeta ~ 1,
                 prior = prior, mcmc = mcmc, 
                 state = state, status = TRUE)
  
  
  ### Extract posterior samples
  get_posterior_DPmeta(outp, nburn_DP = 4000)
  
}

### Mark end time
end_time <- Sys.time()

difftime(end_time, start_time, unit = "min")

comp_time <- end_time - start_time

comp_time



###'######################################################################
###'
###' Implement the parallel simulation in parallel with %dopar%
###' 
###' - Check points:
###'   (1) Do the GitHub-installed hhsim and DPpackage works well? 
###'      => DPpackage: YES
###'         hhsim: 
###'   
###'   (2) Do I need to define/assign local functions/objects 
###'       within the foreach loop? 
###'      => NO. It seems that it submits local subjects and functions
###'         to the cluster. 
###'
###' - Takes only 10 minutes << 141 minutes (when computed sequentially)
###'   14 times fater!!!
###'   
###'

### Mark start time
start_time <- Sys.time()


### Implement a loop
temp_list <- foreach (i = 1:100, .options.azure = opt) %dopar% {
  
  ### Call libraries within the active nodes
  library(tidyverse)
  library(DPpackage)
  library(hhsim)
  
  ### Generate data
  list_meta <- Generate_Mixture(n = K_max, rgt = rgt, rsr = rsr, 
                                delta = delta, eps = eps, ups = ups)
  
  df_meta <- list_meta %>% list_to_df_G()
  priquan <- Q2CDF_Mixture(tail_p)
  
  
  ### Prepare dataset as a matrix form (with y and sigma2)
  mat_DPmeta <- df_meta %>% select(Y, sigma2) %>% as.matrix()
  
  
  ### Fit the model with the DP prior
  outp <- DPmeta(formula = mat_DPmeta ~ 1,
                 prior = prior, mcmc = mcmc, 
                 state = state, status = TRUE)
  
  
  ### Extract posterior samples
  get_posterior_DPmeta(outp, nburn_DP = 4000)
  
}

### Mark end time
end_time <- Sys.time()

difftime(end_time, start_time, unit = "min")



###'######################################################################
###'
###' Clean up resources
###'
###'

stopCluster(cluster)




