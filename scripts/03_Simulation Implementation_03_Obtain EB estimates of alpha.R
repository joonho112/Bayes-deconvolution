
###'######################################################################
###'
###' Category: Simulation Implementation
###' 
###' Task: Obtain Empirical Bayes estimates of the precision parameter
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
data_dir <- file.path(work_dir, "datasets", "05_Pre-generated simulated data")


### Call libraries
library(tidyverse)
library(cowplot)
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
###' Call simulation variables - Only data specific ones
###'
###'

setwd(work_dir)

sim_parm <- read.csv(file = "tables/simulation parameter_2_data-generation only_fname.csv")

sim_parm <- sim_parm %>%
  select(-X)



###'######################################################################
###'
###' Create a Batch pool    
###'
###'

### Generate a cluster configuration JSON file in the working directory
setwd(work_dir)
generateClusterConfig("cluster_EB_estimate.json")


### Create your cluster if it does not exist; this takes a few minutes
cluster <- makeCluster("cluster_EB_estimate.json") 


### Register your parallel backend 
registerDoAzureParallel(cluster) 


### Check that the nodes are running 
getDoParWorkers()



###'######################################################################
###'
###' Start loops over simulation conditions
###'
###'

### Set data containing working directory 
setwd(data_dir)


### The number of simulation conditions 
N_cond <- nrow(sim_parm)   


### The number of Monte Carlo iterations to perform
N_sim <- 100


### Mark start time
start_time <- Sys.time()


for (j in 5:N_cond){
  

  ###'######################################################################
  ###'
  ###' Extract simulation condition elements
  ###' 
  ###'

  gen_name <- paste0(sim_parm$truth[j])
  
  N <- sim_parm$N[j]
  
  ICC <- sim_parm$ICC[j]
  
  rsr <- sim_parm$rsr[j]
  
  nu <- sim_parm$df[j]
  
  fname <- sim_parm$fname[j]
  
  
  ###'######################################################################
  ###'
  ###' Load the pre-generated 100 random samples
  ###'
  ###'
  
  setwd(data_dir)
  list_data <- readRDS(file = paste0(sim_parm$fname[j], ".rds"))
  
  
  ###'######################################################################
  ###'
  ###' Loop over 100 random samples
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
    
    
    ### Print progress
    cat(paste0(sprintf("%03d", j)), 
        gen_name, paste0("N = ", N), 
        paste0("ICC = ", ICC), paste0("rsr = ", rsr), 
        paste0("i = ", i), "\n")
    
    
    ### Extract the simulated data
    df_G <- list_data[[i]]$df_data
    
    
    ###'######################################################################
    ###'
    ###' Prepare DP model fitting
    ###'
    ###'
    
    ### Set initial state (the current value of the parameters)
    state <- NULL
    
    
    ### Set MCMC parameters
    nburn <- 2000    # the number of burn-in scans
    nsave <- 2000    # the total number of scans to be saved
    nskip <- 20      # the thinning interval
    ndisplay <- 100  # the number of saved scans to be displayed on screen
    
    mcmc <- list(nburn = nburn, nsave = nsave,
                 nskip = nskip, ndisplay = ndisplay)
    
    
    ### Prepare dataset as a matrix form (with y and sigma2)
    mat_DPmeta <- df_G %>% select(Y, sigma2) %>% as.matrix()
    
    
  
    ###'######################################################################
    ###'
    ###' Construct a `repeat` loop for the alpha convergence
    ###'
    ###'
    
    ### Set the initial value of alpha and tolerance level
    init_range <- c(1/log(N), N/log(N))
    alpha0 <- (init_range[1] + init_range[2])/2
    tol <- 0.01
    # tol <- 1e-3
    ir <- 1   # Loop indicator
    vec_alpha <- vector()   # allocate the space to store alpha values
    
    
    ### Start a repeat loop
    repeat {
      
      ## Set prior parameters
      prior <- list(alpha = alpha0,  # fixed precision parameter
                    tau1 = 1,   # G0 variance: shape param.
                    tau2 = 1,   # G0 variance: rate param. 
                    mub = 0, 
                    Sb = 100)   # G0 mean: mean 
      
      
      ## Estimate the Dirichlet Process model using DPpackage 
      outp <- DPmeta(formula = mat_DPmeta ~ 1,
                     prior = prior, mcmc = mcmc, 
                     state = state, status = TRUE)
      
      
      ### Get the posterior samples
      df_posterior <- get_posterior_DPmeta(outp, nburn_DP = nburn)
      
      
      ###' Compute the value of alpha that satisfies 
      ###' Kbar = sum_{i=1}^{N}(alpha/(alpha + i - 1))
      K_postmean <- df_posterior$Ncluster %>%
        mean(na.rm = TRUE)
      
      i_vec <- seq(from = 1, to = N, by = 1)
      
      fn_K <- function(alpha, N, K_postmean){
        
        seq(from = 1, to = N, by = 1)
        abs(K_postmean - sum(alpha/(alpha + i_vec - 1))) 
        
      }
      
      solution <- optimize(fn_K, c(0, 100), N = N, K_postmean, tol = 0.0000001)
      
      vec_alpha[ir] <- alpha1 <- solution$minimum
      
      
      ### Test convergence: Close enough? or exceed 50 iter?
      
      if(abs(alpha1 - alpha0) < tol | ir > 100) {  
        break
      } else {
        
        # Assign new values
        ir <- ir + 1
        alpha0 <- alpha1
        
        # Print progress
        cat(paste0("ir = ", ir, ", alpha = ", alpha0), "\n")
      } 
    }
    
    ### Tidy up and save the alpha vector
    df_alpha <- data.frame(vec_alpha) %>%
      rownames_to_column() %>% 
      rename(iteration = rowname, alpha = vec_alpha) %>%
      mutate(iteration = as.numeric(iteration))
    
    df_alpha
    
  } # End of loops over i (Monte Carlo iteration)
  
  
  ###'######################################################################
  ###'
  ###' Save the generated dataframes for EB estimates
  ###'
  ###'
  
  setwd(work_dir)
  
  saveRDS(list_collect, file = file.path("datasets",
                                         "06_EB estimates of alpha", 
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



