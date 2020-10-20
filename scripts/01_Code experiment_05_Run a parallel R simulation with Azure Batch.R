
###'######################################################################
###'
###' Category: Code expriment
###' 
###' Task: Run a parallel R simulation with Azure Batch
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
library(foreach)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)



###'######################################################################
###'
###' Install doAzureParallel package
###'
###'

# Install the devtools package  
install.packages("devtools") 

# Install rAzureBatch package
devtools::install_github("Azure/rAzureBatch") 

# Install the doAzureParallel package 
devtools::install_github("Azure/doAzureParallel") 

# Load the doAzureParallel library 
library(doAzureParallel)



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

# Generate a configuration file
setwd(work_dir)
generateCredentialsConfig("credentials.json")


# Set the credentials for my current R session
setCredentials("credentials.json")



###'######################################################################
###'
###' Create a Batch pool    
###'
###'

# Generate a cluster configuration JSON file in the working directory
setwd(work_dir)
generateClusterConfig("cluster.json")


# Create your cluster if it does not exist; this takes a few minutes
cluster <- makeCluster("cluster.json") 


# Register your parallel backend 
registerDoAzureParallel(cluster) 


# Check that the nodes are running 
getDoParWorkers()



###'######################################################################
###'
###' Prepare a parallel simulation
###'
###'

# Set parameters for the Monte Carlo simulation
mean_change <- 1.001 
volatility <- 0.01 
opening_price <- 100


# Define a function to simulate closing prices
getClosingPrice <- function() { 

  days <- 1825 # ~ 5 years 
  movement <- rnorm(days, mean = mean_change, sd = volatility) 
  path <- cumprod(c(opening_price, movement)) 
  closingPrice <- path[days] 
  return(closingPrice) 

}



###'######################################################################
###'
###' Run 10,000 simulations locally using a standard foreach loop
###' 
###' with %do%  
###'
###'

# Run 10,000 simulations in series 
start_s <- Sys.time() 

closingPrices_s <- foreach(i = 1:10, .combine='c') %do% { 
  replicate(1000, getClosingPrice()) 
} 

end_s <- Sys.time()


# Plot the closing prices in a histogram
hist(closingPrices_s)


# Estimate computation time
difftime(end_s, start_s)
1000*difftime(end_s, start_s, unit = "min")



###'######################################################################
###'
###' Run the code using foreach with %dopar%
###'
###'

# Optimize runtime. 
# Chunking allows running multiple iterations on a single R instance.
opt <- list(chunkSize = 10) 


# Run simulation
start_p <- Sys.time()  

closingPrices_p <- foreach(i = 1:100, .combine='c', .options.azure = opt) %dopar% { 
  replicate(100000, getClosingPrice()) 
} 

end_p <- Sys.time()


# Plot the closing prices in a histogram
hist(closingPrices_p)


# Estimate computation time
difftime(end_p, start_p)
1000*difftime(end_s, start_s, unit = "min")



###'######################################################################
###'
###' Clean up resources
###'
###'

stopCluster(cluster)



