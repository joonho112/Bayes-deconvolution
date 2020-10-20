
###'######################################################################
###'
###' Category: Code Experiment
###' 
###' Task: HETOP
###'       
###'       (1) Try fh_hetop() function
###'       
###' Data: Simulated data
###' 
###' Data: 2020-03-07
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
library(HETOP)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)


### Set ggplot themes
### Theme settings
theme_preset <- 
  theme_bw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_blank())



###'######################################################################
###' 
###' Generate data
###'
###'

set.seed(1001)

### Define mean-centered covariates
G  <- 12
z1 <- sample(c(0, 1), size = G, replace = TRUE)
z2 <- 0.5*z1 + rnorm(G)
Z  <- cbind(z1 - mean(z1), z2 = z2 - mean(z2))


### Define true parameters dependent on covariates
beta_m    <- c(0.3,  0.8)
beta_s    <- c(0.1, -0.1)
mug       <- Z[, 1]*beta_m[1] + Z[, 2]*beta_m[2] + rnorm(G, sd = 0.3)
sigmag    <- exp(0.3 + Z[,1]*beta_s[1] + Z[,2]*beta_s[2] + 0.2*rt(G, df = 7))
cutpoints <- c(-1.0, 0.0, 1.2)


### Generate data
ng   <- rep(200, G)
ngk  <- gendata_hetop(G, K = 4, ng, mug, sigmag, cutpoints)
print(ngk)



###'######################################################################
###' 
###' Fit FH-HETOP model including covariates  
###' 
###' NOTE: using an extremely small number of iterations for testing, 
###'       so that convergence is not expected
###' 
###' 

### Fit FH-HETOP model
m <- fh_hetop(ngk,  # data 
              fixedcuts = c(-1.0, 0.0),  # the first two cutpoints
              
              # Arguments needed to parameterize Efron prior
              p = c(10,10),   # degrees of freedom for cubic spline basis  
              m = c(100, 100),  # number of grid points  
              gridL = c(-5.0, log(0.10)),  # lower bounds for grids
              gridU = c(5.0, log(5.0)),    # upper bounds for grids
              
              # Covariates
              Xm = Z,  
              Xs = Z,
              
              # MCMC options
              n.iter = 100, 
              n.burnin = 50)



###'######################################################################
###'
###' Print results
###'
###'

### Entire summary
print(m)


### Look into additional information stored as a list
print(names(m$fh_hetop_extras))
m$fh_hetop_extras

s <- m$BUGSoutput$summary
print(data.frame(truth = c(beta_m, beta_s), s[grep("beta", rownames(s)),]))

print(cor(mug,    s[grep("mu",    rownames(s)),"mean"]))
print(cor(sigmag, s[grep("sigma", rownames(s)),"mean"]))


### Manual calculation of WAIC (see help file for waic_hetop)
tmp <- waic_hetop(ngk, m$BUGSoutput$sims.matrix)
identical(tmp, m$fh_hetop_extras$waicinfo)




