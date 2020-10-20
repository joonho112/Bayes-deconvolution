
###'######################################################################
###'
###' Category: Code Validation
###' 
###' Task: Compare packages to implement the Empirical Bayes deconvolution
###'       
###' Data: Simulated data
###' 
###' Date: 2020-03-02
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
library(hhsim)
library(DPpackage)
library(bspmma)
library(deconvolveR)
library(ashr)
library(EbayesThresh)


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
###' Generate a simulated dataset 
###' 
###' from a mixture of two Gaussian components
###' 
###' (1) Get true PDF of thetas
###'
###'

### Set quantiles (tail probabilities)
tail_p <- c(0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95) 


### Simulation conditions
N <- 2000
ICC <- 0.9
rsr <- 5

delta <- 5   # distance between two means
eps <- 0.3   # proportion of the small component
ups <- 2     # ratio between two variances

list_meta <- Gen_Mixture(n = N, ICC = ICC, rsr = rsr, 
                         delta = delta, eps = eps, ups = ups)

df_meta <- list_meta %>% list_to_df_G()
priquan <- Q2CDF_Mixture(tail_p)



###'######################################################################
###'
###' Parameter estimation using the `deconvolveR` package 
###' 
###'
###'

### Prepare vectors for deconv fit
X_vec <- df_meta$theta
tau <- seq(from = -6, to = 6, by = 0.01)


### Fit the model deconvolveR
fit_deconv <- deconvolveR::deconv(
  tau = tau,            # a vector of discrete support points for theta 
  X = X_vec,            # a vector of sample values
  family = "Normal")      


### Plot estimation result
gData <- as.data.frame(fit_deconv$stats[, c("theta", "g")])

gData$g <- gData$g/sum(gData$g)

ggplot(data = gData) +
  geom_line(mapping = aes(x = theta, y = g)) +
  geom_line(mapping = aes(x = theta, y = dnorm(theta, mean = -3)),
            color = "red") +
  labs(x = expression(theta), y = expression(g(theta)))






g_vec <- fit_deconv$stats[, 2]


ggplot() + 
  geom_histogram(data = df_meta, aes(x = theta, y = ..count.. / sum(..count..)), 
                 color = "blue", fill = "red", alpha = 0.5) + 
  geom_histogram(data = df_meta, aes(x = Y, y = ..count.. / sum(..count..)), 
                 color = "brown", alpha = 0.5) + 
  geom_line(aes(x = tau, y = g_vec))




if (require("ggplot2")) {
  ggplot() +
    geom_histogram(mapping = aes(x = disjointTheta, y = ..count.. / sum(..count..) ),
                   color = "blue", fill = "red", bins = 40, alpha = 0.5) +
    geom_histogram(mapping = aes(x = z, y = ..count.. / sum(..count..) ),
                   color = "brown", bins = 40, alpha = 0.5) +
    geom_line(mapping = aes(x = tau, y = g), color = "black") +
    labs(x = paste(expression(theta), "and x"), y = paste(expression(g(theta)), " and f(x)"))
}







