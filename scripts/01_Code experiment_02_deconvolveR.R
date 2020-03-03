
###'######################################################################
###'
###' Category: Code Experiment
###' 
###' Task: deconvolveR - Empirical Bayes methods for learning prior distribution
###'       from data
###'       
###' Data: Simulated data
###' 
###' Data: 2020-03-02
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
library(hhsim)
library(DPpackage)
library(bspmma)
library(deconvolveR)


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
###' A simulated examole
###'
###'

set.seed(238923)  # for reproducibility
N <- 1000
theta <- rchisq(N, df = 10)
X <- rpois(n = N, lambda = theta)
tau <- seq(1, 32)
result <- deconv(tau = tau, X = X, ignoreZero = FALSE)
print(result$stats)

result

str(result)


###'######################################################################
###'
###' Twin Towers Example 
###' 
###' See Brad Efron: Bayes, Oracle Bayes and Empirical Bayes
###' disjointTheta is provided by deconvolveR package
###' 
###' 

theta <- disjointTheta 
N <- length(disjointTheta)
z <- rnorm(n = N, mean = disjointTheta)
tau <- seq(from = -4, to = 5, by = 0.2)
result <- deconv(tau = tau, X = z, family = "Normal", pDegree = 6)
g <- result$stats[, "g"]
if (require("ggplot2")) {
  ggplot() +
    geom_histogram(mapping = aes(x = disjointTheta, y = ..count.. / sum(..count..) ),
                   color = "blue", fill = "red", bins = 40, alpha = 0.5) +
    geom_histogram(mapping = aes(x = z, y = ..count.. / sum(..count..) ),
                   color = "brown", bins = 40, alpha = 0.5) +
    geom_line(mapping = aes(x = tau, y = g), color = "black") +
    labs(x = paste(expression(theta), "and x"), y = paste(expression(g(theta)), " and f(x)"))
}



###'######################################################################
###'
###' bardWordCount
###'
###'

data("bardWordCount")

bardWordCount

