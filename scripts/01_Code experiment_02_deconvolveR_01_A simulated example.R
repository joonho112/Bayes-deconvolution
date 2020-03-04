
###'######################################################################
###'
###' Category: Code Experiment
###' 
###' Task: deconvolveR - Empirical Bayes methods for learning prior distribution
###'       from data
###'       
###'       (1) A Simulated Example
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
###' A simulation example
###'
###'

### Generate theta (true site-specific effects)
set.seed(238923) ## for reproducibility
N <- 1000
Theta <- rchisq(N, df = 10)


###' Generate observed data from theta (observed site-specific effects)
nSIM <- 1000
data <- sapply(seq_len(nSIM), function(x) rpois(n = N, lambda = Theta))
dim(data)


### Generate a vector of (implicitly m) discrite support points for theta
tau <- seq(1, 32)


### Apply the deconv function to estimate G
library(deconvolveR)
results <- apply(data, 2, function(x) deconv(tau = tau, X = x, ignoreZero = FALSE, c0 = 1))
str(results, 1)


### Construct a table of values for various values of theta
g <- sapply(results, function(x) x$stats[, "g"])
mean <- apply(g, 1, mean)

SE.g <- sapply(results, function(x) x$stats[, "SE.g"])
sd <- apply(SE.g, 1, mean)

Bias.g <- sapply(results, function(x) x$stats[, "Bias.g"])
bias <- apply(Bias.g, 1, mean)

gTheta <- pchisq(tau, df = 10) - pchisq(c(0, tau[-length(tau)]), df = 10)
gTheta <- gTheta / sum(gTheta)

simData <- data.frame(theta = tau, gTheta = gTheta,
                      Mean = mean, StdDev = sd, Bias = bias,
                      CoefVar = sd / mean)

table1 <- transform(simData,
                    gTheta = 100 * gTheta,
                    Mean = 100 * Mean,
                    StdDev = 100 * StdDev,
                    Bias = 100 * Bias)

knitr::kable(table1, row.names=FALSE)


### Plot
theme_set(theme_get() +
            theme(panel.grid.major = element_line(colour = "gray90", size = 0.2),
                  panel.grid.minor = element_line(colour = "gray98", size = 0.5)))


p1 <- ggplot(data = as.data.frame(results[[1]]$stats)) +
  geom_line(mapping = aes(x = theta, y = SE.g), 
  color = "black", linetype = "solid") +
  geom_line(mapping = aes(x = simData$theta, y = simData$StdDev), 
  color = "red", linetype = "dashed") +
  labs(x = expression(theta), y = "Std. Dev")


p2 <- ggplot(data = as.data.frame(results[[1]]$stats)) +
  geom_line(mapping = aes(x = theta, y = Bias.g), 
            color = "black", linetype = "solid") +
  geom_line(mapping = aes(x = simData$theta, y = simData$Bias), 
            color = "red", linetype = "dashed") +
  labs(x = expression(theta), y = "Bias")


plot_grid(plotlist = list(p1, p2), ncol = 2)

