
###'######################################################################
###'
###' Category: Code Experiment
###' 
###' Task: deconvolveR - Empirical Bayes methods for learning prior distribution
###'       from data
###'       
###'       (3) Normal example
###'       
###' Data: Simulated data
###' 
###' Data: 2020-03-03
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
###' Generate a dataset which follows big chunk of zero
###'
###'

set.seed(129023)
N <- 10000
pi0 <- .90

data <- local({
  nullCase <- (runif(N) <= pi0)
  muAndZ <- t(sapply(nullCase, function(isNull) {
    if (isNull) {
      mu <- 0
      c(mu, rnorm(1))
    } else {
      mu <- rnorm(1, mean = -3)
      c(mu, rnorm(1, mean = mu))
    }
  }))
  data.frame(nullCase = nullCase, mu = muAndZ[, 1], z = muAndZ[, 2])
})



###'######################################################################
###'
###' Plot true and observed data
###'
###'

p1 <- ggplot(mapping = aes(x = data$z)) +
  geom_histogram(mapping = aes(y  = ..count.. / sum(..count..) ),
                 color = "brown", bins = 60, alpha = 0.5) +
  labs(x = "z", y = "Density")

p2 <- ggplot(mapping = aes(x = data$mu)) +
  geom_histogram(mapping = aes(y  = ..count.. / sum(..count..) ),
                 color = "brown", bins = 60, alpha = 0.5) +
  labs(x = expression(theta), y = "Density")

plot_grid(plotlist = list(p1, p2), ncol = 2)



###'######################################################################
###'
###' Deconvolve the data z
###'
###'

summary(data$z)
summary(data$mu)
summary(data$nullCase)

tau <- seq(from = -6, to = 3, by = 0.25)

atomIndex <- which(tau == 0)

result <- deconv(tau = tau, X = data$z, deltaAt = 0, family = "Normal", pDegree = 5)


knitr::kable(result$stats)



gData <- as.data.frame(result$stats[-atomIndex, c("theta", "g")])

gData$g <- gData$g / sum(gData$g)

ggplot(data = gData) +
  geom_line(mapping = aes(x = theta, y = g)) +
  geom_line(mapping = aes(x = theta, y = dnorm(theta, mean = -3)),
            color = "red") +
  labs(x = expression(theta), y = expression(g(theta)))








