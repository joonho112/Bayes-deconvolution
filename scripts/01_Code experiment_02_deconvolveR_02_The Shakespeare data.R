
###'######################################################################
###'
###' Category: Code Experiment
###' 
###' Task: deconvolveR - Empirical Bayes methods for learning prior distribution
###'       from data
###'       
###'       (2) The Shakespeare data
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
###' Load dataset
###'
###'

data(bardWordCount)
str(bardWordCount)
length(bardWordCount)



###'######################################################################
###'
###' Deconvolve the data to get G
###'
###'

### Take the support set \tau
lambda <- seq(-4, 4.5, .025)
tau <- exp(lambda)
length(tau)


### Deconvolve the data to get G
result <- deconv(tau = tau, y = bardWordCount, n = 100, c0 = 2)
stats <- result$stats


### Plot the EB deconvolution estimates for the Shakespeare word counts
ggplot() +
  geom_line(mapping = aes(x = lambda, y = stats[, "g"])) +
  labs(x = expression(log(theta)), y = expression(g(theta)))


### Extract the quantity R(\alpha) in the paper (Efron, Biometrika 2015) 
print(result$S)


###'######################################################################
###'
###' g vs. tg, when there is truncation at zero
###'
###'

### Plot G vs. TG
d <- data.frame(lambda = lambda, g = stats[, "g"], tg = stats[, "tg"], SE.g = stats[, "SE.g"])

indices <- seq(1, length(lambda), 5)

ggplot(data = d) +
  geom_line(mapping = aes(x = lambda, y = g)) +
  geom_errorbar(data = d[indices, ],
                mapping = aes(x = lambda, ymin = g - SE.g, ymax = g + SE.g),
                width = .01, color = "blue") +
  labs(x = expression(log(theta)), y = expression(g(theta))) +
  ylim(0, 0.006) +
  geom_line(mapping = aes(x = lambda, y = tg), linetype = "dashed", color = "red")



### Plot the posterior estimates
gPost <- sapply(seq_len(100), function(i) local({tg <- d$tg * result$P[i, ]; tg / sum(tg)}))

plots <- lapply(c(1, 2, 4, 8), function(i) {
  ggplot() +
    geom_line(mapping = aes(x = tau, y = gPost[, i])) +
    labs(x = expression(theta), y = expression(g(theta)),
         title = sprintf("x = %d", i))
})

plots <- Map(f = function(p, xlim) p + xlim(0, xlim), plots, list(6, 8, 14, 20))
plot_grid(plotlist = plots, ncol = 2)



###'######################################################################
###'
###' Bootstrap comparison
###'
###'

set.seed(1783)

B <- 200

fHat <- as.numeric(result$P %*% d$g)
fHat <- fHat / sum(fHat)
yStar <- rmultinom(n = B, size = sum(bardWordCount), prob = fHat)

gBoot <- apply(yStar, 2,
               function(y) deconv(tau = tau, y = y, n = 100, c0 = 2)$stats[, "g"])

seG <- apply(gBoot, 1, sd)

ggplot(data = d) +
  geom_line(mapping = aes(x = lambda, y = SE.g,
                          color = "Theoretical", linetype = "Theoretical")) +
  geom_line(mapping = aes(x = lambda, y = seG,
                          color = "Bootstrap", linetype = "Bootstrap")) +
  scale_color_manual(name = "Legend",
                     values = c("Bootstrap" = "black", "Theoretical" = "red")) +
  scale_linetype_manual(name = "Legend",
                        values = c("Bootstrap" = "solid", "Theoretical" = "dashed")) +
  theme(legend.title = element_blank()) +
  labs(x = expression(log(theta)), y = expression(sigma(hat(g))))



###'######################################################################
###'
###' Predict ratio of new distinct words
###'
###'

gHat <- stats[, "g"]

Rfn <- function(t) {
  sum( gHat * (1 - exp(-tau * t)) / (exp(tau) - 1) )
}

r <- sapply(0:10, Rfn)

ggplot() +
  geom_line(mapping = aes(x = 0:10, y = r)) +
  labs(x = "time multiple t", y = expression(R(t)))

print(uniroot(f = function(x) Rfn(x) - 1, interval = c(2, 4))$root)







