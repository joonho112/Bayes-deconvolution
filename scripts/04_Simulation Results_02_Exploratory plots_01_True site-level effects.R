
###'######################################################################
###'
###' Category: Simulation Results
###' 
###' Task: Generate TABLES and FIGURES from the simulation results  
###' 
###'       (1) True site-level effect distributions
###'       
###' Data: Results from the simulation
###' 
###' Data: 2020-04-09
###' 
###' Author: JoonHo Lee (`joonho@berkeley.edu`)
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


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "functions/03_Simulation Implementation_04_Plot helpers.R")
source(file = "functions/03_Simulation Implementation_05_Data management helpers.R")



###'######################################################################
###'
###' Generate theoretical PDF with very large N (N = 20000)
###'
###'

N <- 2000000
ICC <- 0.5
rsr <- 5


### (1) Gaussian
Gaussian <- Gen_Gaussian(N, ICC, rsr)$theta


### (2) T (df = 5)
T5 <- Gen_T(N, ICC, rsr, nu = 5)$theta


### (3) Asymmetric Laplace disribution (p = 0.1)
ALD <- Gen_ALD(N, ICC, rsr)$theta


### (4) Bimodal distribution
delta <- 4   # distance between two means
eps <- 0.3   # proportion of the small component
ups <- 1     # ratio between two variances

Bimodal <- Gen_Mixture(N, ICC, rsr, delta, eps, ups)$theta 


### (5) Skewed Normal
slant <- 10   # control the degree of skewness

Skew <- Gen_SkewN(N, ICC, rsr, mean = 0, var = 1, slant = slant)$theta


### (6) Mixed distribution
delta <- 5   # distance between two means
eps <- 0.3  # proportion of the small component
ups <- 2     # ratio between two variances

Mixed <- Gen_Mixture(N, ICC, rsr, delta, eps, ups)$theta 



###'######################################################################
###'
###' Generate a plot comparing true site-level distributions
###'
###'

### Collect samples
vec_levels <- c("Gaussian", "T5", "ALD", "Bimodal", "Skew", "Mixed")

df_smp <- data.frame(Gaussian, T5, ALD, Bimodal, Skew, Mixed) %>%
  gather(key = DGM, value = value) %>%
  mutate(DGM = factor(DGM, levels = vec_levels))


### Plot by grid
p <- ggplot(data = df_smp, aes(x = value)) + 
  geom_density() + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") + 
  facet_wrap(DGM ~ .) + 
  xlim(-6, 6) + 
  theme_bw() + 
  labs(title = "Models for true site-level effects", 
       subtitle = "with mean 0 and variance 1", 
       x = NULL, y = "Density", 
       caption = "Range truncated at [-6, 6]")

p

ggsave(filename = "figures/True site-level distributions.pdf", p, 
       width = 12, height = 6)



###'######################################################################
###'
###' ALD distribution with p = 0.1
###'
###'

p <- ggplot(data = data.frame(ALD), aes(x = ALD)) + 
  geom_density() + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") + 
  theme_bw() + 
  labs(title = "Models for true site-level effects", 
       subtitle = "Asymmetric Laplace distribution (p = 0.1) with mean 0, variance 1", 
       x = NULL, y = "Density")

p

ggsave(filename = "figures/True site-level distributions-ALD.pdf", p, 
       width = 10, height = 6)



###'######################################################################
###'
###' [Publication Figure] Gaussian, Mixed, and ALD
###'
###'

### Plot grid setting
theme_temp <- theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))


### (A) Gaussian
p1 <- ggplot(data = data.frame(Gaussian), aes(x = Gaussian)) + 
  geom_density() + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") + 
  ylim(0, 1) + 
  theme_temp + 
  labs(title = "Gaussian", 
       x = NULL, y = "Density")


### (B) Mixed
p2 <- ggplot(data = data.frame(Mixed), aes(x = Mixed)) + 
  geom_density() + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") + 
  ylim(0, 1) + 
  theme_temp + 
  labs(title = "Gaussian Mixture", 
       x = NULL, y = "Density")


### (c) ALD
p3 <- ggplot(data = data.frame(ALD), aes(x = ALD)) + 
  geom_density() + 
  geom_vline(xintercept = 0, color = "red", linetype = "dashed") + 
  theme_temp + 
  ylim(0, 1) + 
  scale_x_continuous(breaks = seq(from = -2.5, to = 12.5, by = 2.5)) + 
  labs(title = "Asymmetric Laplace", 
       x = NULL, y = "Density")


### Combine (A) and (B) 
p_12 <- plot_grid(p1, p2, labels = "AUTO", rel_widths = c(1, 1))


### Combine (A) + (B) with (C)
p_combine <- plot_grid(p_12, p3, labels = c("", "C"), ncol = 1)


### Save the resulting plot 
ggsave(filename = "figures/Figure01_True site-level distributions.pdf", 
       p_combine, 
       width = 8, height = 6)





