
###'######################################################################
###'
###' Category: Publication Tables and Figures 
###' 
###' Task: Comparing effects of (1) G prior choice and (2) posterior summary methods
###' 
###'       [Figure 05]. Benefits of employing each method
###'           
###' Data: Simulated data
###' 
###' Data: 2020-04-23
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


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "functions/03_Simulation Implementation_04_Plot helpers.R")
source(file = "functions/03_Simulation Implementation_05_Data management helpers.R")



###'######################################################################
###'
###' Generate variables: Proportional reduction compared to ML estimates
###'
###' => Encapsulates the "Benefit" of using the method 
###'
###'

### Load the processed collection of estimates
setwd("~/Bayes-deconvolution/datasets/10_All collections_with DP-inform")
df <- readRDS(file = "All collection of loss estimates.rds")


### Subset only the selected cases
df_sub  <- df %>%
  filter(quantity %in% c("SSEL", "ISEL", "SELrank")) %>%
  filter(truth %in% c("Gaussian", "ALD", "Mixed")) %>%
  filter(assumed %in% c("Gaussian", "DP-diffuse", "DP-Inform"))


### Calculate proportional reductions
df_ML <- df_sub %>%
  filter(estimator == "ML") %>%
  select(-assumed, -estimator) %>%
  distinct() %>%
  rename(ML = value) 

df_bnf <- df_sub %>%
  filter(estimator != "ML") %>%
  full_join_track(df_ML, by = c("quantity", "truth", "N", "ICC", "rsr")) %>%
  mutate(benefit = ((ML - value)/ML)*100)
  


###'######################################################################
###'
###' Plot benefits using line graph  
###'
###'

### Prepare a dataframe to plot
df_plot <- df_bnf %>%
  filter(quantity == "SSEL") %>%
  filter(N == "N = 50") %>%
  filter(rsr == "R = 1")


### Plot
p <- four_dim_plot(df_plot, 
                   y = benefit, 
                   x = estimator,
                   group = assumed, 
                   facet_row = ICC, 
                   facet_col = truth, 
                   set_expand = c(0.3, 0.3), 
                   sprintf = "%.0f")














