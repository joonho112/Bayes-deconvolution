
###'######################################################################
###'
###' Category: Simulation Results
###' 
###' Task: Generate exploratory plots from the simulation results  
###' 
###'       (3) Visualize loss estimates
###'       
###'       3-3. Mean squared error loss (MSEL) for 
###'            the ranks
###'       
###' Data: Results from the simulation
###' 
###' Data: 2020-04-11
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
data_dir <- file.path(work_dir, "datasets", "10_All collections_with DP-inform")


### Call libraries
library(tidyverse)
library(cowplot)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "functions/03_Simulation Implementation_04_Plot helpers.R")
source(file = "functions/03_Simulation Implementation_05_Data management helpers.R")



###'######################################################################
###'
###' Load the processed collection of estimates
###'
###'

setwd(data_dir)

df_loss_temp <- readRDS(file = "All collection of loss estimates.rds")



###'######################################################################
###'
###' [ON/OFF] Filter the loaded data for more succinct presentation
###' 
###'

df_loss <- df_loss_temp %>%
  filter(truth %in% c("Gaussian", "ALD", "Mixed")) %>%
  filter(assumed %in% c("Gaussian", "DP-diffuse", "DP-Inform")) %>%
  filter(!(N %in% c("N = 100")))

# df_loss <- df_loss_temp



###'######################################################################
###' 
###'  (1) Fix: True DGM and N
###'      
###'      By I and R 
###'
###'

### Prepare parameters
vec_truth <- unique(df_loss$truth) 
vec_N <- unique(df_loss$N)
tbl_loop <- expand.grid(vec_truth, vec_N)


### Define labels
lab_elem0 <- c("MSEL_for_ranks")
title_vec <- "Mean Squared Error Loss (MSEL) for the rank estimates"
caption_vec <- c("PM: Posterior Mean, CB: Constrained Bayes, GR: Triple-goal")


### Loop over elements

for (i in seq(nrow(tbl_loop))){
  
  ### Extract parameters
  elem0 <- c("SELrank")
  elem1 <- as.character(tbl_loop[i, 1])
  elem2 <- as.character(tbl_loop[i, 2])
  
  
  ### Prepare a dataframe to plot
  df_plot <- df_loss %>%
    filter(estimator != "ML") %>%
    filter(quantity == elem0) %>%
    filter(truth == elem1) %>%
    filter(N == elem2) %>%
    mutate(value = value*1000)
  
  
  ### Plot!
  p <- four_dim_plot(df_plot, 
                     y = value, 
                     x = estimator,
                     group = assumed, 
                     facet_row = ICC, 
                     facet_col = rsr, 
                     set_expand = c(0.3, 0.3), 
                     sprintf = "%.0f")
  
  
  ### Add labels and save as .pdf
  p <- p + labs(title = title_vec, 
                subtitle = paste0(elem1, ",  ", elem2),
                x = "Estimator", y = "MSEL", 
                caption = caption_vec)
  
  
  figure_name <- paste(sprintf("%02d", i), 
                       lab_elem0, 
                       elem1, 
                       elem2,  
                       sep = "_")
  
  setwd(file.path(work_dir, "figures"))
  ggsave(paste0(figure_name, ".pdf"), p, width = 8, height = 8)
  
}


###'######################################################################
###' 
###'  (2) Fix: ICC & R 
###'  
###'      By TrueDGM & N
###'
###'

### Prepare parameters
vec_ICC <- unique(df_loss$ICC) 
vec_R <- unique(df_loss$rsr)
tbl_loop <- expand.grid(vec_ICC, vec_R)


### Define labels
lab_elem0 <- c("MSEL_for_ranks")
title_vec <- "Mean Squared Error Loss (MSEL) for the rank estimates"
caption_vec <- c("PM: Posterior Mean, CB: Constrained Bayes, GR: Triple-goal")


### Loop over elements

for (i in seq(nrow(tbl_loop))){
  
  ### Extract parameters
  elem0 <- c("SELrank")
  elem1 <- as.character(tbl_loop[i, 1])
  elem2 <- as.character(tbl_loop[i, 2])
  
  
  ### Prepare a dataframe to plot
  df_plot <- df_loss %>%
    filter(estimator != "ML") %>%
    filter(quantity == elem0) %>%
    filter(ICC == elem1) %>%
    filter(rsr == elem2) %>%
    mutate(value = value*1000)
  
  
  ### Plot!
  p <- four_dim_plot(df_plot, 
                     y = value, 
                     x = estimator,
                     group = assumed, 
                     facet_row = truth, 
                     facet_col = N, 
                     set_scale = "fixed", 
                     set_expand = c(0.3, 0.3), 
                     sprintf = "%.0f")
  
  
  ### Add labels and save as .pdf
  p <- p + labs(title = title_vec, 
                subtitle = paste0(elem1, ",  ", elem2),
                x = "Estimator", y = "MSEL", 
                caption = caption_vec)
  
  
  figure_name <- paste(sprintf("%02d", i), 
                       lab_elem0, 
                       elem1, 
                       elem2,  
                       sep = "_")
  
  setwd(file.path(work_dir, "figures"))
  ggsave(paste0(figure_name, ".pdf"), p, width = 8, height = 8)
  
}

