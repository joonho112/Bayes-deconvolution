
###'######################################################################
###'
###' Category: Simulation Results
###' 
###' Task: Generate TABLES and FIGURES from the simulation results  
###'       
###' Data: Results from the simulation
###' 
###' Data: 2020-04-01
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
gc()            # force R to release memory it is no longer using
rm(list=ls())   # delete all the objects in the workspace


### Set working directory 
work_dir <- c("~/Bayes-deconvolution")
setwd(work_dir)


### Set a data directory
# data_dir <- file.path(work_dir, "datasets", "04_Simulation results")
data_dir <- c("D:/Data/DP_Simulation")


### Call libraries
library(tidyverse)
library(cowplot)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "functions/03_Simulation Implementation_04_Plot helpers.R")
source(file = "functions/03_Simulation Implementation_05_Data management helpers.R")



###'######################################################################
###'
###' Load the table containing the simulation parameters
###'
###'

setwd(work_dir)
sim_parm <- read.csv(file = "tables/simulation parameters.csv")



###'######################################################################
###'
###' Conditionally select a subgroup of folders
###'
###'

###' Parameter selection
df_select <- sim_parm %>%
  filter(n == 50) 


### Return a vector for the selected subgroup of folders
folder_list <- folder_select(df_select, data_dir)



###'######################################################################
###'
###' Collect the selected .csv files within the selelcted folders
###'
###'

df_collect <- collect_selected_files(folder_dir = data_dir, 
                                     folder_list = folder_list, 
                                     file_name = "09_Collection of loss estimates.csv", 
                                     sim_parm = sim_parm)

### Extract true variances
df_truevar <- df_collect %>%
  filter(estimator == "true")

df_collect <- df_collect %>% 
  filter(estimator != "true")



###'######################################################################
###'
###' Define factors 
###'
###'

### Check classes
classmode(df_collect, everything())


### (1) truth
tabdf(df_collect, truth)

vec_level <- c("Gaussian", "T", "Mix", "Skew")
vec_label <- c("Gaussian", "T (df = 5)", "Gaussian Mixture", "Skewed Normal")

df_collect <- df_collect %>%
  mutate(truth = factor(truth, levels = vec_level, labels = vec_label))


### (2) assumed
tabdf(df_collect, assumed)

vec_level <- c("Gaussian", "T", "DP")
vec_label <- c("Gaussian", "T_5", "DP")

df_collect <- df_collect %>%
  mutate(assumed = factor(assumed, levels = vec_level, labels = vec_label))


### (3) n
tabdf(df_collect, n)

vec_level <- unique(df_collect$n)
vec_label <- paste0("N = ", unique(df_collect$n))

df_collect <- df_collect %>%
  mutate(n = factor(n, levels = vec_level, labels = vec_label))


### (4) rgt
tabdf(df_collect, rgt)

vec_level <- unique(df_collect$rgt)
vec_label <- paste0("rgt = ", unique(df_collect$rgt))

df_collect <- df_collect %>%
  mutate(rgt = factor(rgt, levels = vec_level, labels = vec_label))


### (5) rsr
tabdf(df_collect, rsr)

vec_level <- unique(df_collect$rsr)
vec_label <- paste0("rsr = ", unique(df_collect$rsr))

df_collect <- df_collect %>%
  mutate(rsr = factor(rsr, levels = vec_level, labels = vec_label))


### (6) estimator
tabdf(df_collect, estimator)

vec_level <- c("ML", "PM", "CB", "GR")
vec_label <- c("ML", "PM", "CB", "GR")

df_collect <- df_collect %>%
  mutate(estimator = factor(estimator, levels = vec_level, labels = vec_label)) %>%
  dplyr::select(-X) %>%
  arrange(truth, assumed, n, rsr, rgt, estimator)



###'######################################################################
###'
###' (1) Sum of Squared Error Loss (SSEL) for 
###'     the individual site-specific effects
###'     
###' - Loop over true data-generating distributions
###'
###'

### Assign necessary objects
vec_truth <- levels(df_collect$truth) 
quantity_val <- c("SSEL")
N_sites <- 50

title_vec <- "Sum of Squared Error Loss (SSEL) for the individual site-specific effects"
caption_vec <- c("PM: Posterior Mean, CB: Constrained Bayes, GR: Triple-goal")


### Implement a loop
for (i in seq_along(vec_truth)){
  
  ### Extract loop elements
  true_dist <- vec_truth[i]
  
  
  ### Define subtitle & figure name
  subtitle_vec <- paste0("True data-generating dist.: ", true_dist, 
                         ", Number of sites = ", N_sites)
  
  figure_name <- paste("01_SSEL for individual effects", 
                       sprintf("%02d", i), 
                       true_dist, sep = "_")
  
  
  ### Prepare a dataframe to plot
  df_plot <- df_collect %>%
    filter(truth == true_dist) %>%
    filter(quantity == quantity_val) %>%
    filter(estimator != "ML")
  
  
  ### Plot!
  p <- plot_lines_grp(df_plot, x = estimator, y = value, group = assumed, 
                      sprintf = "%.3f") + 
    facet_wrap(rsr ~ rgt, scales = "free") + 
    labs(title = title_vec, 
         subtitle = subtitle_vec,
         x = "Estimator", y = quantity_val, 
         caption = caption_vec)
  
  setwd(file.path(work_dir, "figures"))
  ggsave(paste0(figure_name, ".pdf"), p, width = 14, height = 10)
  
}



###'######################################################################
###'
###' (2) Integrated SEL (ISEL) for 
###'     the estimated "empirical distribution function (EDF)"
###'     
###' - Loop over true data-generating distributions
###'
###'

### Assign necessary objects
vec_truth <- levels(df_collect$truth) 
quantity_val <- c("ISEL")
N_sites <- 50

title_vec <- "Integrated Squared Error Loss (SSEL) for the empirical distribution function (EDF)"
caption_vec <- c("PM: Posterior Mean, CB: Constrained Bayes, GR: Triple-goal")


### Implement a loop
for (i in seq_along(vec_truth)){
  
  ### Extract loop elements
  true_dist <- vec_truth[i]
  
  
  ### Define subtitle & figure name
  subtitle_vec <- paste0("True data-generating dist.: ", true_dist, 
                         ", Number of sites = ", N_sites)
  
  figure_name <- paste("02_ISEL for EDF", 
                       sprintf("%02d", i), 
                       true_dist, sep = "_")
  
  
  ### Prepare a dataframe to plot
  df_plot <- df_collect %>%
    filter(truth == true_dist) %>%
    filter(quantity == quantity_val) %>%
    filter(estimator != "ML")
  
  
  ### Plot!
  p <- plot_lines_grp(df_plot, x = estimator, y = value, group = assumed, 
                      sprintf = "%.3f") + 
    facet_wrap(rsr ~ rgt, scales = "free") + 
    labs(title = title_vec, 
         subtitle = subtitle_vec,
         x = "Estimator", y = quantity_val, 
         caption = caption_vec)
  
  setwd(file.path(work_dir, "figures"))
  ggsave(paste0(figure_name, ".pdf"), p, width = 14, height = 10)
  
}



###'######################################################################
###'
###' (3) Squared Error Loss (SEL) for 
###'     the estimated ranks
###'     
###' - Loop over true data-generating distributions
###'
###'

### Assign necessary objects
vec_truth <- levels(df_collect$truth) 
quantity_val <- c("SELrank")
N_sites <- 50

title_vec <- "Sum of Squared Error Loss (SSEL) for the rank estimates"
caption_vec <- c("PM: Posterior Mean, CB: Constrained Bayes, GR: Triple-goal")


### Implement a loop
for (i in seq_along(vec_truth)){
  
  ### Extract loop elements
  true_dist <- vec_truth[i]
  
  
  ### Define subtitle & figure name
  subtitle_vec <- paste0("True data-generating dist.: ", true_dist, 
                         ", Number of sites = ", N_sites)
  
  figure_name <- paste("03_SSEL for rank", 
                       sprintf("%02d", i), 
                       true_dist, sep = "_")
  
  
  ### Prepare a dataframe to plot
  df_plot <- df_collect %>%
    filter(truth == true_dist) %>%
    filter(quantity == quantity_val) %>%
    filter(estimator != "ML")
  
  
  ### Plot!
  p <- plot_lines_grp(df_plot, x = estimator, y = value, group = assumed, 
                      sprintf = "%.3f") + 
    facet_wrap(rsr ~ rgt, scales = "free") + 
    labs(title = title_vec, 
         subtitle = subtitle_vec,
         x = "Estimator", y = quantity_val, 
         caption = caption_vec)
  
  setwd(file.path(work_dir, "figures"))
  ggsave(paste0(figure_name, ".pdf"), p, width = 14, height = 10)
  
}


