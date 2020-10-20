
###'######################################################################
###'
###' Category: Simulation Results
###' 
###' Task: Generate exploratory plots from the simulation results  
###' 
###'       (3) Visualize loss estimates
###'       
###'       3-4. Quantile estimates
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
library(ggpubr)
library(ggthemes)


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

df_pct_temp <- readRDS(file = "All percentile estimates of G.rds")



###'######################################################################
###'
###' [ON/OFF] Filter the loaded data for more succinct presentation
###' 
###'

df_pct <- df_pct_temp %>%
  filter(truth %in% c("Gaussian", "ALD", "Mixed")) %>%
  filter(assumed %in% c("Gaussian", "DP-diffuse", "DP-Inform")) %>%
  filter(!(N %in% c("N = 100")))

# df_pct <- df_pct_temp



###'######################################################################
###' 
###'  (1) Fix: True DGM, Nsites, and rsr 
###'  
###'  - row: ICC
###'
###'

### Prepare parameters
vec_truth <- unique(df_pct$truth) 
vec_N <- unique(df_pct$N)
vec_rsr <- unique(df_pct$rsr) 
tbl_loop <- expand.grid(vec_truth, vec_N, vec_rsr)


### Define labels
lab_elem0 <- c("Bias_in_percentile_estimates")
title_vec <- "Bias in the percentile estimates"
caption_vec <- c("PM: Posterior Mean, CB: Constrained Bayes, GR: Triple-goal")


### Loop over elements

for (i in seq(nrow(tbl_loop))){
  
  ### Extract parameters
  elem1 <- as.character(tbl_loop[i, 1])
  elem2 <- as.character(tbl_loop[i, 2])
  elem3 <- as.character(tbl_loop[i, 3])
  
  
  ### Prepare a dataframe to plot 
  df_plot <- df_pct %>%
    filter(estimator != "ML") %>%
    filter(truth == elem1) %>%
    filter(N == elem2) %>% 
    filter(rsr == elem3) %>%
    mutate(diff = diff*100) %>%
    mutate(tail_p = factor(tail_p, 
                           levels = levels(tail_p), 
                           labels = c("5", "10", "25", "50", "75", "90", "95")))

  df_plot <- df_plot %>%
    rename(yvar = diff, xvar = tail_p, groupvar = estimator, 
           facet_row_var = ICC, 
           facet_col_var = assumed)
  
  
  ### Generate a lolipop plot
  p <- ggplot(df_plot) + 
    aes(x = xvar, y = yvar, group = groupvar) + 
    geom_segment(aes(x = xvar, xend = xvar, 
                     y = 0, yend = yvar, 
                     color = groupvar)) + 
    geom_point(aes(x = xvar, y = yvar, color = groupvar), 
               size = 2) + 
    geom_text(aes(label = round(yvar, 0), vjust = ifelse(yvar >= 0, -1.0, 1.5)), 
              color = "black", size = 3)
  
  
  ### Add horizontal line layer
  p <- p + geom_hline(aes(yintercept = 0), color = "black", linetype = "solid")
  
    
  ### Faceting
  p <- p + facet_grid(facet_row_var ~ facet_col_var + groupvar, 
                      scales = "fixed")
  
  
  ### Themes, temporary labels, and manual colors
  p <- p + theme_trend + temp_labels + 
    scale_y_continuous(expand = c(0.3, 0.3)) + 
    scale_color_manual(values = rev(color_palette[seq(unique(df_plot$groupvar))])) + 
    scale_shape_manual(values = rev(shape_palette[seq(unique(df_plot$groupvar))]))
  
  
  ### Add labels and save as .pdf
  p <- p + labs(title = title_vec, 
                subtitle = paste0(elem1, ",  ", elem2, ", ", elem3), 
                x = "True percentile", y = "Bias in percentile estimates", 
                caption = caption_vec)
  
  
  figure_name <- paste(sprintf("%02d", i), 
                       lab_elem0, 
                       elem1, 
                       elem2,
                       elem3, 
                       sep = "_")
  
  setwd(file.path(work_dir, "figures"))
  ggsave(paste0(figure_name, ".pdf"), p, width = 13, height = 9)

}



###'######################################################################
###' 
###'  (2) Fix: True DGM, Nsites, and ICC
###'  
###'  - row: rsr
###'
###'

### Prepare parameters
vec_truth <- unique(df_pct$truth) 
vec_N <- unique(df_pct$N)
vec_ICC <- unique(df_pct$ICC) 
tbl_loop <- expand.grid(vec_truth, vec_N, vec_ICC)


### Define labels
lab_elem0 <- c("Bias_in_percentile_estimates")
title_vec <- "Bias in the percentile estimates"
caption_vec <- c("PM: Posterior Mean, CB: Constrained Bayes, GR: Triple-goal")


### Loop over elements

for (i in seq(nrow(tbl_loop))){
  
  ### Extract parameters
  elem1 <- as.character(tbl_loop[i, 1])
  elem2 <- as.character(tbl_loop[i, 2])
  elem3 <- as.character(tbl_loop[i, 3])
  
  
  ### Prepare a dataframe to plot 
  df_plot <- df_pct %>%
    filter(estimator != "ML") %>%
    filter(truth == elem1) %>%
    filter(N == elem2) %>% 
    filter(ICC == elem3) %>%
    mutate(diff = diff*100) %>%
    mutate(tail_p = factor(tail_p, 
                           levels = levels(tail_p), 
                           labels = c("5", "10", "25", "50", "75", "90", "95")))
  
  df_plot <- df_plot %>%
    rename(yvar = diff, xvar = tail_p, groupvar = estimator, 
           facet_row_var = rsr, 
           facet_col_var = assumed)
  
  
  ### Generate a lolipop plot
  p <- ggplot(df_plot) + 
    aes(x = xvar, y = yvar, group = groupvar) + 
    geom_segment(aes(x = xvar, xend = xvar, 
                     y = 0, yend = yvar, 
                     color = groupvar)) + 
    geom_point(aes(x = xvar, y = yvar, color = groupvar), 
               size = 2) + 
    geom_text(aes(label = round(yvar, 0), vjust = ifelse(yvar >= 0, -1.0, 1.5)), 
              color = "black", size = 3)
  
  
  ### Add horizontal line layer
  p <- p + geom_hline(aes(yintercept = 0), color = "black", linetype = "solid")
  
  
  ### Faceting
  p <- p + facet_grid(facet_row_var ~ facet_col_var + groupvar, 
                      scales = "fixed")
  
  ### Themes, temporary labels, and manual colors
  p <- p + theme_trend + temp_labels + 
    scale_y_continuous(expand = c(0.3, 0.3)) + 
    scale_color_manual(values = rev(color_palette[seq(unique(df_plot$groupvar))])) + 
    scale_shape_manual(values = rev(shape_palette[seq(unique(df_plot$groupvar))]))
  
  
  ### Add labels and save as .pdf
  p <- p + labs(title = title_vec, 
                subtitle = paste0(elem1, ",  ", elem2, ", ", elem3), 
                x = "True percentile", y = "Bias in percentile estimates", 
                caption = caption_vec)
  
  
  figure_name <- paste(sprintf("%02d", i), 
                       lab_elem0, 
                       elem1, 
                       elem2,
                       elem3, 
                       sep = "_")
  
  setwd(file.path(work_dir, "figures"))
  ggsave(paste0(figure_name, ".pdf"), p, width = 13, height = 9)
  
}



###'######################################################################
###' 
###'  (3) Fix: ICC, R, Nsites
###'  
###'   - row: TrueDGM
###'
###'

### Prepare parameters
vec_ICC <- unique(df_pct$ICC) 
vec_R <- unique(df_pct$rsr)
vec_N <- unique(df_pct$N) 
tbl_loop <- expand.grid(vec_ICC, vec_R, vec_N)


### Define labels
lab_elem0 <- c("Bias_in_percentile_estimates")
title_vec <- "Bias in the percentile estimates"
caption_vec <- c("PM: Posterior Mean, CB: Constrained Bayes, GR: Triple-goal")


### Loop over elements

for (i in seq(nrow(tbl_loop))){
  
  ### Extract parameters
  elem1 <- as.character(tbl_loop[i, 1])
  elem2 <- as.character(tbl_loop[i, 2])
  elem3 <- as.character(tbl_loop[i, 3])
  
  
  ### Prepare a dataframe to plot 
  df_plot <- df_pct %>%
    filter(estimator != "ML") %>%
    filter(ICC == elem1) %>%
    filter(rsr == elem2) %>% 
    filter(N == elem3) %>%
    mutate(diff = diff*100) %>%
    mutate(tail_p = factor(tail_p, 
                           levels = levels(tail_p), 
                           labels = c("5", "10", "25", "50", "75", "90", "95")))
  
  df_plot <- df_plot %>%
    rename(yvar = diff, xvar = tail_p, groupvar = estimator, 
           facet_row_var = truth, 
           facet_col_var = assumed)
  
  
  ### Generate a lolipop plot
  p <- ggplot(df_plot) + 
    aes(x = xvar, y = yvar, group = groupvar) + 
    geom_segment(aes(x = xvar, xend = xvar, 
                     y = 0, yend = yvar, 
                     color = groupvar)) + 
    geom_point(aes(x = xvar, y = yvar, color = groupvar), 
               size = 2) + 
    geom_text(aes(label = round(yvar, 0), vjust = ifelse(yvar >= 0, -1.0, 1.5)), 
              color = "black", size = 3)
  
  
  ### Add horizontal line layer
  p <- p + geom_hline(aes(yintercept = 0), color = "black", linetype = "solid")
  
  
  ### Faceting
  p <- p + facet_grid(facet_row_var ~ facet_col_var + groupvar, 
                      scales = "fixed")
  
  
  ### Themes, temporary labels, and manual colors
  p <- p + theme_trend + temp_labels + 
    scale_y_continuous(expand = c(0.3, 0.3)) + 
    scale_color_manual(values = rev(color_palette[seq(unique(df_plot$groupvar))])) + 
    scale_shape_manual(values = rev(shape_palette[seq(unique(df_plot$groupvar))]))
  
  
  ### Add labels and save as .pdf
  p <- p + labs(title = title_vec, 
                subtitle = paste0(elem1, ",  ", elem2, ", ", elem3), 
                x = "True percentile", y = "Bias in percentile estimates", 
                caption = caption_vec)
  
  
  figure_name <- paste(sprintf("%02d", i), 
                       lab_elem0, 
                       elem1, 
                       elem2,
                       elem3, 
                       sep = "_")
  
  setwd(file.path(work_dir, "figures"))
  ggsave(paste0(figure_name, ".pdf"), p, width = 16, height = 9)
  
}



###'######################################################################
###' 
###'  (4) Fix: ICC, R, TrueDGM
###'  
###'   - row: Nsites
###'
###'

### Prepare parameters
vec_ICC <- unique(df_pct$ICC) 
vec_R <- unique(df_pct$rsr)
vec_truth <- unique(df_pct$truth) 
tbl_loop <- expand.grid(vec_ICC, vec_R, vec_truth)


### Define labels
lab_elem0 <- c("Bias_in_percentile_estimates")
title_vec <- "Bias in the percentile estimates"
caption_vec <- c("PM: Posterior Mean, CB: Constrained Bayes, GR: Triple-goal")


### Loop over elements

for (i in seq(nrow(tbl_loop))){
  
  ### Extract parameters
  elem1 <- as.character(tbl_loop[i, 1])
  elem2 <- as.character(tbl_loop[i, 2])
  elem3 <- as.character(tbl_loop[i, 3])
  
  
  ### Prepare a dataframe to plot 
  df_plot <- df_pct %>%
    filter(estimator != "ML") %>%
    filter(ICC == elem1) %>%
    filter(rsr == elem2) %>% 
    filter(truth == elem3) %>%
    mutate(diff = diff*100) %>%
    mutate(tail_p = factor(tail_p, 
                           levels = levels(tail_p), 
                           labels = c("5", "10", "25", "50", "75", "90", "95")))
  
  df_plot <- df_plot %>%
    rename(yvar = diff, xvar = tail_p, groupvar = estimator, 
           facet_row_var = N, 
           facet_col_var = assumed)
  
  
  ### Generate a lolipop plot
  p <- ggplot(df_plot) + 
    aes(x = xvar, y = yvar, group = groupvar) + 
    geom_segment(aes(x = xvar, xend = xvar, 
                     y = 0, yend = yvar, 
                     color = groupvar)) + 
    geom_point(aes(x = xvar, y = yvar, color = groupvar), 
               size = 2) + 
    geom_text(aes(label = round(yvar, 0), vjust = ifelse(yvar >= 0, -1.0, 1.5)), 
              color = "black", size = 3)
  
  
  ### Add horizontal line layer
  p <- p + geom_hline(aes(yintercept = 0), color = "black", linetype = "solid")
  
  
  ### Faceting
  p <- p + facet_grid(facet_row_var ~ facet_col_var + groupvar, 
                      scales = "fixed")
  
  
  ### Themes, temporary labels, and manual colors
  p <- p + theme_trend + temp_labels + 
    scale_y_continuous(expand = c(0.3, 0.3)) + 
    scale_color_manual(values = rev(color_palette[seq(unique(df_plot$groupvar))])) + 
    scale_shape_manual(values = rev(shape_palette[seq(unique(df_plot$groupvar))]))
  
  
  ### Add labels and save as .pdf
  p <- p + labs(title = title_vec, 
                subtitle = paste0(elem1, ",  ", elem2, ", ", elem3), 
                x = "True percentile", y = "Bias in percentile estimates", 
                caption = caption_vec)
  
  
  figure_name <- paste(sprintf("%02d", i), 
                       lab_elem0, 
                       elem1, 
                       elem2,
                       elem3, 
                       sep = "_")
  
  setwd(file.path(work_dir, "figures"))
  ggsave(paste0(figure_name, ".pdf"), p, width = 16, height = 9)
  
}


