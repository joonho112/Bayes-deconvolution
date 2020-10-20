
###'######################################################################
###'
###' Category: Simulation Results
###' 
###' Task: Generate exploratory plots from the simulation results  
###' 
###'       (3) Visualize loss estimates
###'       
###'       3-5. Histogram bars
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
library(hhsim)


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

df_hist_temp <- readRDS(file = "All histogram bars.rds")



###'######################################################################
###'
###' [ON/OFF] Filter the loaded data for more succinct presentation
###' 
###'

df_hist <- df_hist_temp %>%
  filter(truth %in% c("Gaussian", "ALD", "Mixed")) %>%
  filter(assumed %in% c("Gaussian", "DP-diffuse", "DP-Inform")) %>%
  filter(!(N %in% c("N = 100")))

# df_hist <- df_hist_temp



###'######################################################################
###'
###' Gen_True_his()
###'
###' - Define a function to generate the true EDF
###'   conforming to the true distribution
###'
###'

Gen_True_his <- function(gen_name,
                         N = 2000000, 
                         N_bins = 50, 
                         alow = -6, ahigh = 6){
  
  if (gen_name == "Gaussian"){ 
    
    mn <- 0; tau <- sqrt(1) 
    true_theta <- rnorm(N, mn, tau) 
    
  } else if (gen_name == "T"){   
    
    mn <- 0; nu <- 5
    true_theta <- rt(N, nu)*sqrt((nu - 2)/nu)
    
  } else if (gen_name == "ALD"){
    
    mean <- 0; var <- 1; p <- 0.1   
    scale <- sqrt((2*p^2*var)/(1 + p^4))
    location <- mean - ((scale*(1/p - p))/sqrt(2))
    true_theta <- LaplacesDemon::ralaplace(N, location, scale, p)
    
  } else if (gen_name == "Bimodal"){
    
    delta <- 4; eps <- 0.3; ups <- 1
    a <- sqrt((1 - eps) + eps*ups^2 + eps*(1 - eps)*delta^2)
    ind <- runif(N) < (1 - eps)
    true_theta <- ind*rnorm(N, -eps*delta/a, sqrt(1/a^2)) + 
      (1 - ind)*rnorm(N, (1 - eps)*delta/a, sqrt(ups^2/a^2))
    
  } else if (gen_name == "Skew"){
    
    mean <- 0; var <- 1; slant <- 10
    delta <- slant/sqrt(1 + slant^2)
    scale <- sqrt(var/(1-(2*delta^2/pi)))
    location <- mean - scale*sqrt(2/pi)*delta
    true_theta <- sn::rsn(n = N, xi = location, omega = scale, 
                          alpha = slant)[1:N]

  } else if (gen_name == "Mixed"){
    
    delta <- 5; eps <- 0.3; ups <- 2
    a <- sqrt((1 - eps) + eps*ups^2 + eps*(1 - eps)*delta^2)
    ind <- runif(N) < (1 - eps)
    true_theta <- ind*rnorm(N, -eps*delta/a, sqrt(1/a^2)) + 
      (1 - ind)*rnorm(N, (1 - eps)*delta/a, sqrt(ups^2/a^2))
    
  }
  
  ### Calculate EDF estimates for each bin 
  vec_his <- rep(0, N_bins)
  est_his <- enshis(N_bins, N, alow, ahigh, vec_his, true_theta)$his
  true_his <- est_his/N 
  
  ### Returen the resulting dataframe 
  true_his %>% as.data.frame() %>%
    rename(Bin = Var1, true_his = Freq) %>%
    mutate(truth = gen_name) %>%
    select(truth, Bin, true_his)
}


###'######################################################################
###'
###' Extract Bin information
###'
###'

### Split Bin into start and end points
df_splitted <- df_hist$Bin %>%
  str_split_fixed(",", n = 2) %>%
  data.frame() %>%
  mutate(start = str_sub(X1, start = 2) %>% as.numeric(), 
         end = str_sub(X2, end = -2) %>% as.numeric(), 
         middle = (start + end)/2) %>%
  dplyr::select(-X1, -X2)


### Combine with the original dataframe
df_temp <- cbind.data.frame(df_hist, df_splitted) 



###'######################################################################
###' 
###'  (1) Fix: True DGM, Nsites, and rsr 
###'  
###'  - row: ICC
###'
###'

### Prepare parameters
vec_truth <- unique(df_temp$truth) 
vec_N <- unique(df_temp$N)
vec_rsr <- unique(df_temp$rsr) 
tbl_loop <- expand.grid(vec_truth, vec_N, vec_rsr)


### Define labels
lab_elem0 <- c("Scaled_EDF_estimates")
title_vec <- "Scaled EDF estimates vs. True effect distribution"
caption_vec <- c("ML: Raw, PM: Posterior Mean, CB: Constrained Bayes, GR: Triple-goal")


### Loop over elements

for (i in seq(nrow(tbl_loop))){
  
  ### Extract parameters
  elem1 <- as.character(tbl_loop[i, 1])
  elem2 <- as.character(tbl_loop[i, 2])
  elem3 <- as.character(tbl_loop[i, 3])
  
  
  ### Prepare a dataframe to plot
  df_plot <- df_temp %>%
    filter(truth == elem1) %>%
    filter(N == elem2) %>%
    filter(rsr == elem3) 
  
  
  ###' Generate the TRUE effect distributions
  df_true_his <- Gen_True_his(gen_name = elem1)
  df_join <- full_join_track(df_plot, df_true_his, by = c("truth", "Bin"))
  
  
  ### Plot 
  p <- ggplot(df_join) + 
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = value), 
              fill = "limegreen",  color = "gray80", size = 0.0001) + 
    geom_line(aes(x = middle, y = true_his), size = 0.6) + 
    facet_grid(ICC ~ assumed + estimator) + 
    theme_trend + temp_labels
  
  
  ### Add labels and save as .pdf
  p <- p + labs(title = title_vec, 
                subtitle = paste0(elem1, ",  ", elem2, ", ", elem3), 
                x = NULL, y = "density", 
                caption = caption_vec)
  
  
  figure_name <- paste(sprintf("%02d", i), 
                       lab_elem0, 
                       elem1, 
                       elem2,
                       elem3, 
                       sep = "_")
  
  setwd(file.path(work_dir, "figures"))
  ggsave(paste0(figure_name, ".pdf"), p, width = 14, height = 6)

}



###'######################################################################
###' 
###'  (2) Fix: True DGM, Nsites, and ICC
###'  
###'  - row: rsr
###'
###'

### Prepare parameters
vec_truth <- unique(df_temp$truth) 
vec_N <- unique(df_temp$N)
vec_ICC <- unique(df_temp$ICC) 
tbl_loop <- expand.grid(vec_truth, vec_N, vec_ICC)


### Define labels
lab_elem0 <- c("Scaled_EDF_estimates")
title_vec <- "Scaled EDF estimates vs. True effect distribution"
caption_vec <- c("ML: Raw, PM: Posterior Mean, CB: Constrained Bayes, GR: Triple-goal")


### Loop over elements

for (i in seq(nrow(tbl_loop))){
  
  ### Extract parameters
  elem1 <- as.character(tbl_loop[i, 1])
  elem2 <- as.character(tbl_loop[i, 2])
  elem3 <- as.character(tbl_loop[i, 3])
  
  
  ### Prepare a dataframe to plot
  df_plot <- df_temp %>%
    filter(truth == elem1) %>%
    filter(N == elem2) %>%
    filter(ICC == elem3) 
  
  
  ###' Generate the TRUE effect distributions
  df_true_his <- Gen_True_his(gen_name = elem1)
  df_join <- full_join_track(df_plot, df_true_his, by = c("truth", "Bin"))
  
  
  ### Plot 
  p <- ggplot(df_join) + 
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = value), 
              fill = "limegreen",  color = "gray80", size = 0.0001) + 
    geom_line(aes(x = middle, y = true_his), size = 0.6) + 
    facet_grid(rsr ~ assumed + estimator) + 
    theme_trend + temp_labels
  
  
  ### Add labels and save as .pdf
  p <- p + labs(title = title_vec, 
                subtitle = paste0(elem1, ",  ", elem2, ", ", elem3), 
                x = NULL, y = "density", 
                caption = caption_vec)
  
  
  figure_name <- paste(sprintf("%02d", i), 
                       lab_elem0, 
                       elem1, 
                       elem2,
                       elem3, 
                       sep = "_")
  
  setwd(file.path(work_dir, "figures"))
  ggsave(paste0(figure_name, ".pdf"), p, width = 14, height = 6)
  
}



###'######################################################################
###' 
###'  (3) Fix: ICC, R, Nsites
###'  
###'   - row: TrueDGM
###'
###'

### Prepare parameters
vec_ICC <- unique(df_temp$ICC) 
vec_R <- unique(df_temp$rsr)
vec_N <- unique(df_temp$N) 
tbl_loop <- expand.grid(vec_ICC, vec_R, vec_N)


### Define labels
lab_elem0 <- c("Scaled_EDF_estimates")
title_vec <- "Scaled EDF estimates vs. True effect distribution"
caption_vec <- c("ML: Raw, PM: Posterior Mean, CB: Constrained Bayes, GR: Triple-goal")


### Loop over elements

for (i in seq(nrow(tbl_loop))){
  
  ### Extract parameters
  elem1 <- as.character(tbl_loop[i, 1])
  elem2 <- as.character(tbl_loop[i, 2])
  elem3 <- as.character(tbl_loop[i, 3])
  
  
  ### Prepare a dataframe to plot
  df_plot <- df_temp %>%
    filter(ICC == elem1) %>%
    filter(rsr == elem2) %>%
    filter(N == elem3) 
  
  
  ###' Generate the TRUE effect distributions
  ###' per each TrueDGM
  vec_truth <- unique(df_plot$truth)
  list_temp <- list()
  
  for (j in seq_along(vec_truth)){
    
    list_temp[[j]] <- Gen_True_his(gen_name = vec_truth[j])
  
  }
  
  df_true_his <- bind_rows(list_temp)
  
  df_join <- full_join_track(df_plot, df_true_his, by = c("truth", "Bin")) %>%
    mutate(truth = factor(truth, levels = vec_truth))
  
  
  ### Plot 
  p <- ggplot(df_join) + 
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = value), 
              fill = "limegreen",  color = "gray80", size = 0.0001) + 
    geom_line(aes(x = middle, y = true_his), size = 0.6) + 
    facet_grid(truth ~ assumed + estimator) + 
    theme_trend + temp_labels
  
  
  ### Add labels and save as .pdf
  p <- p + labs(title = title_vec, 
                subtitle = paste0(elem1, ",  ", elem2, ", ", elem3), 
                x = NULL, y = "density", 
                caption = caption_vec)
  
  
  figure_name <- paste(sprintf("%02d", i), 
                       lab_elem0, 
                       elem1, 
                       elem2,
                       elem3, 
                       sep = "_")
  
  setwd(file.path(work_dir, "figures"))
  ggsave(paste0(figure_name, ".pdf"), p, width = 14, height = 6)
  
}



###'######################################################################
###' 
###'  (4) Fix: ICC, R, TrueDGM
###'  
###'   - row: Nsites
###'
###'

### Prepare parameters
vec_ICC <- unique(df_temp$ICC) 
vec_R <- unique(df_temp$rsr)
vec_truth <- unique(df_temp$truth) 
tbl_loop <- expand.grid(vec_ICC, vec_R, vec_truth)


### Define labels
lab_elem0 <- c("Scaled_EDF_estimates")
title_vec <- "Scaled EDF estimates vs. True effect distribution"
caption_vec <- c("ML: Raw, PM: Posterior Mean, CB: Constrained Bayes, GR: Triple-goal")


### Loop over elements

for (i in seq(nrow(tbl_loop))){
  
  ### Extract parameters
  elem1 <- as.character(tbl_loop[i, 1])
  elem2 <- as.character(tbl_loop[i, 2])
  elem3 <- as.character(tbl_loop[i, 3])
  
  
  ### Prepare a dataframe to plot
  df_plot <- df_temp %>%
    filter(ICC == elem1) %>%
    filter(rsr == elem2) %>%
    filter(truth == elem3) 
  
  
  ###' Generate the TRUE effect distributions
  df_true_his <- Gen_True_his(gen_name = elem3)
  df_join <- full_join_track(df_plot, df_true_his, by = c("truth", "Bin"))
  
  
  ### Plot 
  p <- ggplot(df_join) + 
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = value), 
              fill = "limegreen",  color = "gray80", size = 0.0001) + 
    geom_line(aes(x = middle, y = true_his), size = 0.6) + 
    facet_grid(N ~ assumed + estimator) + 
    theme_trend + temp_labels
  
  
  ### Add labels and save as .pdf
  p <- p + labs(title = title_vec, 
                subtitle = paste0(elem1, ",  ", elem2, ", ", elem3),
                x = NULL, y = "density", 
                caption = caption_vec)
  
  
  figure_name <- paste(sprintf("%02d", i), 
                       lab_elem0, 
                       elem1, 
                       elem2,
                       elem3, 
                       sep = "_")
  
  setwd(file.path(work_dir, "figures"))
  ggsave(paste0(figure_name, ".pdf"), p, width = 14, height = 6)
  
}

