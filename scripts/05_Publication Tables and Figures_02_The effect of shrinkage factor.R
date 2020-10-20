
###'######################################################################
###'
###' Category: Publication Tables and Figures 
###' 
###' Task: Present shrinkage factor effect 
###' 
###'  [Table 01]. The effect of shrinkage factor
###'  [Figure 04]. The scaled EDF by the level of shrinkage factor
###'           
###' Data: Simulated data
###' 
###' Data: 2020-04-22
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
###' [Table 01]. The effect of shrinkage factor
###'
###' - Varying factor: 
###' 
###'   (1) shrinkage factor I
###'   (2) True data-generating models
###'   (3) Data-analytic models
###' 
###' - Fixed factor: 
###'   
###'   (1) N = 50
###'   (2) R = 1
###' 
###'

### Load the processed collection of estimates
setwd("~/Bayes-deconvolution/datasets/10_All collections_with DP-inform")
df <- readRDS(file = "All collection of loss estimates.rds")


### Subset only the selected cases
df_sub  <- df %>%
  filter(quantity %in% c("SSEL", "ISEL", "SELrank")) %>%
  filter(truth %in% c("Gaussian", "ALD", "Mixed")) %>%
  filter(assumed %in% c("Gaussian", "DP-diffuse", "DP-Inform")) %>%
  filter(N == "N = 50") %>%
  filter(rsr == "R = 1") %>%
  mutate(value = round(1000*value, 0)) 
  

### Unite column variables 
df_united <- df_sub %>%
  unite(column, assumed, estimator) 

col_levels <- unique(df_united$column)
qt_levels <- c("SSEL", "ISEL", "SELrank")

df_united <- df_united %>%
  mutate(column = factor(column, levels = col_levels), 
         quantity = factor(quantity, levels = qt_levels))


### Convert to the wide format
df_wide <- df_united %>%  
  spread(key = column, value = value) %>%
  select(N, rsr, quantity, ICC, truth, everything()) %>%
  arrange(N, rsr, quantity, ICC, truth) %>%
  select(-"DP-diffuse_ML", -"DP-Inform_ML")


### Save as .csv file
setwd(file.path(work_dir, "figures", "17_Publication plots"))
write.csv(df_wide, file = "Table01_The effect of shrinkage factor.csv")



###'######################################################################
###'
###' [Figure 04]. The scaled EDF by the level of shrinkage factor
###'
###' - Varying factor: 
###' 
###'   (1) Shrinkage factor I
###'   (2) True data-generating models
###'   (3) Data-analytic models
###' 
###' - Fixed factor: 
###'   
###'   (1) N = 50
###'   (2) R = 1
###' 
###'

### Load the processed collection of estimates
setwd("~/Bayes-deconvolution/datasets/10_All collections_with DP-inform")
df <- readRDS(file = "All histogram bars.rds")


### Subset only the selected cases
df_hist  <- df %>%
  filter(assumed %in% c("Gaussian", "DP-diffuse", "DP-Inform")) %>%
  filter(N == "N = 50") %>%
  filter(rsr == "R = 1") %>%
  filter(truth == "Mixed") %>%
  filter(estimator != c("CB")) %>%
  filter(!(assumed %in% c("DP-diffuse", "DP-Inform") & estimator == "ML"))


### Split Bin into start and end points
df_splitted <- df_hist$Bin %>%
  str_split_fixed(",", n = 2) %>%
  data.frame() %>%
  mutate(start = str_sub(X1, start = 2) %>% as.numeric(), 
         end = str_sub(X2, end = -2) %>% as.numeric(), 
         middle = (start + end)/2) %>%
  dplyr::select(-X1, -X2)


### Combine with the original dataframe
df_plot <- cbind.data.frame(df_hist, df_splitted) 


### Define labels
title_vec <- "Scaled EDF estimates vs. True effect distribution (A mixture of two Gaussians)"
subtitle_vec <- "N = 50, R = 1"
caption_vec <- c("PM: Posterior Mean, GR: Triple-goal")


### Generate the TRUE effect distributions
df_true_his <- Gen_True_his(gen_name = "Mixed")
df_join <- full_join_track(df_plot, df_true_his, by = c("truth", "Bin"))


### Plot 
p <- ggplot(df_join) + 
  geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = value), 
            fill = "limegreen",  color = "gray80", size = 0.0001) + 
  geom_line(aes(x = middle, y = true_his), size = 0.4) + 
  facet_grid(ICC ~ assumed + estimator) + 
  theme_minimal() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        strip.text.x = element_text(size = 12), 
        strip.text.y = element_text(size = 12)) + 
  labs(title = title_vec, 
       subtitle = subtitle_vec, 
       x = NULL, y = "Proportions", 
       caption = caption_vec)


### Save the resulting plot
setwd(file.path(work_dir, "figures"))
ggsave("Figure04_Scaled EDF estimates_for a mixture.pdf",
       p, width = 14, height = 5)


