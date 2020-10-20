
###'######################################################################
###'
###' Category: Simulation Results
###' 
###' Task: Generate exploratory plots from the simulation results  
###' 
###'       (4) Fit meta-models for simulation results  
###'       
###'       4-1. Develop a function
###'       
###'       
###' Data: Results from the simulation
###' 
###' Data: 2020-04-24
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
data_dir <- file.path(work_dir, "datasets", 
                      "11_All collected and replicated loss estimates")


### Call libraries
library(tidyverse)
library(cowplot)
library(ggpubr)
library(ggthemes)
library(jtools)
library(sjPlot)
library(multcomp)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "functions/03_Simulation Implementation_04_Plot helpers.R")
source(file = "functions/03_Simulation Implementation_05_Data management helpers.R")



###'######################################################################
###'
###' Load the processed collection of estimates
###'
###'

### Load the dataset
setwd(data_dir)

df <- readRDS(file = "All collected and replicated loss estimates.rds") 


### Convert the cluster_ID into factor
df <- df %>%
  mutate(cluster_ID = factor(cluster_ID)) %>%
  mutate(SSEL = SSEL*1000, 
         ISEL = ISEL*1000, 
         SELrank = SELrank*1000)

classmode(df, everything())



###'######################################################################
###'
###' Filter only conditional cases
###'
###'

df_sub <- df %>%
  filter(estimator != c("ML")) %>%
  filter(truth %in% c("ALD")) %>%
  filter(assumed != c("T")) %>%
  filter(!(N %in% c("N = 100"))) %>%
  rename(I = ICC, R = rsr, anl = assumed, est = estimator)

  
names(df_sub)



###'######################################################################
###'
###' Simplify the factor labels
###'
###'

### Define new labels
lev_N <- unique(df_sub$N) 
lab_N <- lev_N %>% str_remove(" = ")

lev_I <- unique(df_sub$I)
lab_I <- lev_I %>% str_remove(" = ") 

lev_R <- unique(df_sub$R)
lab_R <- lev_R %>% str_remove(" = ") 

lev_anl <- unique(df_sub$anl)
lab_anl <- lev_anl %>% str_replace("-", "_")


### Mutate factors
df_sub <- df_sub %>%
  mutate(N = factor(N, levels = lev_N, labels = lab_N), 
         I = factor(I, levels = lev_I, labels = lab_I),
         R = factor(R, levels = lev_R, labels = lab_R),
         anl = factor(anl, levels = lev_anl, labels = lab_anl))



###'######################################################################
###'
###' get_lincom_results()
###' 
###' - Define a function to return lincom results
###' 
###'

get_lincom_results <- function(fit, df_sub, 
                               var1, var2, 
                               ref1, ref2, 
                               vec_ref){
  
  ###' Generate a dataframe for regression terms
  ###' Tag reference categories
  ###' Generate hypotheses
  var1_levels <- unique(df_sub[var1])[[var1]]
  var2_levels <- unique(df_sub[var2])[[var2]]
  
  df_vars <- expand.grid(var1_levels, var2_levels) %>% 
    data.frame() %>% set_names(c("ovar1", "ovar2")) %>%
    mutate(evar1 = paste0(var1, ovar1), 
           evar2 = paste0(var2, ovar2)) %>%
    unite(int_terms, evar1, evar2, sep = ":", remove = FALSE) %>%
    dplyr::select(ovar1, ovar2, evar1, evar2, int_terms) %>%
    mutate(tag1 = (ovar1 == ref1) %>% as.numeric(), 
           tag2 = (ovar2 == ref2) %>% as.numeric(), 
           tag_mod = if_else(tag1 == 0 & tag2 == 0, 1, 0), 
           elem1 = if_else(tag1 == 1, "", as.character(evar1)), 
           elem2 = if_else(tag2 == 1, "", as.character(evar2)), 
           elem3 = if_else(tag_mod == 0, "", as.character(int_terms)), 
           hypo = case_when(
             tag_mod == 1 ~ paste0(elem1, " + ", elem3, " = 0"), 
             tag1 == 0 & tag2 == 1 ~ paste0(elem1, " = 0"),
             tag1 == 1 & tag2 == 0 ~ paste0(elem2, " = 0"), 
             tag1 == 1 & tag2 == 1 ~ NA_character_, 
             TRUE ~ NA_character_)) %>%
    drop_na()
  
  
  ### Get lincom test results
  obj_ht <- glht(fit, linfct = df_vars$hypo)
  
  df_glht <- tidy(summary(obj_ht)) %>%
    full_join_track(tidy(confint(obj_ht)), by = c("lhs", "rhs", "estimate")) 
  
  df_lincom <- cbind.data.frame(df_vars, df_glht) %>%
    dplyr::select(-lhs, -rhs)
  
  
  ### Plot!
  df_plot <- df_lincom %>%
    dplyr::select(ovar1, ovar2, hypo, estimate, conf.low, conf.high) %>%
    filter(ovar1 != ref1)
  
  p <- ggplot(data = df_plot, aes(x = ovar2, y = estimate, 
                                  group = ovar1, color = ovar1, shape = ovar1)) +
    geom_point(aes(y = estimate), size = 3, position = position_dodge(width = 0)) + 
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                  position = position_dodge(width = 0), width = 0.1) + 
    geom_line(aes(y = estimate)) + 
    geom_hline(yintercept = 0, color = "gray30", linetype = "dashed") + 
    # geom_text(aes(label = round(estimate, 0)), color = "black", 
    #           position = position_stack(vjust = 0.1, reverse = TRUE), size = 3) + 
    scale_y_continuous(labels = comma) + 
    scale_color_manual(values = color_palette[seq(unique(df_plot$ovar1))]) + 
    scale_shape_manual(values = shape_palette[seq(unique(df_plot$ovar1))]) + 
    theme_trend +
    labs(subtitle = vec_ref)
  
  ### Return a list object
  list(df_lincom, df_plot, p)
}



###'######################################################################
###'
###' Test a meta-model fitting
###'
###'

### Set reference category for each factor
ref_N <- "N25"
ref_I <- "I0.1"
ref_R <- "R1"
ref_anl <- "Gaussian"
ref_est <- "PM"

vec_ref <- paste(ref_N, ref_I, ref_R, ref_anl, ref_est, sep = ", ")


df_sub_relevel <- df_sub %>%
  mutate(N = relevel(N, ref = ref_N), 
         I = relevel(I, ref = ref_I), 
         R = relevel(R, ref = ref_R), 
         anl = relevel(anl, ref = ref_anl), 
         est = relevel(est, ref = ref_est))


### Fit the OLS model with clustered SE

fit1 <- lm(ISEL ~ N + I + R + anl + est + 
             N:I + N:R + N:anl + N:est + 
             I:R + I:anl + I:est + 
             R:anl + R:est + 
             anl:est, data = df_sub_relevel)


# fit1 <- lm(SSEL ~ N*I + N*R + N*anl + N*est + 
#              I*R + I*anl + I*est + 
#              R*anl + R*est + 
#              anl*est, data = df_sub)


### Generate a table summary with cluster SEs
df_fit1 <- get_robust_se(fit1, cluster = "cluster_ID")$coefs %>% 
  data.frame() %>% round(1)

df_fit1


### Get lincom results
lincom1 <- get_lincom_results(fit1, df_sub_relevel, 
                              "anl", "est", ref_anl, ref_est, 
                              vec_ref)

lincom2 <- get_lincom_results(fit1, df_sub_relevel, 
                              "N", "I", ref_N, ref_I, 
                              vec_ref)

lincom3 <- get_lincom_results(fit1, df_sub_relevel, 
                              "I", "R", ref_I, ref_R, 
                              vec_ref)

lincom4 <- get_lincom_results(fit1, df_sub_relevel, 
                              "I", "anl", ref_I, ref_anl, 
                              vec_ref)




