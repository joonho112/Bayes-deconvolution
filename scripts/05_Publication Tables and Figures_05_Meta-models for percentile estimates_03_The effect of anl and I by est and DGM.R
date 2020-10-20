
###'######################################################################
###'
###' Category: Publication Tables and Figures 
###' 
###' Task: Interaction plots
###' 
###'  [Figure 09]. Interaction plots for Meta-model regression results 
###'  
###'   - Raw biases for percentile estimates
###'   - The effect of DP priors by DGM, est, I
###'     
###'     * group factor: anl (DP priors)
###'     * x-axis: I shrinkage factor
###'     * facet factor 1: estimators (PM, CB, GR)
###'     * facet factor 2: DGM (ALD, Mixture) 
###'           
###' Data: Simulated data
###' 
###' Data: 2020-04-30
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
                      "12_All collected and replicated quantile estimates")


### Call libraries
library(tidyverse)
library(cowplot)
library(multcomp)
library(jtools)
library(broom)


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
df <- readRDS(file = "All collected and replicated quantile est.rds") 


### Convert the cluster_ID into factor
df <- df %>%
  mutate(cluster_ID = factor(cluster_ID)) %>%
  mutate(tail_p = 100*tail_p, 
         quantile = 100*quantile, 
         diff = 100*diff)

classmode(df, everything())



###'######################################################################
###'
###' Simplify the factor labels
###'
###'

### Rename variables
df <- df %>%
  rename(I = ICC, R = rsr, anl = assumed, est = estimator)


### Define new labels
lev_N <- unique(df$N) 
lab_N <- lev_N %>% str_remove(" = ")

lev_I <- unique(df$I)
lab_I <- lev_I %>% str_remove(" = ") 

lev_R <- unique(df$R)
lab_R <- lev_R %>% str_remove(" = ") 

lev_anl <- unique(df$anl)
lab_anl <- lev_anl %>% str_replace("-", "_")


### Mutate factors
df <- df %>%
  mutate(N = factor(N, levels = lev_N, labels = lab_N), 
         I = factor(I, levels = lev_I, labels = lab_I),
         R = factor(R, levels = lev_R, labels = lab_R),
         anl = factor(anl, levels = lev_anl, labels = lab_anl))



###'######################################################################
###'
###' Subset only conditional cases
###'
###'

### Tabulate existing variables
names(df)
tabdf(df, truth)
tabdf(df, N)
tabdf(df, anl)
tabdf(df, I)
tabdf(df, est)


### Filter only conditional cases
df_cond <- df %>%
  filter(truth %in% c("Gaussian", "ALD", "Mixed")) %>%
  filter(est != c("ML")) %>%
  filter(anl != c("T")) %>%
  filter(!(N %in% c("N100"))) 



###'######################################################################
###'
###' [Figure 08]. 
###' 
###' [Panel A] Effects of DP priors on the bias of percentile estimates
###' 
###' By estimator, shrinkage factor I, and DGM
###' 
###'

### Prepare a loop
list_lincom <- list()
list_df_plot <- list()


### Generate a table for looping
vec_truth <- c("ALD", "Mixed")
vec_est <- c("PM", "CB", "GR")
vec_tail_p <- unique(df_cond$tail_p)
df_grid <- expand.grid(vec_tail_p, vec_est, vec_truth) %>%
  data.frame() %>%
  dplyr::select(Var3, Var2, Var1) %>%
  set_names(c("truth", "est", "tail_p")) 


### Start loop 
for (i in seq(nrow(df_grid))){
  
  ### Extract loop elements
  elem_truth <- df_grid$truth[i]
  elem_est <- df_grid$est[i] %>% as.character()
  elem_tail_p <- df_grid$tail_p[i]
  
  
  ### Filter only conditional cases
  df_sub <- df_cond %>%
    filter(truth %in% elem_truth) %>% 
    filter(tail_p %in% elem_tail_p) 
  
  
  ### Set reference category for each factor
  ref_N <- "N25"
  ref_I <- "I0.1"
  ref_R <- "R1"
  ref_anl <- "Gaussian"
  ref_est <- elem_est
  
  vec_ref <- paste(elem_truth, ref_N, ref_I, ref_R, 
                   ref_anl, ref_est, sep = ", ")
  
  df_sub_relevel <- df_sub %>%
    mutate(N = relevel(N, ref = ref_N), 
           I = relevel(I, ref = ref_I), 
           R = relevel(R, ref = ref_R), 
           anl = relevel(anl, ref = ref_anl), 
           est = relevel(est, ref = ref_est))
  
  
  ###' (1) SSEL for individual site-specific effects
  ###' "N effect by I level": Obtain lincom results 
  fit <- lm(diff ~ N + anl + R + I + est + 
              N:anl + N:R + N:I + N:est + 
              anl:R + anl:I + anl:est + 
              R:I + R:est + 
              I:est, 
            data = df_sub_relevel)
  
  lincom <- get_lincom_results(fit, df_sub_relevel, 
                               "anl", "I", ref_anl, ref_I, 
                               vec_ref)
  
  ### Save list objects 
  list_lincom[[i]] <- lincom
  list_df_plot[[i]] <- lincom[[2]] %>%
    mutate(elem1 = elem_truth, elem2 = elem_est, elem3 = elem_tail_p) %>%
    dplyr::select(elem1, elem2, elem3, everything())
  
}


### Bind rows and convert elements to factors
df_plot <- bind_rows(list_df_plot) %>%
  mutate(elem1 = factor(elem1, levels = vec_truth), 
         elem2 = factor(elem2, levels = vec_est), 
         elem3 = factor(elem3, levels = vec_tail_p))

classmode(df_plot, everything())


### Relabel ovars
tabdf(df_plot, ovar1)
tabdf(df_plot, ovar2)

df_plot <- df_plot %>%
  mutate(ovar1 = factor(ovar1, 
                        labels = c("DP-diffuse", "DP-inform")), 
         ovar2 = factor(ovar2, 
                        labels = c("I = 0.1", "I = 0.5", "I = 0.9")),
         elem1 = factor(elem1, 
                        labels = c("ALD", "Gaussian Mixture")), 
         elem2 = factor(elem2, 
                        labels = c("PM", "CB", "GR")))


### Plot A: The effect of DP priors on the bias of 50 percentile
df_plot_sub <- df_plot %>% 
  filter(elem3 == 50)

pA <- ggplot(data = df_plot_sub, 
             aes(x = ovar2, y = estimate, 
                 group = ovar1, color = ovar1, shape = ovar1)) +
  geom_point(aes(y = estimate), size = 2.5, position = position_dodge(width = 0)) + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(width = 0), width = 0.2) + 
  geom_line(aes(y = estimate)) + 
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed") + 
  facet_grid(elem1 ~ elem2) + 
  # geom_text(aes(label = round(estimate, 0)), color = "black", 
  #           position = position_stack(vjust = 0.1, reverse = TRUE), size = 3) + 
  # scale_y_continuous(labels = comma, 
  #                    limits = c(-800, 200), 
  #                    breaks = seq(-800, 200, by = 200)) + 
  scale_color_manual(values = color_palette[seq(unique(df_plot_sub$ovar1))]) +
  scale_shape_manual(values = shape_palette[seq(unique(df_plot_sub$ovar1))]) +
  theme_trend +
  theme(strip.background = element_rect(fill = "gray100")) + 
  theme(legend.position = "right", legend.direction = "vertical") +  
  labs(title = "The effect of DP priors on the bias of 50 percentile estimate", 
       # subtitle = "Estimator: PM", 
       y = "Effect estimate", x = NULL)



### Plot B: The effect of N on the bias of 90 percentile (GR)
df_plot_sub <- df_plot %>% 
  filter(elem3 == 90)

pB <- ggplot(data = df_plot_sub, 
             aes(x = ovar2, y = estimate, 
                 group = ovar1, color = ovar1, shape = ovar1)) +
  geom_point(aes(y = estimate), size = 2.5, position = position_dodge(width = 0)) + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(width = 0), width = 0.2) + 
  geom_line(aes(y = estimate)) + 
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed") + 
  facet_grid(elem1 ~ elem2) + 
  # geom_text(aes(label = round(estimate, 0)), color = "black", 
  #           position = position_stack(vjust = 0.1, reverse = TRUE), size = 3) + 
  # scale_y_continuous(labels = comma, 
  #                    limits = c(-800, 200), 
  #                    breaks = seq(-800, 200, by = 200)) + 
  scale_color_manual(values = color_palette[seq(unique(df_plot_sub$ovar1))]) +
  scale_shape_manual(values = shape_palette[seq(unique(df_plot_sub$ovar1))]) +
  theme_trend +
  theme(strip.background = element_rect(fill = "gray100")) + 
  theme(legend.position = "right", legend.direction = "vertical") +  
  labs(title = "The effect of DP priors on the bias of 90 percentile estimate", 
       # subtitle = "Estimator: PM", 
       y = "Effect estimate", x = NULL)



###'######################################################################
###'
###' Combine all two plots into one 
###' 
###'

p_all <- cowplot::plot_grid(pA, pB, ncol = 1, labels = "AUTO")

setwd(work_dir)
ggsave("figures/Figure09_Effects of DP priors on percentile estimates.pdf", 
       p_all, width = 8, height = 9)



