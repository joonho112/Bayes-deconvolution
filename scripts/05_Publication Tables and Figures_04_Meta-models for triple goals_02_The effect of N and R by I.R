
###'######################################################################
###'
###' Category: Publication Tables and Figures 
###' 
###' Task: Interaction plots
###' 
###'  [Figure 05]. Interaction plots for Meta-model regression results 
###'  
###'   - MSEL for individual site effects & ISEL for EDF (MSELR is missing)
###'  
###'   - The effect of N and R by I levels
###'   
###'           
###' Data: Simulated data
###' 
###' Data: 2020-04-27
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
# data_dir <- file.path(work_dir, "datasets", 
#                       "11_All collected and replicated loss estimates")

data_dir <- file.path(work_dir, "datasets", 
                      "11-1_All collected and replicated loss estimates_updated-MSELR")



### Call libraries
library(tidyverse)
library(cowplot)
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
         SELrank = SELrank)

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
###' [Figure 05]. 
###' 
###' [Panel A] Effects of increasing N on the MSEL of individual effects (PM)
###' 
###' By shrinkage factor I
###' 
###'

### Prepare a loop
list_lincom <- list()
list_df_plot <- list()

vec_truth <- c("ALD", "Mixed")
vec_anl <- c("Gaussian", "DP_diffuse", "DP_Inform")
df_grid <- expand.grid(vec_truth, vec_anl) %>%
  data.frame() %>% set_names(c("truth", "anl"))


### Start loop 
for (i in seq(nrow(df_grid))){
  
  ### Extract loop elements
  elem_truth <- df_grid$truth[i]
  elem_anl <- df_grid$anl[i] %>% as.character()
  
  
  ### Filter only conditional cases
  df_sub <- df %>%
    filter(truth %in% elem_truth) %>%
    filter(est != c("ML")) %>%
    filter(anl != c("T")) %>%
    filter(!(N %in% c("N100"))) 
  
  
  ### Set reference category for each factor
  ref_N <- "N25"
  ref_I <- "I0.1"
  ref_R <- "R1"
  ref_anl <- elem_anl
  ref_est <- "PM"
  
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
  fit <- lm(SSEL ~ N + I + R + anl + est + 
              N:I + N:R + N:anl + N:est + 
              I:R + I:anl + I:est + 
              R:anl + R:est + 
              anl:est, 
            data = df_sub_relevel)
  
  lincom <- get_lincom_results(fit, df_sub_relevel, 
                               "N", "I", ref_N, ref_I, 
                                vec_ref)
  
  ### Save list objects 
  list_lincom[[i]] <- lincom
  list_df_plot[[i]] <- lincom[[2]] %>%
    mutate(elem1 = elem_truth, elem2 = elem_anl) %>%
    dplyr::select(elem1, elem2, everything())
  
}


### Bind rows and convert elements to factors
df_plot <- bind_rows(list_df_plot) %>%
  mutate(elem1 = factor(elem1, levels = vec_truth), 
         elem2 = factor(elem2, levels = vec_anl))
  
classmode(df_plot, everything())


### Relabel ovars
tabdf(df_plot, ovar1)
tabdf(df_plot, ovar2)

df_plot <- df_plot %>%
  mutate(ovar1 = factor(ovar1, 
                        labels = c("N = 50", "N = 200")), 
         ovar2 = factor(ovar2, 
                        labels = c("I = 0.1", "I = 0.5", "I = 0.9")), 
         elem2 = factor(elem2, 
                        labels = c("Gaussian", "DP-diffuse", "DP-inform")))


### Plot!
df_plot_sub <- df_plot %>% 
  filter(elem1 == "Mixed")

pA <- ggplot(data = df_plot_sub, 
             aes(x = ovar2, y = estimate, 
                 group = ovar1, color = ovar1, shape = ovar1)) +
  geom_point(aes(y = estimate), size = 2.5, position = position_dodge(width = 0)) + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(width = 0), width = 0.2) + 
  geom_line(aes(y = estimate)) + 
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed") + 
  facet_grid(. ~ elem2) + 
  # geom_text(aes(label = round(estimate, 0)), color = "black", 
  #           position = position_stack(vjust = 0.1, reverse = TRUE), size = 3) + 
  scale_y_continuous(labels = comma, 
                     limits = c(-800, 200), 
                     breaks = seq(-800, 200, by = 200)) + 
  scale_color_manual(values = color_palette[seq(unique(df_plot_sub$ovar1))]) +
  scale_shape_manual(values = shape_palette[seq(unique(df_plot_sub$ovar1))]) +
  theme_trend +
  theme(strip.background = element_rect(fill = "gray100")) + 
  theme(legend.position = "right", legend.direction = "vertical") +  
  labs(title = "The effect of N on the MSEL of site-specific effects (PM)", 
       # subtitle = "Estimator: PM", 
       y = "Effect estimate", x = NULL)



###'######################################################################
###'
###' [Figure 05]. 
###' 
###' [Panel B] Effects of increasing N on the ISEL of EDF (GR)
###' 
###' By shrinkage factor I
###' 
###'

### Prepare a loop
list_lincom <- list()
list_df_plot <- list()

vec_truth <- c("ALD", "Mixed")
vec_anl <- c("Gaussian", "DP_diffuse", "DP_Inform")
df_grid <- expand.grid(vec_truth, vec_anl) %>%
  data.frame() %>% set_names(c("truth", "anl"))


### Start loop 
for (i in seq(nrow(df_grid))){
  
  ### Extract loop elements
  elem_truth <- df_grid$truth[i]
  elem_anl <- df_grid$anl[i] %>% as.character()
  
  
  ### Filter only conditional cases
  df_sub <- df %>%
    filter(truth %in% elem_truth) %>%
    filter(est != c("ML")) %>%
    filter(anl != c("T")) %>%
    filter(!(N %in% c("N100"))) 
  
  
  ### Set reference category for each factor
  ref_N <- "N25"
  ref_I <- "I0.1"
  ref_R <- "R1"
  ref_anl <- elem_anl
  ref_est <- "GR"
  
  vec_ref <- paste(elem_truth, ref_N, ref_I, ref_R, 
                   ref_anl, ref_est, sep = ", ")
  
  df_sub_relevel <- df_sub %>%
    mutate(N = relevel(N, ref = ref_N), 
           I = relevel(I, ref = ref_I), 
           R = relevel(R, ref = ref_R), 
           anl = relevel(anl, ref = ref_anl), 
           est = relevel(est, ref = ref_est))
  
  
  ###' (2) ISEL for EDFs
  ###' "N effect by I level": Obtain lincom results 
  fit <- lm(ISEL ~ N + I + R + anl + est + 
              N:I + N:R + N:anl + N:est + 
              I:R + I:anl + I:est + 
              R:anl + R:est + 
              anl:est, 
            data = df_sub_relevel)
  
  lincom <- get_lincom_results(fit, df_sub_relevel, 
                               "N", "I", ref_N, ref_I, 
                               vec_ref)
  
  ### Save list objects 
  list_lincom[[i]] <- lincom
  list_df_plot[[i]] <- lincom[[2]] %>%
    mutate(elem1 = elem_truth, elem2 = elem_anl) %>%
    dplyr::select(elem1, elem2, everything())
  
}


### Bind rows and convert elements to factors
df_plot <- bind_rows(list_df_plot) %>%
  mutate(elem1 = factor(elem1, levels = vec_truth), 
         elem2 = factor(elem2, levels = vec_anl))

classmode(df_plot, everything())


### Relabel ovars
tabdf(df_plot, ovar1)
tabdf(df_plot, ovar2)

df_plot <- df_plot %>%
  mutate(ovar1 = factor(ovar1, 
                        labels = c("N = 50", "N = 200")), 
         ovar2 = factor(ovar2, 
                        labels = c("I = 0.1", "I = 0.5", "I = 0.9")), 
         elem2 = factor(elem2, 
                        labels = c("Gaussian", "DP-diffuse", "DP-inform")))


### Plot!
df_plot_sub <- df_plot %>% 
  filter(elem1 == "Mixed")

pB <- ggplot(data = df_plot_sub, 
             aes(x = ovar2, y = estimate, 
                 group = ovar1, color = ovar1, shape = ovar1)) +
  geom_point(aes(y = estimate), size = 2.5, position = position_dodge(width = 0)) + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(width = 0), width = 0.2) + 
  geom_line(aes(y = estimate)) + 
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed") + 
  facet_grid(. ~ elem2) + 
  # geom_text(aes(label = round(estimate, 0)), color = "black", 
  #           position = position_stack(vjust = 0.1, reverse = TRUE), size = 3) + 
  scale_y_continuous(labels = comma, 
                     limits = c(-100, 20), 
                     breaks = seq(-100, 20, by = 20)) + 
  scale_color_manual(values = color_palette[seq(unique(df_plot_sub$ovar1))]) +
  scale_shape_manual(values = shape_palette[seq(unique(df_plot_sub$ovar1))]) +
  theme_trend +
  theme(strip.background = element_rect(fill = "gray100")) + 
  theme(legend.position = "right", legend.direction = "vertical") + 
  labs(title = "The effect of N on the ISEL of EDF (GR)", 
       # subtitle = "Estimator: GR", 
       y = "Effect estimate", x = NULL)



###'######################################################################
###'
###' [Figure 05]. 
###' 
###' [Panel C] Effects of R on the MSEL of individual effects (PM)
###' 
###' By shrinkage factor I
###' 
###'

### Prepare a loop
list_lincom <- list()
list_df_plot <- list()

vec_truth <- c("ALD", "Mixed")
vec_anl <- c("Gaussian", "DP_diffuse", "DP_Inform")
df_grid <- expand.grid(vec_truth, vec_anl) %>%
  data.frame() %>% set_names(c("truth", "anl"))


### Start loop 
for (i in seq(nrow(df_grid))){
  
  ### Extract loop elements
  elem_truth <- df_grid$truth[i]
  elem_anl <- df_grid$anl[i] %>% as.character()
  
  
  ### Filter only conditional cases
  df_sub <- df %>%
    filter(truth %in% elem_truth) %>%
    filter(est != c("ML")) %>%
    filter(anl != c("T")) %>%
    filter(!(N %in% c("N100"))) 
  
  
  ### Set reference category for each factor
  ref_N <- "N25"
  ref_I <- "I0.1"
  ref_R <- "R1"
  ref_anl <- elem_anl
  ref_est <- "PM"
  
  vec_ref <- paste(elem_truth, ref_N, ref_I, ref_R, 
                   ref_anl, ref_est, sep = ", ")
  
  df_sub_relevel <- df_sub %>%
    mutate(N = relevel(N, ref = ref_N), 
           I = relevel(I, ref = ref_I), 
           R = relevel(R, ref = ref_R), 
           anl = relevel(anl, ref = ref_anl), 
           est = relevel(est, ref = ref_est))
  
  
  ###' (1) SSEL for individual site-specific effects
  ###' "R effect by I level": Obtain lincom results 
  fit <- lm(SSEL ~ N + R + I + anl + est + 
              N:R + N:I + N:anl + N:est + 
              R:I + R:anl + R:est +
              I:anl + I:est + 
              anl:est, 
            data = df_sub_relevel)
  
  lincom <- get_lincom_results(fit, df_sub_relevel, 
                               "R", "I", ref_R, ref_I, 
                               vec_ref)
  
  ### Save list objects 
  list_lincom[[i]] <- lincom
  list_df_plot[[i]] <- lincom[[2]] %>%
    mutate(elem1 = elem_truth, elem2 = elem_anl) %>%
    dplyr::select(elem1, elem2, everything())
  
}


### Bind rows and convert elements to factors
df_plot <- bind_rows(list_df_plot) %>%
  mutate(elem1 = factor(elem1, levels = vec_truth), 
         elem2 = factor(elem2, levels = vec_anl))

classmode(df_plot, everything())


### Relabel ovars
tabdf(df_plot, ovar1)
tabdf(df_plot, ovar2)

df_plot <- df_plot %>%
  mutate(ovar1 = factor(ovar1, 
                        labels = c("R = 5", "R = 10")), 
         ovar2 = factor(ovar2, 
                        labels = c("I = 0.1", "I = 0.5", "I = 0.9")), 
         elem2 = factor(elem2, 
                        labels = c("Gaussian", "DP-diffuse", "DP-inform")))


### Plot!
df_plot_sub <- df_plot %>% 
  filter(elem1 == "Mixed")

pC <- ggplot(data = df_plot_sub, 
             aes(x = ovar2, y = estimate, 
                 group = ovar1, color = ovar1, shape = ovar1)) +
  geom_point(aes(y = estimate), size = 2.5, position = position_dodge(width = 0)) + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(width = 0), width = 0.2) + 
  geom_line(aes(y = estimate)) + 
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed") + 
  facet_grid(. ~ elem2) + 
  # geom_text(aes(label = round(estimate, 0)), color = "black", 
  #           position = position_stack(vjust = 0.1, reverse = TRUE), size = 3) + 
  scale_y_continuous(labels = comma, 
                     limits = c(-800, 200), 
                     breaks = seq(-800, 200, by = 200)) + 
  scale_color_manual(values = color_palette[seq(unique(df_plot_sub$ovar1))]) +
  scale_shape_manual(values = shape_palette[seq(unique(df_plot_sub$ovar1))]) +
  theme_trend +
  theme(strip.background = element_rect(fill = "gray100")) + 
  theme(legend.position = "right", legend.direction = "vertical") + 
  labs(title = "The effect of R on the MSEL of site-specific effects (PM)", 
       # subtitle = "Estimator: PM", 
       y = "Effect estimate", x = NULL)



###'######################################################################
###'
###' [Figure 05]. 
###' 
###' [Panel D] Effects of R on the ISEL of EDF (GR)
###' 
###' By shrinkage factor I
###' 
###'

### Prepare a loop
list_lincom <- list()
list_df_plot <- list()

vec_truth <- c("ALD", "Mixed")
vec_anl <- c("Gaussian", "DP_diffuse", "DP_Inform")
df_grid <- expand.grid(vec_truth, vec_anl) %>%
  data.frame() %>% set_names(c("truth", "anl"))


### Start loop 
for (i in seq(nrow(df_grid))){
  
  ### Extract loop elements
  elem_truth <- df_grid$truth[i]
  elem_anl <- df_grid$anl[i] %>% as.character()
  
  
  ### Filter only conditional cases
  df_sub <- df %>%
    filter(truth %in% elem_truth) %>%
    filter(est != c("ML")) %>%
    filter(anl != c("T")) %>%
    filter(!(N %in% c("N100"))) 
  
  
  ### Set reference category for each factor
  ref_N <- "N25"
  ref_I <- "I0.1"
  ref_R <- "R1"
  ref_anl <- elem_anl
  ref_est <- "GR"
  
  vec_ref <- paste(elem_truth, ref_N, ref_I, ref_R, 
                   ref_anl, ref_est, sep = ", ")
  
  df_sub_relevel <- df_sub %>%
    mutate(N = relevel(N, ref = ref_N), 
           I = relevel(I, ref = ref_I), 
           R = relevel(R, ref = ref_R), 
           anl = relevel(anl, ref = ref_anl), 
           est = relevel(est, ref = ref_est))
  
  
  ###' (2) ISEL for EDFs
  ###' "R effect by I level": Obtain lincom results 
  fit <- lm(ISEL ~ N + R + I + anl + est + 
              N:R + N:I + N:anl + N:est + 
              R:I + R:anl + R:est +
              I:anl + I:est + 
              anl:est, 
            data = df_sub_relevel)
  
  lincom <- get_lincom_results(fit, df_sub_relevel, 
                               "R", "I", ref_R, ref_I, 
                               vec_ref)
  
  ### Save list objects 
  list_lincom[[i]] <- lincom
  list_df_plot[[i]] <- lincom[[2]] %>%
    mutate(elem1 = elem_truth, elem2 = elem_anl) %>%
    dplyr::select(elem1, elem2, everything())
  
}


### Bind rows and convert elements to factors
df_plot <- bind_rows(list_df_plot) %>%
  mutate(elem1 = factor(elem1, levels = vec_truth), 
         elem2 = factor(elem2, levels = vec_anl))

classmode(df_plot, everything())


### Relabel ovars
tabdf(df_plot, ovar1)
tabdf(df_plot, ovar2)

df_plot <- df_plot %>%
  mutate(ovar1 = factor(ovar1, 
                        labels = c("R = 5", "R = 10")), 
         ovar2 = factor(ovar2, 
                        labels = c("I = 0.1", "I = 0.5", "I = 0.9")), 
         elem2 = factor(elem2, 
                        labels = c("Gaussian", "DP-diffuse", "DP-inform")))


### Plot!
df_plot_sub <- df_plot %>% 
  filter(elem1 == "Mixed")

pD <- ggplot(data = df_plot_sub, 
             aes(x = ovar2, y = estimate, 
                 group = ovar1, color = ovar1, shape = ovar1)) +
  geom_point(aes(y = estimate), size = 2.5, position = position_dodge(width = 0)) + 
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                position = position_dodge(width = 0), width = 0.2) + 
  geom_line(aes(y = estimate)) + 
  geom_hline(yintercept = 0, color = "gray30", linetype = "dashed") + 
  facet_grid(. ~ elem2) + 
  # geom_text(aes(label = round(estimate, 0)), color = "black", 
  #           position = position_stack(vjust = 0.1, reverse = TRUE), size = 3) + 
  scale_y_continuous(labels = comma, 
                     limits = c(-100, 20), 
                     breaks = seq(-100, 20, by = 20)) + 
  scale_color_manual(values = color_palette[seq(unique(df_plot_sub$ovar1))]) +
  scale_shape_manual(values = shape_palette[seq(unique(df_plot_sub$ovar1))]) +
  theme_trend +
  theme(legend.position = "right", legend.direction = "vertical") + 
  theme(strip.background = element_rect(fill = "gray100")) + 
  labs(title = "The effect of R on the ISEL of EDF (GR)", 
       # subtitle = "Estimator: GR", 
       y = "Effect estimate", x = NULL)



###'######################################################################
###'
###' Combine all four plots into one 
###' 
###'

p_all <- cowplot::plot_grid(pA, pB, pC, pD, ncol = 1, labels = "AUTO")

setwd(work_dir)
ggsave("figures/Figure05_Effects of N and R on loss estimates_01_Mixed.pdf", 
       p_all, width = 8, height = 10)


