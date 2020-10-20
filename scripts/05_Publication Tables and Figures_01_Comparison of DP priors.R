
###'######################################################################
###'
###' Category: Publication Tables and Figures 
###' 
###' Task: [Figure 02]. Comparison of DP priors
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
library(gmp)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "functions/03_Simulation Implementation_04_Plot helpers.R")
source(file = "functions/03_Simulation Implementation_05_Data management helpers.R")


### Set ggplot themes
### Theme settings
theme_preset <- theme_bw() + 
  theme(panel.background = element_blank(),
        panel.grid = element_blank(), 
        legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_blank())



###'######################################################################
###'
###' get_induced_K_prior():
###'
###' - Define a function to obtain the alpha-induced prior of K 
###'   (the number of clusters)
###'
###'

### Define a function to generate the alpha induced prior of K
get_induced_K_prior <- function(alpha, N, a, b){
  
  ### Define an empty vector to store pobability mass
  temp_vec <- vector()
  
  for (k in seq(N)){
    
    ### Code the constant part using "unsigned" Sirling number of the first kind
    sign <- (-1)^(N - k)
    constant <- sign*(b^a*as.numeric(gmp::Stirling1(n = N, k = k)))/gamma(a)
    
    ### Code the numerical integration part
    integrand <- function(alpha, N, k, a, b){
      
      log_comp1 <- (k + a - 1)*log(alpha)
      log_comp2 <- -b*alpha
      log_comp3 <- lgamma(alpha)
      log_comp4 <- -lgamma(alpha + N)
      
      ifelse(alpha == 0, 0, exp(log_comp1 + log_comp2 + log_comp3 + log_comp4))
    }
    
    integral <- integrate(f = integrand, lower = 0, upper = Inf, 
                          N = N, k = k, a = a, b = b)$value
    
    ### Store in the empty vector
    temp_vec[k] <- constant*integral
  }
  
  return(temp_vec)
}


### Get example quantity
# N <- 50; a <- 0.467; b <- 0.007
N <- 50; a <- 2.5; b <- 0.1 
induced_K_prior_ex <- get_induced_K_prior(alpha = alpha, N = N, a = a, b = b)
plot(induced_K_prior_ex)



###'######################################################################
###'
###' Simulate an alpha vector 
###' 
###' alpha ~ Gamma(0.467, 0.007)
###' 
###' a, b values come from Dorazio (2009) Table 1
###' 
###'

### Set shape and rate parameters for Gamma distribution
a <- 0.467  # shape parameter
b <- 0.007  # rate parameter


### Calculate a priori values for the mean and variance of alpha
pri_mean <- a/b
pri_var <- a/b^2


### Generate alpha from the Gamma distribution
set.seed(1234567)
alpha <- rgamma(n = 2000, shape = a, rate = b)


### Plot the empirical distribution of alpha
ggplot(data = data.frame(alpha), aes(x = alpha)) + 
  geom_density(fill = "red", alpha = 0.2) + 
  theme_minimal()


### Compare simulated alpha's mean and variance with those of a priori
mean(alpha); pri_mean
var(alpha); pri_var



###'######################################################################
###'
###' [Panel A]. Comparison of 
###'   
###' (1) Chi-squared prior on K and 
###' (2) alpha-induced K (solved by KLD)   
###'
###'

### Generate (1) and (2)
N <- 50
a <- 1.60
b <- 1.22
df <- 5

uniform <- rep(1/N, N)

induced <- get_induced_K_prior(alpha = alpha, N, a, b)

a_priori <- dchisq(seq(N), df = df, ncp = 0)

df_temp <- data.frame(seq(N), uniform, induced, a_priori)

names(df_temp) <- c("K", "Uniform_K",  "Induced_K", "A_priori_K")


### Convert to the long format
df_plot <- df_temp %>%
  gather(key = K_Prior, value = value, -K) %>%
  mutate(K_Prior = factor(K_Prior, 
                          levels = c("Uniform_K", "A_priori_K", "Induced_K"), 
                          labels = c("Uniform K", "a priori K", 
                                     "Induced K by Gamma(1.60, 1.22)")))


### Plot!
p1 <- ggplot(data = df_plot, 
            aes(x = K, y = value, 
                group = K_Prior, color = K_Prior, linetype = K_Prior)) + 
  geom_line(size = 0.7) + 
  # geom_vline(xintercept = df, linetype = "longdash", color = "gray70") + 
  scale_x_continuous(breaks = seq(from = 0, to = N, by = 5)) + 
  scale_color_manual(values = rev(color_palette[seq(unique(df_plot$K_Prior))])) +
  scale_linetype_manual(values = c("dotted", "solid", "dashed")) +
  theme_preset + 
  labs(title = "The probability mass function for the expected number of clusters (K)", 
       subtitle = expression(paste("A priori ", chi^2, "(5) vs. the prior of K induced by the informative Gamma(1.60, 1.22) prior for ", alpha)), 
       x = "The expected number of clusters (K)", 
       y = "Probability")



###'######################################################################
###'
###' [Panel B]. AB solution by grid search method
###'
###'
###'

### Import AB grid search results & tidy up the results
setwd("~/Bayes-deconvolution/datasets/09_Informative prior ab solutions")
df_result <- read.csv(file = "df_ab_solution_informative_N50.csv")


### Filter out -Inf
df_result <- df_result %>%
  filter(KLD != -Inf) %>%
  filter(KLD >= 0)


### Get ab solution with KLD measure
KLD_min <- min(df_result$KLD, na.rm = TRUE)

vec_ab_solution <- df_result %>%
  filter(KLD == KLD_min) %>%
  unlist()


### Plot the solution 
df_plot <- df_result %>%
  filter(a < 3.0 & b < 3.0)

p2 <- ggplot(data = df_plot, aes(x = a, y = b, z = KLD)) +
  geom_raster(aes(fill = KLD)) +
  scale_fill_viridis_c(direction = 1) +  
  geom_contour(color = "white") +
  geom_text_contour(colour = "white") + 
  stat_subset(aes(subset = abs(KLD) < abs(KLD_min)*2), color = "red", size = 0.5) + 
  scale_x_continuous(breaks = seq(0, 3, by = 0.5)) + 
  scale_y_continuous(breaks = seq(0, 3, by = 0.5)) + 
  theme_bw() + 
  labs(title = paste0("The Kullback-Leibler divergence measures (N = ", N, ")"), 
       subtitle = paste0("Solution: ", 
                         "a = 1.60, ", 
                         "b = 1.22, ", 
                         "KLD = 0.002"), 
       x = "Shape parameter (a) of the Gamma prior", 
       y = "Rate parameter (b) of the Gamma prior", 
       caption = NULL)

p2


###'######################################################################
###'
###' Combine [Panel A] & [Panel B]
###'
###'

### Combine
p_row1 <- plot_grid(p2, p1, nrow = 1, labels = c("A", "B"), rel_widths = c(1, 1.2))

p_row1


### Save the resulting plot 
setwd(work_dir)
ggsave(filename = "figures/Figure03_Derivation of informative prior on alpha.pdf", 
       p_row1, width = 12, height = 5)



###'######################################################################
###'
###' Informative vs. Diffuse prior on the K
###'
###'

### Generate a dataframe
uniform <- rep(1/N, N)

diffuse <- get_induced_K_prior(alpha = alpha, N = 50, 
                               a = 2.5, b = 0.1)

informative <- dchisq(seq(N), df = df, ncp = 0)

df_temp <- data.frame(seq(N), uniform, diffuse, informative)

names(df_temp) <- c("K", "Uniform",  "Diffuse", "Informative")


### Convert to the long format
df_plot <- df_temp %>%
  gather(key = K_Prior, value = value, -K) %>%
  mutate(K_Prior = factor(K_Prior, 
                          levels = c("Uniform", "Diffuse", "Informative"), 
                          labels = c("Uniform", "Diffuse", "Informative")))


### Plot!
p3 <- ggplot(data = df_plot, 
             aes(x = K, y = value, 
                 group = K_Prior, color = K_Prior, linetype = K_Prior)) + 
  geom_line(size = 0.5) + 
  # geom_vline(xintercept = df, linetype = "longdash", color = "gray70") + 
  scale_x_continuous(breaks = seq(from = 0, to = N, by = 5)) + 
  scale_color_manual(values = c("forestgreen", "firebrick1", "dodgerblue1")) +
  scale_linetype_manual(values = c("dotted", "longdash", "solid")) +
  theme_preset + 
  labs(title = "The probability mass function for the expected number of clusters (K)", 
       subtitle = expression(paste("The prior of K induced by the diffuse Gamma(2.5, 0.1) prior for ", alpha, " vs. the informative ", chi^2, "(5)")), 
       x = "The expected number of clusters (K)", 
       y = "Probability")


### Save the resulting plot 
setwd(work_dir)
ggsave(filename = "figures/Figure02_Diffuse vs Informative K priors.pdf", 
       p3, width = 8, height = 5)



###'######################################################################
###'
###' Combine ALL
###' 
###'

### Combine
p_all <- plot_grid(p_row1, p3, ncol = 1, labels = c("", "A"))

p_all


### Save the resulting plot 
setwd(work_dir)
ggsave(filename = "figures/Figure02_Derivation of informative prior on alpha.pdf", 
       p_row1, width = 12, height = 5)

