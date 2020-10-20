
###'######################################################################
###'
###' Category: 
###' 
###' Task: (1) Dirichlet Process Mixtures
###'           
###'        Selecting a Prior for precision parameter alpha
###'        
###'        using DORO method
###'           
###' Data: Simulated data
###' 
###' Data: 2020-03-19
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
library(hhsim)
library(DPpackage)
library(bspmma)
library(rstan)
library(gmp)
library(broom)
library(metR)
library(parallel)
library(foreach)
library(doAzureParallel)


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
###' Generate a configuration file called "credentials.json"
###' 
###'  in the working directory 
###'
###' - And populate this file with 
###'   my Batch and storage account names and keys
###'
###'

### Generate a configuration file
setwd(work_dir)
generateCredentialsConfig("cluster/credentials.json")


### Set the credentials for my current R session
setCredentials("cluster/credentials.json")



###'######################################################################
###'
###' Create a Batch pool    
###'
###'

### Generate a cluster configuration JSON file in the working directory
setwd(work_dir)
generateClusterConfig("cluster/cluster_grid_search.json")


### Create your cluster if it does not exist; this takes a few minutes
cluster <- makeCluster("cluster/cluster_grid_search.json") 


### Register your parallel backend 
registerDoAzureParallel(cluster) 


### Check that the nodes are running 
getDoParWorkers()



###'######################################################################
###'
###' Generate a simulated dataset 
###' 
###' from a mixture of two Gaussian homoscedastic components
###'
###'

### Set parameters
delta <- 4    # distance between two mixtures
ups <- 1      # want components to have equal variance
eps <- 0.2    # mixing proportion (smaller side)


### Generate sample draws from the mixture model
set.seed(12345)

df_G <- Gen_Mixture(n = 50, ICC = 0.5, rsr = 5, delta, eps, ups) %>%
  data.frame() %>%
  rownames_to_column("K") %>%
  mutate_at(.vars = c("K"), .funs = as.numeric)



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
###' Obtain the alpha induced prior of K (the number of clusters)
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
###' (Visually) Compare to the uniform prior for the K 
###' 
###' 

df_temp <- data.frame(seq(N), rep(1/N, N), induced_K_prior_ex)
names(df_temp) <- c("K", "Uniform", "Induced")

df_plot <- df_temp %>%
  gather(key = K_Prior, value = value, -K) %>%
  mutate(K_Prior = factor(K_Prior, levels = c("Uniform", "Induced")))

vec_subtitle <- paste0("The priors of K induced by Gamma(", 
                       a, ", ", b, ") for the precision parameter")


p <- ggplot(data = df_plot, aes(x = K, y = value, 
                                group = K_Prior, color = K_Prior)) + 
  geom_line(size = 1) + 
  scale_x_continuous(breaks = seq(from = 0, to = N, by = 2)) + 
  labs(title = "The probability mass function for the uniform vs. the alpha-induced priors for K", 
       subtitle = vec_subtitle, 
       x = "The number of expected number of clusters (K)", 
       y = "Probability") + theme_preset

vec_filename <- paste0("figures/Uniform vs Induced K priors_Gamma_", a, "_", b, ".pdf")

ggsave(filename = vec_filename, p, 
       width = 14, height = 6)



###'######################################################################
###'
###' (Visually) Compare to the INFORMATIVE prior for the K 
###' 
###' A decreasing function
###' 
###' 

### Try chi-squared density functions
dens <- dchisq(seq(N), df = 5, ncp = 0)
plot(dens)


### Generate a dataframe combining uniform, informative, induced
df <- 5

df_temp <- data.frame(seq(N), 
                      rep(1/N, N), 
                      dchisq(seq(N), df = df, ncp = 0), 
                      induced_K_prior_ex)

names(df_temp) <- c("K", "Uniform", "Informative", "Induced")

df_plot <- df_temp %>%
  gather(key = K_Prior, value = value, -K) %>%
  mutate(K_Prior = factor(K_Prior, levels = c("Uniform", "Informative", "Induced")))


vec_subtitle <- paste0("The priors of K induced by Gamma(", 
                       a, ", ", b, ") for the precision parameter")


p <- ggplot(data = df_plot, aes(x = K, y = value, 
                                group = K_Prior, color = K_Prior, linetype = K_Prior)) + 
  geom_line(size = 1) + 
  scale_x_continuous(breaks = seq(from = 0, to = N, by = 2)) + 
  labs(title = "The probability mass function for the uniform vs. the alpha-induced priors for K", 
       subtitle = vec_subtitle, 
       x = "The number of expected number of clusters (K)", 
       y = "Probability") + theme_preset

vec_filename <- paste0("figures/Uniform_Informative_Induced K priors_Gamma_", a, "_", b, ".pdf")

ggsave(filename = vec_filename, p, 
       width = 14, height = 6)



###'######################################################################
###'
###' Calculate the Kullback-Leibler divergence measure
###' 
###' between (1) The prior of K and (2) the alpha-induced prior for K
###'
###' pi_x: the “true” distribution of data, observations, or theoretical distribution
###' pi_y: a theory, model, or approximation of pi_x
###'

### Assign uniform and induced priors
pi_x <- df_temp$Informative
pi_y <- df_temp$Induced

### Define a function to calculate the Kullback-Leibler divergence measure 
get_KLD <- function(pi_x, pi_y){
  sum(pi_x*(log(pi_x) - log(pi_y)))
}


###' Get KLD measure 
###' (be careful about the order of priors: pi_x should be the first)
get_KLD(pi_x, pi_y)
get_KLD(pi_y, pi_x)


### KLD measure from the reduced form
-log(N) -(1/N)*sum(log(pi_y))



###'######################################################################
###'
###' Define a function to be used in Gamma(a, b) optimization problem
###' 
###' => Compare to the informative Chi-squared distribution
###'
###'

ab_KLD_info <- function(N, df, a, b){
  
  ### Get draws from Gamma prior & set the number of cases
  alpha <- rgamma(n = 2000, shape = a, rate = b)
  
  ### Obtain the alpha induced prior of K (the number of clusters)
  pi_y <- get_induced_K_prior(alpha = alpha, N = N, a = a, b = b)
  
  ### Obtain the informative chi-squared prior of K
  pi_x <- dchisq(seq(N), df = df, ncp = 0)
  
  ### Return KLD measure
  sum(pi_x*(log(pi_x) - log(pi_y)))
}



###'######################################################################
###'
###' Grid search
###' 
###' - Multidimensional optimization functions don't work 
###' 
###' - Let's work with input grids
###'
###'

### Set the number of sites & the (assumed) number of clusters
N <- 170
df <- 5


### Define an input grid
n_breaks <- 501
x0 <- seq(0, 10, length = n_breaks)
x <- x0[-1] 

a <- rep(x, each = length(x))
b <- rep(x, times = length(x))
N_vec <- rep(N, length(a))
df_vec <- rep(df, length(a))

input_list <- list(N_vec, df_vec, a, b) 
input_mat <- cbind(N_vec, df_vec, a, b)


###' Optimize runtime. 
###' Chunking allows running multiple iterations on a single R instance.
opt <- list(chunkSize = 1000) 

start_time <- Sys.time()  # Mark start time

list_result <- 
  foreach (i = seq(nrow(input_mat)), .options.azure = opt, .errorhandling = "pass") %dopar% {
  
    library(gmp)
    ab_KLD_info(input_mat[i, 1], input_mat[i, 2], 
                input_mat[i, 3], input_mat[i, 4])
  
}

### Check the computation time
end_time <- Sys.time()  
end_time - start_time


### Save the dataframe
saveRDS(list_result,
        file = paste0("datasets/list_ab_solution_informative_N", N, ".rds"))



###'######################################################################
###'
###' Get ab solution with KLD measure
###' 
###' 

### Tidy up the results
df_KLD_values <- as.numeric(sapply(list_result, `[[`, 1))

df_result <- cbind.data.frame(input_mat, df_KLD_values) %>%
  arrange(df_KLD_values)

names(df_result) <- c("N", "df", "a", "b", "KLD")

write.csv(df_result, file = paste0("datasets/df_ab_solution_informative_N", N, ".csv"))


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

p <- ggplot(data = df_plot, aes(x = a, y = b, z = KLD)) +
  geom_raster(aes(fill = KLD)) +
  scale_fill_viridis_c(direction = 1) +  
  geom_contour(color = "white") +
  geom_text_contour(colour = "white") + 
  stat_subset(aes(subset = abs(KLD) < abs(KLD_min)*3), color = "red") + 
  theme_bw() + 
  labs(title = paste0("The Kullback-Leibler divergence measure (N = ", N, ")"), 
       subtitle = "Assuming an informative prior on the number of clusters (K)", 
       x = "Shape parameter (a) of the Gamma prior", 
       y = "Rate parameter (b) of the Gamma prior", 
       caption = paste0("Solution: ", 
                        "a = ", vec_ab_solution[3], ", ", 
                        "b = ", vec_ab_solution[4], ", ", 
                        "KLD = ", round(vec_ab_solution[5], 3)))

p

ggsave(filename = paste0("figures/KLD Contour plot_informative_N", N, ".pdf"), p, 
       width = 6, height = 6)

