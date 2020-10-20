
###'######################################################################
###'
###' Category: 
###' 
###' Task: (1) Dirichlet Process Mixtures
###'           
###'        Selecting a Prior for precision parameter alpha
###'        
###'        1. Adopting a prior for alpha
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
gc()            # force R to release memory it is no longer using
rm(list=ls())   # delete all the objects in the workspace


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


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)


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
###' Generate a simulated dataset 
###' 
###' from a mixture of two Gaussian homoscedastic components
###' 
###' (1) Get true PDF of thetas
###'
###'

### Set parameters
tailp <- c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95)  # tail probabilities
delta <- 4    # distance between two mixtures
ups <- 1      # want components to have equal variance
eps <- 0.2    # mixing proportion (smaller side)


### Generate a random vector indicating components 
ind <- runif(20000) < (1 - eps)
prop.table(table(ind))


### Define a normalizing factor `a`
a <- sqrt((1 - eps) + eps*ups^2 + eps*(1 - eps)*delta^2)
a


### Get true PDF of thetas from mixture
quant <- ind*rnorm(20000, -eps*delta/a, sqrt(1/a^2)) + 
  (1 - ind)*rnorm(20000, (1 - eps)*delta/a, sqrt(ups^2/a^2))


### Get true tail quantiles
mdppriquan <- quantile(quant, tailp)
mdppriquan


### Check mean, SD, and quantiles of the true PDF
c(mean(quant), sd(quant))



###'######################################################################
###'
###' Generate a simulated dataset 
###' 
###' from a mixture of two Gaussian homoscedastic components
###' 
###' (2) Define a function to generate random draws from the mixture distribution
###'
###'

GenerateMdp <- function(n, rgt, rsr){
  
  ### Set parameters 
  data_obj <- list()
  data_obj <- NULL
  delta <- 4    # distance between two mixtures
  ups <- 1      # want components to have equal variance
  eps <- 0.2    # mixing proportion (smaller side)
  
  
  ### Define a normalizing factor `a`
  a <- sqrt((1 - eps) + eps*ups^2 + eps*(1 - eps)*delta^2)
  
  
  ### Simulate a mixture of 2 normals with mean 0 and var 1  with this parameterization
  ind <- runif(n) < (1 - eps)
  
  data_obj$tau_k <- ind*rnorm(n, -eps*delta/a, sqrt(1/a^2)) + 
    (1 - ind)*rnorm(n, (1 - eps)*delta/a, sqrt(ups^2/a^2))
  
  
  ### Simulate true distribution for quantiles etc.
  ind <- runif(n) < (1 - eps)
  
  data_obj$tau_k <- ind*rnorm(n, -eps*delta/a, sqrt(1/a^2)) + 
    (1 - ind)*rnorm(n, (1 - eps)*delta/a, sqrt(ups^2/a^2))
  
  tailp <- c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95)
  
  
  ### Generate within-site variances (sds) for K sites
  sigma2max <- rgt*rsr
  sigma2min <- sigma2max/(rsr^2)
  sigma2 <- exp(seq(from = log(sigma2min), to = log(sigma2max), length = n))
  sd <- sqrt(sigma2)
  
  
  ### Generate data object (observed Y)
  data_obj$tau_k_hat <- rnorm(n, data_obj$tau_k, sd)
  data_obj$se_k2 <- sigma2
  data_obj$priquan <- mdppriquan
  return(data_obj)
}



###'######################################################################
###'
###' Simulate an example data
###'
###'

### Generate sample draws from the mixture model
set.seed(12345)

df_G <- GenerateMdp(n = 50, rgt = 0.33, rsr = 5)[c(1:3)] %>%
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
sd(alpha); pri_var



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
N <- 50; a <- 0.467; b <- 0.007
induced_K_prior_ex <- get_induced_K_prior(alpha = alpha, N = N, a = a, b = b)



###'######################################################################
###'
###' (Visually) Compare to the uniform prior for the K (1/N)
###' 
###' 

df_temp <- data.frame(seq(N), rep(1/N, N), induced_K_prior_ex)
names(df_temp) <- c("K", "Uniform", "Induced")

df_plot <- df_temp %>%
  gather(key = K_Prior, value = value, -K) %>%
  mutate(K_Prior = factor(K_Prior, levels = c("Uniform", "Induced")))


p <- ggplot(data = df_plot, aes(x = K, y = value, 
                                group = K_Prior, color = K_Prior)) + 
  geom_line(size = 1) + 
  scale_x_continuous(breaks = seq(from = 0, to = 50, by = 2)) + 
  labs(title = "The probability mass function for the uniform vs. the alpha-induced priors for K", 
       subtitle = "The priors of K induced by Gamma(0.467, 0.007) for the precision parameter", 
       x = "The number of expected number of clusters (K)", 
       y = "Probability") + theme_preset

ggsave(filename = "figures/Uniform vs Induced K priors.pdf", p, 
       width = 8, height = 6)



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
pi_x <- df_temp$Uniform
pi_y <- df_temp$Induced

### Define a function to calculate the Kullback-Leibler divergence measure 
get_KLD <- function(pi_x, pi_y){
  sum(pi_x*(log(pi_x) - log(pi_y)))
  }


### Get KLD measure (be careful about the order of priors)
get_KLD(df_temp$Uniform, df_temp$Induced)
get_KLD(df_temp$Induced, df_temp$Uniform)


### KLD measure from the reduced form
-log(N) -(1/N)*sum(log(pi_y))



###'######################################################################
###'
###' Define a function to be used in Gamma(a, b) optimization problem
###'
###'

N <- 200

ab_KLD <- function(a, b){
  
  ### Get draws from Gamma prior & set the number of cases
  alpha <- rgamma(n = 2000, shape = a, rate = b)
  N <- 200
  
  ### Obtain the alpha induced prior of K (the number of clusters)
  pi_y <- get_induced_K_prior(alpha = alpha, N = N, a = a, b = b)
  
  ### Return KLD measure from the reduced form
  -log(N) -(1/N)*sum(log(pi_y)) 
}



###'######################################################################
###'
###' Multidimensional optimization functions don't work 
###' 
###' Let's first work with input grids
###'
###'

### Define an input grid
n_breaks <- 501
x <- seq(0, 2, length = n_breaks)

N_vec <- rep(N, length(x)*n_breaks)
a <- rep(x, each = n_breaks)
b <- rep(x, times = n_breaks)

input_list <- list(N_vec, a, b) 
input_mat <- cbind(a, b)


### Apply the function across possible pairs
safe_ab_KLD <- safely(ab_KLD, otherwise = NA_real_)

list_result <- map2(input_mat[,1], input_mat[,2], safe_ab_KLD)

saveRDS(list_result,
        file = paste0("datasets/list_ab_solution_N", N, ".rds"))


### Tidy up the results
df_KLD_values <- map(list_result, "result") %>% unlist()

df_result <- cbind.data.frame(input_mat, df_KLD_values) %>%
  arrange(df_KLD_values)

names(df_result) <- c("a", "b", "KLD")

write.csv(df_result, file = paste0("datasets/df_ab_solution_N", N, ".csv"))


### Get ab solution with KLD measure
KLD_min <- min(df_result$KLD, na.rm = TRUE)

vec_ab_solution <- df_result %>%
  filter(KLD == KLD_min) %>%
  unlist()


### Plot the solution 
p <- ggplot(data = df_result, aes(x = a, y = b, z = KLD)) +
  geom_raster(aes(fill = KLD)) +
  scale_fill_viridis_c() +  
  geom_contour(color = "white") +
  geom_text_contour(colour = "white") + 
  stat_subset(aes(subset = KLD < KLD_min*3), color = "red") + 
  theme_bw() + 
  labs(title = paste0("The Kullback-Leibler divergence measure (N = ", N, ")"), 
       subtitle = "Assuming a uniform prior on the number of clusters (K)", 
       x = "Shape parameter (a) of the Gamma prior", 
       y = "Rate parameter (b) of the Gamma prior", 
       caption = paste0("Solution: ", 
                        "a = ", vec_ab_solution[1], ", ", 
                        "b = ", vec_ab_solution[2], ", ", 
                        "KLD = ", round(vec_ab_solution[3], 3)))

p

ggsave(filename = paste0("figures/KLD Contour plot_N", N, ".pdf"), p, 
       width = 6, height = 6)

