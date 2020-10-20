
###'######################################################################
###'
###' Category: 
###' 
###' Task: (1) Dirichlet Process Mixtures
###'           
###'        Selecting a Prior for precision parameter alpha
###'        
###'        3. EB estimate of alpha using importance sampling
###'           
###' Data: Simulated data
###' 
###' Data: 2020-03-23
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
library(diagis)


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

N <- 50
rgt <- 0.33
rsr <- 5

df_G <- GenerateMdp(n = N, rgt = rgt, rsr = rsr)[c(1:3)] %>%
  data.frame() %>%
  rownames_to_column("K") %>%
  mutate_at(.vars = c("K"), .funs = as.numeric)


### Prepare dataset as a matrix form (with y and sigma2)
mat_DPmeta <- df_G %>% select(tau_k_hat, se_k2) %>% as.matrix()



###'######################################################################
###'
###' Preset necessary objects 
###'
###'

### Set initial state (the current value of the parameters)
state <- NULL


### Set MCMC parameters
nburn <- 4000    # the number of burn-in scans
nsave <- 4000    # the total number of scans to be saved
nskip <- 20      # the thinning interval
ndisplay <- 100  # the number of saved scans to be displayed on screen

mcmc <- list(nburn = nburn, nsave = nsave,
             nskip = nskip, ndisplay = ndisplay)


### Define a function to tidy up the DPmeta objects
get_posterior_DPmeta <- function(outp = outp, nburn = 4000){
  
  # (1) Tidy up the "theta" output object
  df_theta <- outp$save.state$randsave %>%
    data.frame() %>% 
    dplyr::select(-Prediction) %>%
    slice(c((nburn/2 + 1):nburn))  # throw away burn-ins
  colnames(df_theta) <- paste0("tau_k[", seq(ncol(df_theta)), "]")
  
  # (2) Tidy up the G0 parameters (m, s2), alpha0, and N of clusters outputs
  df_hyperparm <- outp$save.state$thetasave %>%
    data.frame() %>%
    dplyr::select(-tau_k_hat) %>%
    rename(G0_mu = mu, G0_s2 = sigma2, Ncluster = ncluster, alpha0 = alpha) %>%
    dplyr::select(G0_mu, G0_s2, alpha0, Ncluster) %>%
    slice(c((nburn/2 + 1):nburn))  # throw away burn-ins
  
  # Return the resulting object
  cbind.data.frame(df_hyperparm, df_theta)
}



###'######################################################################
###'
###' Set different hyperparameter choices
###'
###'

### Set various hyperparameter options
list_hyperparm <- list(c(4, 4), 
                       c(2, 1),
                       c(4, 1), 
                       c(0.392, 0.004), 
                       c(100, 0.1))

get_character_Gamma <- function(x){paste0("Gamma(", x[1], ", ", x[2], ")")}

vec_hyperparm <- unlist(map(list_hyperparm, get_character_Gamma))  


### The number of iteration per each option
N_iter <- 30


### Prepare an empty object to save the results
array_result <- array(data = NA, dim = c(N_iter, 14, length(list_hyperparm)))

temp <- array(data = NA, dim = c(nsave/2, N_iter, length(list_hyperparm)))
array_posterior <- array_prior <- array_K_posterior <- array_weights <- temp 



###'######################################################################
###'
###' Implement the loop over
###' 
###' j: hyperparameter options - 5 choices
###' i: iteration - 30 iterations
###'
###'

# for (j in seq_along(list_hyperparm)){
for (j in 2:5){
  
  for (i in seq(N_iter)){
    
    ### Extract loop elements 
    a <- list_hyperparm[[j]][1]
    b <- list_hyperparm[[j]][2]
    array_result[i, 1, j] <- a
    array_result[i, 2, j] <- b
    
    
    ###' Calculate the theoretical mean and variance of 
    ###' the assumed Gamma distribution and Number of clusters
    array_result[i, 3, j] <- a/b
    array_result[i, 4, j] <- a/b^2
    
    
    ###'######################################################################
    ###'
    ###' (1) Obtain marginal posterior density of alpha 
    ###'     given data, and parameters for the Gamma distribution
    ###'     
    ###'     \pi(\alpha | data, a, b)
    ###'
    ###'
    
    ### Set the prior
    prior <- list(a0 = a,     # alpha0: shape param.
                  b0 = b,     # alpha0: rate param.  
                  tau1 = 1,   # G0 variance: shape param.
                  tau2 = 1,   # G0 variance: rate param. 
                  mub = 0, 
                  Sb = 100)   # G0 mean: mean 
    
    
    ### Estimate the Dirichlet Process model using DPpackage 
    outp <- DPmeta(formula = mat_DPmeta ~ 1,
                   prior = prior, mcmc = mcmc, 
                   state = state, status = TRUE)
    
    
    ### Tidy up and save the DPmeta object
    df_posterior <- get_posterior_DPmeta(outp, nburn)
    
    
    ### Get alpha and K given data and parameters for the Gamma distribution
    array_posterior[ , i, j] <- alpha_posterior <- df_posterior$alpha0
    array_K_posterior[ , i, j] <- K_posterior <- df_posterior$Ncluster

    
    ### Calculate the means and variances of alpha and K
    array_result[i, 5, j] <- alpha_posterior_mean <- mean(alpha_posterior)
    array_result[i, 6, j] <- alpha_posterior_var <- var(alpha_posterior)
    array_result[i, 7, j] <- K_posterior_mean <- mean(K_posterior)
    array_result[i, 8, j] <- K_posterior_var <- var(K_posterior)
    
    
    ###'######################################################################
    ###'
    ###' (2) Get weights
    ###'
    ###'
    
    ### 1. Denominator: alpha prior with the specific a, b parameters
    # set.seed(12345)
    N_smp <- length(alpha_posterior)
    i_vec <- seq(N_smp)
    
    array_prior[ , i, j] <- alpha_prior <- rgamma(n = N_smp, shape = a, rate = b)
    
    array_result[i, 9, j] <- alpha_prior_mean <- mean(alpha_prior)
    array_result[i, 10, j] <- alpha_prior_var <- var(alpha_prior)
    
    array_result[i, 11, j] <- K_prior_mean <- 
      sum(alpha_prior/(alpha_prior + i_vec - 1))
    array_result[i, 12, j] <- K_prior_var <-  
      sum((alpha_prior*(i_vec - 1))/(alpha_prior + i_vec -1)^2)
    
    
    ###' 2. Numerator: a uniform prior for \alpha on the interval (0, C)
    ###'    with C > N to ensure we assign mass to all resonable values of \alpha
    share_fct <- 0.01
    C <- N + (share_fct*N)
    alpha_uniform <- runif(n = N_smp, min = 0, max = C)
    
    
    ### Calculate and save the weight
    alpha_weight <- alpha_uniform/alpha_prior
    array_weights[ , i, j] <- alpha_weight
    
    
    ###'######################################################################
    ###'
    ###' (3) Importance sampling 
    ###' => The MLE of alpha is very sensitive to the a, b parameter settting. 
    ###'
    ###'
    
    ### Calculate and save the MLE for alpha
    array_result[i, 13, j] <- diagis::weighted_mean(alpha_posterior, alpha_weight)
    array_result[i, 14, j] <- diagis::weighted_var(alpha_posterior, alpha_weight)

    
    ###'######################################################################
    ###'
    ###' Print progress
    ###'
    ###'
    
    cat(paste0("j = ", j, ", i = ", i, "\n"))
  }
}



###'######################################################################
###'
###' Tidy up the simulation results
###'
###' (1) array_result
###' 
###' 

### Convert the array into a dataframe 
list_temp <- array_branch(array_result, 3)
list_result <- map(list_temp, data.frame)
df_result <- bind_rows(list_result)


### Assign variable names & generate a factor variable 
names(df_result) <- c("a", "b", "A_mean_theory", "A_var_theory", 
                      "A_mean_poster", "A_var_poster", "K_mean_poster", "K_var_poster", 
                      "A_mean_prior", "A_var_prior", "K_mean_prior", "K_var_prior", 
                      "A_mean_IS", "A_var_IS")

df_result <- df_result %>%
  mutate(hyperparm = paste0("Gamma(", a, ", ", b, ")")) %>%
  mutate(hyperparm = factor(hyperparm, levels = vec_hyperparm)) %>%
  dplyr::select(hyperparm, everything()) 


### Save the resulting dataframe 
saveRDS(df_result, file = "datasets/Importance sampling simulation.rds")



###'######################################################################
###'
###' Tidy up the simulation results
###'
###' (2) array_weights
###' (3) array_priors
###' 
###' 

### Convert the array into a dataframe 
list_temp <- array_branch(array_weights, 3)
list_weights <- map(list_temp, data.frame)
names(list_weights) <- vec_hyperparm
saveRDS(list_weights, file = "datasets/Importance sampling weights.rds")


### Convert the array into a dataframe 
list_temp <- array_branch(array_prior, 3)
list_priors <- map(list_temp, data.frame)
names(list_priors) <- vec_hyperparm
saveRDS(list_priors, file = "datasets/Importance sampling priors.rds")



###'######################################################################
###'
###' Generate a table and a plot summarizing the simulation result 
###'
###'

### Table
df_tbl <- df_result %>%
  gather(key = key, value = value, -hyperparm, factor_key = TRUE) %>%
  group_by(hyperparm, key) %>%
  summarise(mean_value = mean(value), 
            var_value = var(value))

write.csv(df_tbl, file = "tables/Importance sampling result.csv")


### Plot 01: Alpha means and variances
df_plot <- df_result %>%
  gather(key = key, value = value, -hyperparm, factor_key = TRUE) %>%
  separate(key, sep = "_", c("parm", "stat", "type")) %>%
  mutate(type = factor(type, levels = c("theory", "prior", "poster", "IS"), 
                       labels = c("Theory","Prior", "Posterior", "IS")), 
         stat = factor(stat, levels = c("mean", "var"), labels = c("Mean", "Variance")))


p <- ggplot(data = df_plot %>% filter(parm == "A"), 
            aes(x = type, y = value)) + 
  geom_boxplot(aes(colour = type)) + 
  facet_wrap(hyperparm ~ stat, scales = "free") + 
  theme_preset + theme(legend.position = "none") + 
  labs(title = "Means and variances of the precision parameter", 
       subtitle = "The number of iteration = 30", 
       x = "Samples", y = "Value")

ggsave(filename = "figures/Means and variances of the precision parameter.pdf", 
       p, width = 11, height = 8)

