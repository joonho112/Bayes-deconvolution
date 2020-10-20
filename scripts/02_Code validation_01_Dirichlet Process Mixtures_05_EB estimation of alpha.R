
###'######################################################################
###'
###' Category: 
###' 
###' Task: (1) Dirichlet Process Mixtures
###'           
###'        Selecting a Prior for precision parameter alpha
###'        
###'        2. Empirical Bayes estimation of alpha
###'           
###' Data: Simulated data
###' 
###' Data: 2020-03-22
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

N <- 50
rgt <- 0.33
rsr <- 5

df_G <- GenerateMdp(n = N, rgt = rgt, rsr = rsr)[c(1:3)] %>%
  data.frame() %>%
  rownames_to_column("K") %>%
  mutate_at(.vars = c("K"), .funs = as.numeric)



###'######################################################################
###'
###' Obtain the Empirical Bayes estimation of the precision parameter alpha
###'
###'

###'######################################################################
###'
###' (1) Set the initial value of alpha
###' 
###' A guess within the range 1/log(n) to n/log(n)
###' (McAuliffe et al., 2006) 
###' 
###' 

init_range <- c(1/log(N), N/log(N))

init_alpha <- (init_range[1] + init_range[2])/2



###'######################################################################
###'
###' (2) Use MCMC to generate posterior smaples given the fixed alpha
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


### Prepare dataset as a matrix form (with y and sigma2)
mat_DPmeta <- df_G %>% select(tau_k_hat, se_k2) %>% as.matrix()


## Set prior parameters
prior <- list(alpha = init_alpha,  # fixed precision parameter
              tau1 = 1,   # G0 variance: shape param.
              tau2 = 1,   # G0 variance: rate param. 
              mub = 0, 
              Sb = 100)   # G0 mean: mean 


## Estimate the Dirichlet Process model using DPpackage 
outp <- DPmeta(formula = mat_DPmeta ~ 1,
               prior = prior, mcmc = mcmc, 
               state = state, status = TRUE)


### Define a function to tidy up the DPmeta objects
get_posterior_DPmeta <- function(outp = outp_DPmeta, nburn = 4000){
  
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


### Get the posterior samples
df_posterior <- get_posterior_DPmeta(outp)



###'######################################################################
###'
###' (3) Approximating the posterior mean of the number of clusters (K)
###'
###'

K_postmean <- df_posterior$Ncluster %>%
  mean(na.rm = TRUE)



###'######################################################################
###'
###' (4) Compute the value of alpha that satisfies 
###' 
###'     Kbar = sum_{i=1}^{N}(alpha/(alpha + i - 1))
###'
###'

K_postmean

i_vec <- seq(from = 1, to = N, by = 1)

sum(alpha/(alpha + seq(from = 1, to = N, by = 1) - 1))


fn_K <- function(alpha, N, K_postmean){
  
  seq(from = 1, to = N, by = 1)
  abs(K_postmean - sum(alpha/(alpha + i_vec - 1))) 

}

solution <- optimize(fn_K, c(0, 100), N = 50, K_postmean, tol = 0.0000001)

new_alpha <- solution$minimum



###'######################################################################
###'
###' Construct a `repeat` loop for the alpha convergence
###'
###'

### Set the initial value of alpha and tolerance level
init_range <- c(1/log(N), N/log(N))
alpha0 <- (init_range[1] + init_range[2])/2
# tol <- 0.01
tol <- 1e-3
i <- 1   # Loop indicator
vec_alpha <- vector()   # allocate the space to store alpha values


### Start a repeat loop
repeat {
  
  ## Set prior parameters
  prior <- list(alpha = alpha0,  # fixed precision parameter
                tau1 = 1,   # G0 variance: shape param.
                tau2 = 1,   # G0 variance: rate param. 
                mub = 0, 
                Sb = 100)   # G0 mean: mean 
  
  
  ## Estimate the Dirichlet Process model using DPpackage 
  outp <- DPmeta(formula = mat_DPmeta ~ 1,
                 prior = prior, mcmc = mcmc, 
                 state = state, status = TRUE)
  
  
  ### Get the posterior samples
  df_posterior <- get_posterior_DPmeta(outp)
  
  
  ###' Compute the value of alpha that satisfies 
  ###' Kbar = sum_{i=1}^{N}(alpha/(alpha + i - 1))
  K_postmean <- df_posterior$Ncluster %>%
    mean(na.rm = TRUE)
  
  fn_K <- function(alpha, N, K_postmean){
    
    seq(from = 1, to = N, by = 1)
    abs(K_postmean - sum(alpha/(alpha + i_vec - 1))) 
    
  }
  
  solution <- optimize(fn_K, c(0, 100), N = 50, K_postmean, tol = 0.0000001)
  
  alpha_vec[i] <- alpha1 <- solution$minimum
  
  
  ### Test convergence
  if(abs(alpha1 - alpha0) < tol) {  ## Close enough?
    break
  } else {
    
    # Assign new values
    i <- i + 1
    alpha0 <- alpha1
    
    # Print progress
    cat(paste0("i = ", i, ", alpha = ", alpha0), "\n")
  } 
}

### Tidy up and plot the alpha vector
df_alpha <- data.frame(alpha_vec) %>%
  rownames_to_column() %>% 
  rename(iteration = rowname, alpha = alpha_vec) %>%
  mutate(iteration = as.numeric(iteration))

p <- plot_trend_xy(df_alpha, iteration, alpha, sprintf = "%#.2f") + 
  labs(title = "Empirical Bayes estimates of the precision parameter (alpha)", 
       subtitle = "N = 50", 
       x = "Iteration", 
       y = "Precision parameter (alpha) estimate", 
       caption = NULL)

p

ggsave(filename = paste0("figures/EB estimates of alpha_N", N, ".pdf"), p, 
       width = 10, height = 6)

