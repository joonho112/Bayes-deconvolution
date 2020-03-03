
###'######################################################################
###'
###' Category: Code Validation
###' 
###' Task: Compare packages to implement Dirichlet Process Mixtures
###'       
###' Data: Simulated data
###' 
###' Data: 2020-03-02
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


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)


### Set ggplot themes
### Theme settings
theme_preset <- 
  theme_bw() + 
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

### Generate sample draws from the mixture model
set.seed(12345)

df_G <- GenerateMdp(n = 50, rgt = 1, rsr = 5)[c(1:3)] %>%
  data.frame() %>%
  rownames_to_column("K") %>%
  mutate_at(.vars = c("K"), .funs = as.numeric)

head(df_G)


### Save the resulting sample
setwd(work_dir)
write.csv(df_G, file = "datasets/df_normal_mixtures.csv")



###'######################################################################
###'
###' Parameter estimation using the `hhsim` package 
###' 
###' This packages is very unstable when being used in RStudio
###' Recommend excute this within original R console
###'
###'

### Import the generated dataset from a normal mixture distribution
df_G <- read.csv(file = "datasets/df_normal_mixtures.csv")


### Set numbers of iterations
nmcmc <- 2000    # number of MCMC iterations per each model fit
nburn <- 2000    # burn-in iterations to discard at each round
nmcmcall <- nmcmc + nburn  


### Extract vectors of observed tau_ks and their SEs
Y <- df_G$tau_k_hat   # a vector of observed tau_ks
sigma2 <- df_G$se_k2  # a vector of known SEs
nDraws <- nmcmcall    # number of MCMC draws


### Estimate the Dirichlet Process model using hhsim package
### This results in aborted R sesssion. Try with the saved outputs
# outp <- GaussianMDP(Y, sigma2, nDraws,
#                     numconfig = 3,  # N of clusters
#                     m = 0,          # mean of G0     
#                     w = 1,          # G0 variance: shape param. 
#                     W = 1,          # G0 variance: rate param. 
#                     par1 = 4,       # alpha0: shape param.
#                     par2 = 4)       # alpha0: rate param.

outp_hhsim <- readRDS("datasets/outp.RDS") 


### Define a function to tidy up the output object
get_posterior_hhsim <- function(output = outp_hhsim, nburn = 100, nmcmcall = 600){
  
  # (1) Tidy up the "theta" output object
  df_theta <- outp[["theta"]] %>%
    unlist() %>%
    matrix(nrow = length(outp[["theta"]]), byrow = TRUE) %>%
    data.frame() %>%
    slice(c((nburn + 1):nmcmcall))  # throw away burn-ins
  colnames(df_theta) <- paste0("tau_k[", seq(ncol(df_theta)), "]")
  
  # (2) Tidy up the G0 parameters (m, s2), alpha0, and N of clusters outputs
  df_hyperparm <- data.frame(outp[["m"]], outp[["s2"]], 
                             outp[["alpha"]], outp[["numconfig"]]) %>%
    slice(c((nburn + 1):nmcmcall))  # throw away burn-ins
  colnames(df_hyperparm) <- c("G0_mu", "G0_s2", "alpha0", "Ncluster")
  
  ### Return the resulting object
  cbind.data.frame(df_hyperparm, df_theta)
}


### Save posterior samples
posterior_hhsim <- get_posterior_hhsim(outp_hhsim, nburn, nmcmcall) %>%
  mutate(Type = "hhsim") %>% 
  dplyr::select(Type, everything())



###'######################################################################
###'
###' Parameter estimation using the `DPpackage` package
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


### Set prior parameters
prior <- list(a0 = 4,     # alpha0: shape param.
              b0 = 4,     # alpha0: rate param.  
              tau1 = 1,   # G0 variance: shape param.
              tau2 = 1,   # G0 variance: rate param. 
              mub = 0, 
              Sb = 100)  # G0 mean: mean 


### Prepare dataset as a matrix form (with y and sigma2)
mat_DPmeta <- df_G %>% select(tau_k_hat, se_k2) %>% as.matrix()


### Estimate the Dirichlet Process model using DPpackage 
outp_DPmeta <- DPmeta(formula = mat_DPmeta ~ 1,
                      prior = prior, mcmc = mcmc, 
                      state = state, status = TRUE)


### Summary with HPD and Credibility intervals
summary(outp_DPmeta)
summary(outp_DPmeta, hpd = FALSE)


### Plot model parameters (to see the plots gradually set ask=TRUE)
plot(outp_DPmeta)
plot(fit, ask = TRUE, nfigr = 2, nfigc = 2)


### Define a function to tidy up the DPmeta objects
get_posterior_DPmeta <- function(output = outp_DPmeta, nburn = 4000){

  # (1) Tidy up the "theta" output object
  df_theta <- outp_DPmeta$save.state$randsave %>%
    data.frame() %>% 
    dplyr::select(-Prediction) %>%
    slice(c((nburn/2 + 1):nburn))  # throw away burn-ins
  colnames(df_theta) <- paste0("tau_k[", seq(ncol(df_theta)), "]")
  
  # (2) Tidy up the G0 parameters (m, s2), alpha0, and N of clusters outputs
  df_hyperparm <- outp_DPmeta$save.state$thetasave %>%
    data.frame() %>%
    dplyr::select(-tau_k_hat) %>%
    rename(G0_mu = mu, G0_s2 = sigma2, Ncluster = ncluster, alpha0 = alpha) %>%
    dplyr::select(G0_mu, G0_s2, alpha0, Ncluster) %>%
    slice(c((nburn/2 + 1):nburn))  # throw away burn-ins
  
  # Return the resulting object
  cbind.data.frame(df_hyperparm, df_theta)
}

### Save posterior samples
posterior_DPmeta <- get_posterior_DPmeta(outp_DPmeta, nburn) %>%
  mutate(Type = "DPmeta") %>% 
  dplyr::select(Type, everything())



###'######################################################################
###'
###' Parameter estimation using the `bspmma` package
###'
###'

### Check up the precision parameter estimates from hhsim and DPmeta
summary(posterior_hhsim$alpha0)   # Mean = 1.2661
summary(posterior_DPmeta$alpha0)  # Mean = 1.0471


### Prepare dataset as a matrix form (with y and sigma2)
mat_bspmma <- df_G %>% select(tau_k_hat, se_k2) %>% as.matrix()


### Estimate the Dirichlet Process model using bspmma 
ncycles <- nburn <- 4000   # Number of cycles of the Markov chain
alpha0 <- 1                # The precision parameter of the DP prior

outp_bspmma <- dirichlet.c(mat_bspmma,  
                           ncycles = nburn,  
                           M = alpha0)           


### Define a function to tidy up the bspmma objects
get_posterior_bspmma <- function(output = outp_DPmeta, nburn = 4000){

  # (1) Tidy up the "theta" output object
  df_theta <- outp_bspmma$chain %>%
    data.frame() %>%
    slice(c((nburn/2 + 1):nburn)) %>% # throw away burn-ins
    dplyr::select(-mu, -tau) 
  colnames(df_theta) <- paste0("tau_k[", seq(ncol(df_theta)), "]")
  
  # (2) Tidy up the G0 parameters (m, s2)
  df_hyperparm <- outp_bspmma$chain %>%
    data.frame() %>%
    slice(c((nburn/2 + 1):nburn)) %>% # throw away burn-ins
    dplyr::select(mu, tau) %>%
    rename(G0_mu = mu, G0_s2 = tau)
  
  # Return the resulting object
  cbind.data.frame(df_hyperparm, df_theta)
}


### Save posterior samples
posterior_bspmma <- get_posterior_bspmma(outp_bspmma, nburn) %>%
  mutate(Type = "bspmma") %>% 
  dplyr::select(Type, everything())



###'######################################################################
###'
###' Comparison of posterior samples obtained from different packages
###'
###' : `hhsim`, `DPmeta`, and `bspmma`
###' 
###'

### Combine dataframes
list_collect <- list(posterior_hhsim, 
                     posterior_DPmeta, 
                     posterior_bspmma)

df_combine <- bind_rows(list_collect) %>%
  mutate(Type = factor(Type, levels = c("hhsim", "DPmeta", "bspmma")))
  


###'######################################################################
###'
###' (1) Hyperparameters $\tau$ and $\sigma_{\tau}$
###'
###'

### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(Type, G0_mu, G0_s2) %>%
  mutate(G0_s2 = log(G0_s2)) %>%    # Log-transformation for G0_s2
  gather(key = variable, value = value, -Type) %>%
  mutate(variable = factor(variable, levels = c("G0_mu", "G0_s2")))


### Define a function to generate plot
plot_compare_density <- function(df_plot, title = title){
  ggplot(data = df_plot, aes(x = value, group = Type, fill = Type)) +
    geom_density(position = "identity", size = 0.1, alpha = 0.3) + 
    geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
    labs(title = title) + theme_preset
}


### Plot!
title <- c("Posterior densities of hyperparameters")
plot_compare_density(df_plot, title = title) + 
  facet_wrap(~ variable, scales = "free")



###'######################################################################
###'
###' (2) Precision parameter (alpha0)
###'
###'

### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(Type, alpha0) %>%
  gather(key = variable, value = value, -Type) %>%
  mutate(variable = factor(variable, levels = c("alpha0")))


### Plot!
title <- c("Posterior densities of precision parameter (alpha0)")
plot_compare_density(df_plot, title = title) + 
  facet_wrap(~ variable, scales = "free")



###'######################################################################
###'
###' (3) Number of clusters
###'
###'

### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(Type, Ncluster) %>%
  gather(key = variable, value = value, -Type) 


### Define a function to generate plot
plot_compare_histogram <- function(df_plot, title = title){
  ggplot(data = df_plot, aes(x = value, group = Type, fill = Type)) +
    geom_histogram(position = "dodge", size = 1, alpha = 0.3) + 
    geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
    labs(title = title) + theme_preset
}


### Plot!
title <- c("Posterior densities of the number of clusters")
plot_compare_histogram(df_plot, title = title) + 
  facet_wrap(~ variable, scales = "free") + 
  scale_x_continuous(breaks = 1:max(df_plot$value, na.rm = TRUE))



###'######################################################################
###'
###' Empirical distribution function (EDF) estimates for $\tau_{k}$s 
###'
###'

### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(Type, contains("tau_k[")) %>%
  gather(key = variable, value = value, -Type) %>%
  group_by(Type, variable) %>%
  summarise(mean = mean(value, na.rm = TRUE), 
            sd = sd(value, na.rm = TRUE)) %>%
  mutate(id = as.numeric(str_extract(variable, "\\d+"))) %>%
  arrange(id, Type)

### Plot!
p1 <- plot_compare_density(df_plot %>% rename(value = mean), 
                           title = "EDF of posterior means")

p2 <- plot_compare_density(df_plot %>% rename(value = sd), 
                           title = "EDF of posterior SDs")

plot_grid(p1, p2, labels = "AUTO")



###'######################################################################
###'
###' Site-specific individual effect estimates $\tau_{k}$s
###'
###'

### Randomly sample 10 sites
set.seed(12345)
rand_var <- paste0("tau_k[", sample(1:50, 10, replace = FALSE), "]")


### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(Type, rand_var) %>%
  gather(key = variable, value = value, -Type)


### Plot
title <- "Posterior densities of site-specific estimates"
plot_compare_density(df_plot, title = title) + 
  facet_wrap(~ variable, scales = "free")



