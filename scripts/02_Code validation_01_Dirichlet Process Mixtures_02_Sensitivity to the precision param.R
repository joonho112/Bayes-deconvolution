
###'######################################################################
###'
###' Category: 
###' 
###' Task: (1) Dirichlet Process Mixtures
###'           
###'        Experiments for DPM (following Sophia's guide)
###'        
###'        => Compare outcomes with different settings of the precision parameter 
###'           
###' Data: Simulated data
###' 
###' Data: 2020-03-10
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
###' Parameter estimation using the `DPpackage` package
###'
###' (1) Fix alpha to a small, mediaum, and large value
###'     and compare them to the base distribution (Normal)
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


### Prepare a loop over fixed alpha
alpha_vec <- c(0.1, 1, 100)
list_collect <- list()


### Start a loop over alpha_vec

for (i in seq_along(alpha_vec)){
  
  ## Extract an alpha element
  alpha_value <- alpha_vec[i]
  
  
  ## Set prior parameters
  prior <- list(alpha = alpha_value,  # fixed precision parameter
                tau1 = 1,   # G0 variance: shape param.
                tau2 = 1,   # G0 variance: rate param. 
                mub = 0, 
                Sb = 100)   # G0 mean: mean 
  
  
  ## Estimate the Dirichlet Process model using DPpackage 
  outp_DPmeta <- DPmeta(formula = mat_DPmeta ~ 1,
                        prior = prior, mcmc = mcmc, 
                        state = state, status = TRUE)
  
  
  ## Save as a list element
  list_collect[[i]] <- outp_DPmeta
  
}  ### End of a loop over alpha_vec



###'######################################################################
###'
###' Parameter estimation using the `DPpackage` package
###' 
###' (2) Collect posterior samples with different settings of 
###' 
###'     the precision parameter alpha
###'
###'

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


### Collect posterior samples per each alpha setting

list_df <- list()

for (i in seq_along(alpha_vec)){
  list_df[[i]] <- get_posterior_DPmeta(list_collect[[i]], nburn)
}

df_collect <- bind_rows(list_df) %>%
  mutate(alpha0 = as.character(alpha0))

classmode(df_collect, everything())
tabdf(df_collect, alpha0)



###'######################################################################
###'
###' Parameter estimation using **Stan**
###'
###'

### Collect data into a list format suitable for Stan
stan_data <- list(K = length(df_G$K), 
                  tau_hat_k = df_G$tau_k_hat, 
                  se_k = sqrt(df_G$se_k2))


### Compile and run the stan model
fit_norm <- stan(file = "scripts/rubin_model_Gaussian_Invgamma prior.stan",
                 data = stan_data,
                 iter = 1000, chains = 4)

print(fit_norm, probs = c(0.1, 0.5, 0.9), digits = 3)


### Extract posterior samples
df_norm <- as.data.frame(fit_norm) %>%
  mutate(alpha0 = "Gaussian") %>%
  rename(G0_mu = tau, G0_s2 = sigma_tau2) %>%
  dplyr::select(-sigma_tau, -predicted_tau_k, -lp__) %>%
  dplyr::select(G0_mu, G0_s2, alpha0, everything())
  

### Append to the results from DPpackage
dim(df_collect); names(df_collect)
dim(df_norm); names(df_norm)

alpha_vec

df_combine <- bind_rows(list(df_norm, df_collect)) %>% 
  mutate(alpha0 = factor(alpha0, levels = c("Gaussian", alpha_vec)))

tabdf(df_combine, alpha0)



###'######################################################################
###'
###'  Compare the outcomes from different settings of the alpha  
###' 
###' (1) Hyperparameters $\tau$ and $\sigma_{\tau}$
###'
###'

### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(alpha0, G0_mu, G0_s2) %>%
  mutate(G0_s2 = log(G0_s2)) %>%    # Log-transformation for G0_s2
  gather(key = variable, value = value, -alpha0) %>%
  mutate(variable = factor(variable, levels = c("G0_mu", "G0_s2")))


### Define a function to generate plot
plot_compare_density <- function(df_plot, title = title, subtitle = subtitle){
  ggplot(data = df_plot, aes(x = value, group = alpha0, fill = alpha0)) +
    geom_density(position = "identity", size = 0.1, alpha = 0.3) + 
    geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
    labs(title = title, subtitle = subtitle) + theme_preset
}


### Plot!
title <- c("Posterior densities of hyperparameters")
subtitle <- c("by different settings of the fixed precision parameter")
plot_compare_density(df_plot, title = title, subtitle = subtitle) + 
  facet_wrap(~ variable, scales = "free")



###'######################################################################
###'
###' (2) Number of clusters
###'
###'

### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(alpha0, Ncluster) %>%
  gather(key = variable, value = value, -alpha0) %>%
  drop_na()


### Define a function to generate plot
plot_compare_histogram <- function(df_plot, title = title, subtitle = subtitle){
  ggplot(data = df_plot, aes(x = value, group = alpha0, fill = alpha0)) +
    geom_histogram(position = "dodge", size = 1, alpha = 0.3) + 
    geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
    labs(title = title, subtitle = subtitle) + theme_preset
}


### Plot!
title <- c("Posterior densities of the (assumed) number of clusters")
plot_compare_histogram(df_plot, title = title, subtitle = subtitle) + 
  facet_wrap(~ variable, scales = "free") + 
  scale_x_continuous(breaks = 1:max(df_plot$value, na.rm = TRUE))



###'######################################################################
###'
###' Empirical distribution function (EDF) estimates for $\tau_{k}$s 
###'
###'

### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(alpha0, contains("tau_k[")) %>%
  gather(key = variable, value = value, -alpha0) %>%
  group_by(alpha0, variable) %>%
  summarise(mean = mean(value, na.rm = TRUE), 
            sd = sd(value, na.rm = TRUE)) %>%
  mutate(id = as.numeric(str_extract(variable, "\\d+"))) %>%
  arrange(id, alpha0) 


### Plot!
p1 <- plot_compare_density(df_plot %>% rename(value = mean), 
                           title = "EDF of posterior means", 
                           subtitle = NULL)

p2 <- plot_compare_density(df_plot %>% rename(value = sd), 
                           title = "EDF of posterior SDs", 
                           subtitle = NULL)

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
  dplyr::select(alpha0, rand_var) %>%
  gather(key = variable, value = value, -alpha0) 


### Plot
title <- "Posterior densities of site-specific estimates"
plot_compare_density(df_plot, title = title, subtitle = subtitle) + 
  facet_wrap(~ variable, scales = "free")


