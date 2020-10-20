
###'######################################################################
###'
###' Category: 
###' 
###' Task: (1) Dirichlet Process Mixtures
###'           
###'        Experiments for DPM (following Sophia's guide)
###'        
###'        => Compare alpha (precision parameter)'s prior and posterior
###'           Does the data provide information on alpha? 
###'           
###' Data: Simulated data
###' 
###' Data: 2020-03-11
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
###' Construct a loop over the informativeness of data
###' 
###' varying `rgt`: geometric mean of uncertainty estimates
###' 
###' Large value indicates relatively less information in the data about \theta's
###'
###'

### Define a vector of varying values of geometric mean 
rgt_vec <- c(0.10, 0.33, 1.00)


###' Define a list of multiple parameter choices for alpha
###' (1) alpha0 ~ Gamma(4, 4): bumpy and multimodal G (5 expected clusters)
###' (2) alpha0 ~ Gamma(10, 0.1): smoother G (70 expected clusters)
ab_list <- list(c(4, 4), c(10, 0.1)) 


### Generate a datafrmae for the loop reference
loop_ref <- expand.grid(seq_along(rgt_vec), seq_along(ab_list))
list_collect <- list()

for (i in seq(nrow(loop_ref))){
  
  ## Extract loop elements 
  rgt <- rgt_vec[loop_ref[i, ][[1]]]
  alpha_param <- ab_list[[loop_ref[i, ][[2]]]]
  
  ### Generate sample draws from the mixture model
  set.seed(12345)
  
  df_G <- GenerateMdp(n = 50, rgt = rgt, rsr = 5)[c(1:3)] %>%
    data.frame() %>%
    rownames_to_column("K") %>%
    mutate_at(.vars = c("K"), .funs = as.numeric)
  
  
  ### Set MCMC parameters
  state <- NULL    # the current value of the parameters
  nburn <- 4000    # the number of burn-in scans
  nsave <- 4000    # the total number of scans to be saved
  nskip <- 20      # the thinning interval
  ndisplay <- 100  # the number of saved scans to be displayed on screen
  
  mcmc <- list(nburn = nburn, nsave = nsave,
               nskip = nskip, ndisplay = ndisplay)
  
  
  ### Set prior parameters
  prior <- list(a0 = alpha_param[1],     # alpha0: shape param.
                b0 = alpha_param[2],     # alpha0: rate param.  
                tau1 = 1,   # G0 variance: shape param.
                tau2 = 1,   # G0 variance: rate param. 
                mub = 0, 
                Sb = 100)  # G0 mean: mean 
  
  
  ### Prepare dataset as a matrix form (with y and sigma2)
  mat_DPmeta <- df_G %>% select(tau_k_hat, se_k2) %>% as.matrix()
  
  
  ### Estimate the Dirichlet Process model using DPpackage 
  outp <- DPmeta(formula = mat_DPmeta ~ 1,
                 prior = prior, mcmc = mcmc, 
                 state = state, status = TRUE)
  
  ### Save the output
  list_collect[[i]] <- outp
  
} ### End of loop 



###'######################################################################
###'
###' Collect posterior samples 
###'
###'

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


### Collect posterior samples per each simulation setting

meta_list <- list()

for (i in seq(nrow(loop_ref))){
  
  ## Extract loop elements 
  rgt <- rgt_vec[loop_ref[i, ][[1]]]
  alpha_param <- ab_list[[loop_ref[i, ][[2]]]]
  
  
  ##  Get posterior samples
  meta_list[[i]] <- get_posterior_DPmeta(list_collect[[i]], nburn) %>%
    mutate(rgt = paste0("rgt = ", sprintf("%.2f", rgt)), 
           alpha_prior = paste0("gamma(", alpha_param[1], ", ", 
                                alpha_param[2], ")")) %>%
    dplyr::select(rgt, alpha_prior, everything())
  
}

df_combine <- bind_rows(meta_list) %>% 
  mutate(rgt = factor(rgt, levels =  paste0("rgt = ", sprintf("%.2f", rgt_vec))), 
         alpha_prior = factor(alpha_prior, 
                              levels = c("gamma(4, 4)", "gamma(10, 0.1)")))
  
dim(df_combine)



###'######################################################################
###'
###' Get alpha priors - It's simple
###' 
###' => draw the same number of samples as posterior (n = 2,000)
###'
###'

### (1) alpha ~ gamma(4, 4)
set.seed(12345)
gm_prior1 <- rgamma(n = 2000, shape = 4, rate = 4)
plot(density(gm_prior1))
mean(gm_prior1); sd(gm_prior1)


### (2) alpha ~ gamma(10, 0.1)
set.seed(12345)
gm_prior2 <- rgamma(n = 2000, shape = 10, rate = 0.1)
plot(density(gm_prior2))
mean(gm_prior2); sd(gm_prior2)



###'######################################################################
###'
###' Compare densities  
###' 
###' (1) the precision parameter (alpha) Prior vs. posterior
###'
###' => six plots by condition
###'
###'

###' Generate a dataframe to plot
###' Bind gamma prior samples
Prior <- c(rep(gm_prior1, times = 3), rep(gm_prior2, times = 3))
length(Prior)  
  
df_plot <- df_combine %>%
  dplyr::select(rgt, alpha_prior, alpha0) %>%
  arrange(alpha_prior, rgt) %>%
  rename(Posterior = alpha0) %>%
  cbind.data.frame(Prior) %>%
  gather(key = sample, value = value, Prior, Posterior) %>%
  mutate(sample = factor(sample, levels = c("Prior", "Posterior")))


### Plot
title <- c("Precision parameter: Prior vs. Posterior")
subtitle <- c("by the informativeness of the data (rgt) & alpha prior")
ggplot(data = df_plot, aes(x = value, group = sample, fill = sample)) +
  geom_density(position = "identity", size = 0.1, alpha = 0.3) + 
  geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
  labs(title = title, subtitle = subtitle, x = NULL) + theme_preset +
  facet_wrap(alpha_prior ~ rgt, scales = "free") + 
  theme_minimal()



###'######################################################################
###'
###'  Compare the outcomes from different settings of 
###'  
###'  The informativeness of the data (rgt) & Prior for the precision param.  
###' 
###' 
###' (1) Hyperparameters $\tau$ and $\sigma_{\tau}$
###'
###'

### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(rgt, alpha_prior, G0_mu, G0_s2) %>%
  mutate(G0_s2 = log(G0_s2)) %>%    # Log-transformation for G0_s2
  gather(key = variable, value = value, G0_mu, G0_s2) %>%
  mutate(variable = factor(variable, levels = c("G0_mu", "G0_s2")))


### Plot!
title <- c("Posterior densities of hyperparameters")
subtitle <- c("by the informativeness of the data (rgt) & alpha prior")

ggplot(data = df_plot, aes(x = value, group = alpha_prior, fill = alpha_prior)) +
  geom_density(position = "identity", size = 0.1, alpha = 0.3) + 
  geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
  labs(title = title, subtitle = subtitle, x = NULL) + theme_preset +
  facet_wrap(variable ~ rgt, scales = "free") + 
  theme_minimal()



###'######################################################################
###'
###' (2) Number of clusters
###'
###'

### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(rgt, alpha_prior, Ncluster) 


### Plot!
title <- c("Posterior densities of the (assumed) number of clusters")
ggplot(data = df_plot, aes(x = Ncluster, group = rgt, fill = rgt)) +
  geom_histogram(position = "dodge", size = 1, alpha = 0.3) + 
  geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
  labs(title = title, subtitle = subtitle) + theme_preset + 
  facet_wrap(~ alpha_prior, scales = "free") + 
  scale_x_continuous(breaks = min(df_plot$Ncluster):max(df_plot$Ncluster))



###'######################################################################
###'
###' Empirical distribution function (EDF) estimates for $\tau_{k}$s 
###'
###'

### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(rgt, alpha_prior, contains("tau_k[")) %>%
  gather(key = variable, value = value, -rgt, -alpha_prior) %>%
  group_by(rgt, alpha_prior, variable) %>%
  summarise(mean = mean(value, na.rm = TRUE), 
            sd = sd(value, na.rm = TRUE)) %>%
  mutate(id = as.numeric(str_extract(variable, "\\d+"))) %>%
  arrange(rgt, alpha_prior, id) 


### Plot!
p1 <- ggplot(data = df_plot %>% rename(value = mean), 
             aes(x = value, group = alpha_prior, fill = alpha_prior)) +
  geom_density(position = "identity", size = 0.1, alpha = 0.3) + 
  geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
  labs(title = "EDF of posterior means", subtitle = subtitle, x = NULL) + theme_preset +
  facet_wrap(. ~ rgt, scales = "free") + 
  theme_minimal()


p2 <- ggplot(data = df_plot %>% rename(value = sd), 
             aes(x = value, group = alpha_prior, fill = alpha_prior)) +
  geom_density(position = "identity", size = 0.1, alpha = 0.3) + 
  geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
  labs(title = "EDF of posterior SDs", subtitle = subtitle, x = NULL) + theme_preset +
  facet_wrap(. ~ rgt, scales = "free") + 
  theme_minimal()

plot_grid(p1, p2, nrow = 2, labels = "AUTO")



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
  dplyr::select(rgt, alpha_prior, rand_var) %>%
  gather(key = variable, value = value, -rgt, -alpha_prior) 


### Plot
title <- "Posterior densities of site-specific estimates"


ggplot(data = df_plot, aes(x = value, group = alpha_prior, fill = alpha_prior)) +
  geom_density(position = "identity", size = 0.1, alpha = 0.3) + 
  geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
  labs(title = title, subtitle = subtitle, x = NULL) + theme_preset +
  facet_wrap(variable ~ rgt, scales = "free") + 
  theme_minimal()


