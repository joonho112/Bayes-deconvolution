
###'######################################################################
###'
###' Category: 
###' 
###' Task: (1) Dirichlet Process Mixtures
###'           
###'       - Comparing posterior distributions of various quantities
###'         generated from different strategies for selecting 
###'         the precision parameter \alpha
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
###' Generate a list containing different prior choices
###'
###'

list_priors <- list() 
prior_vec <- c("DORO", "Diffuse", "EB", "IS")


### (1) Dorazio (2009)'s method, chosen using the Kullback-Leibler divergence
list_priors[[1]] <- list(a0 = 0.392,  # alpha0: shape param.
                         b0 = 0.004,  # alpha0: rate param.  
                         tau1 = 1,    # G0 variance: shape param.
                         tau2 = 1,    # G0 variance: rate param. 
                         mub = 0, 
                         Sb = 100)    # G0 mean: mean 


### (2) Diffuse prior with respect to \alpha
b <- 0.1
alpha_mean <- N/2
a <- alpha_mean*b
alpha_var <- a/b^2  

list_priors[[2]] <- list(a0 = a,      # alpha0: shape param.
                         b0 = b,      # alpha0: rate param.  
                         tau1 = 1,    # G0 variance: shape param.
                         tau2 = 1,    # G0 variance: rate param. 
                         mub = 0, 
                         Sb = 100)    # G0 mean: mean 


###' (3) Empirical Bayes estimate of alpha (Dorazio et al., 2008)
list_priors[[3]] <- list(alpha = 4.05,  # fixed precision parameter
                         tau1 = 1,   # G0 variance: shape param.
                         tau2 = 1,   # G0 variance: rate param. 
                         mub = 0, 
                         Sb = 100)   # G0 mean: mean 


###' (4) alpha chosen via the importance sampling approach (Antonelli et al., 2016)
list_priors[[4]] <- list(alpha = 3.90,  # fixed precision parameter
                         tau1 = 1,   # G0 variance: shape param.
                         tau2 = 1,   # G0 variance: rate param. 
                         mub = 0, 
                         Sb = 100)   # G0 mean: mean 



###'######################################################################
###'
###' Construct a loop over the alpha selection methods
###' 
###' 

### Set MCMC parameters
state <- NULL    # the current value of the parameters
nburn <- 4000    # the number of burn-in scans
nsave <- 4000    # the total number of scans to be saved
nskip <- 20      # the thinning interval
ndisplay <- 100  # the number of saved scans to be displayed on screen

mcmc <- list(nburn = nburn, nsave = nsave,
             nskip = nskip, ndisplay = ndisplay)


### Prepare dataset as a matrix form (with y and sigma2)
mat_DPmeta <- df_G %>% select(tau_k_hat, se_k2) %>% as.matrix()


### Start a loop

list_collect <- list()

for (i in seq_along(list_priors)){
  
  ### Extract the loop element
  prior <- list_priors[[i]]
  
  ### Estimate the Dirichlet Process model using DPpackage 
  outp <- DPmeta(formula = mat_DPmeta ~ 1,
                 prior = prior, mcmc = mcmc, 
                 state = state, status = TRUE)
  
  ### Save the output
  list_collect[[i]] <- outp
}



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

for (i in seq_along(list_priors)){
  
  ##  Get posterior samples
  meta_list[[i]] <- get_posterior_DPmeta(list_collect[[i]], nburn) %>%
    mutate(alpha_choice = prior_vec[i]) %>%
    dplyr::select(alpha_choice, everything())
  
}

df_combine <- bind_rows(meta_list) %>% 
  mutate(alpha_choice = factor(alpha_choice, levels = prior_vec))

dim(df_combine)



###'######################################################################
###'
###' Get alpha priors - It's simple
###' 
###' => draw the same number of samples as posterior (n = 2,000)
###' 
###' (1) DORO
###' (2) Diffuse 
###'
###'

### Generate the prior distribution: (1) DORO
set.seed(12345)
DORO_prior <- rgamma(n = 2000, shape = 0.392, rate = 0.004)
plot(density(DORO_prior))
mean(DORO_prior); sd(DORO_prior)


### Generate the prior distribution: (2) Diffuse
set.seed(12345)
Diffuse_prior <- rgamma(n = 2000, shape = a, rate = b)
plot(density(Diffuse_prior))
mean(Diffuse_prior); sd(Diffuse_prior)

df_priors <- data.frame(DORO_prior, Diffuse_prior)


### Plot the prior distributions: (1) DORO 
mean_value <- mean(DORO_prior, na.rm = TRUE)
var_value <- var(DORO_prior, na.rm = TRUE)

p <- ggplot(data = df_priors, aes(x = DORO_prior)) +
  geom_density(fill = "red", alpha = 0.1) + 
  geom_vline(xintercept = mean_value, 
             linetype = "dashed", color = "red") + 
  geom_text(aes(x = mean_value, y = 0.001, label = round(mean_value, 2)), 
            color = "red") + 
  labs(title = "Prior distribution of the precision parameter", 
       subtitle = "Gamma(0.392, 0.004) chosen by the Dorazio (2009) method", 
       x = "alpha", 
       caption = paste0("mean = ", round(mean_value, 2), 
                        ", variance = ", round(var_value, 2))) + 
  theme_preset

ggsave("figures/Prior distribution_01_DORO.pdf", p, 
       width = 8, height = 5)


### Plot the prior distributions: (2) Diffuse 
mean_value <- mean(Diffuse_prior, na.rm = TRUE)
var_value <- var(Diffuse_prior, na.rm = TRUE)

p <- ggplot(data = df_priors, aes(x = Diffuse_prior)) +
  geom_density(fill = "blue", alpha = 0.1) + 
  geom_vline(xintercept = mean_value, 
             linetype = "dashed", color = "red") + 
  geom_text(aes(x = mean_value, y = 0.001, label = round(mean_value, 2)), 
            color = "red") + 
  labs(title = "Prior distribution of the precision parameter - Gamma(2.5, 0.1)", 
       subtitle = "Chosen based on the center between 1 and N with a large variance", 
       x = "alpha", 
       caption = paste0("mean = ", round(mean_value, 2), 
                        ", variance = ", round(var_value, 2))) + 
  theme_preset

ggsave("figures/Prior distribution_02_Diffuse.pdf", p, 
       width = 8, height = 5)


### DORO vs. Diffuse priors
df_temp <- df_priors %>% gather(key = "method", value = value)

p <- ggplot(data = df_temp, aes(x = value, group = method, fill = method)) +
  geom_density(alpha = 0.5) + 
  labs(title = "Prior distributions of the precision parameter", 
       subtitle = "DORO vs. Diffuse priors", 
       x = "alpha") + 
  theme_preset

ggsave("figures/Prior distribution_03_DORO vs Diffuse.pdf", p, 
       width = 8, height = 5)



###'######################################################################
###'
###' Prior vs. Posterior
###' 
###' Compare the densities of the precision parameter
###'
###'

### Generate a dataframe to plot - Bind gamma prior samples
Prior <- c(DORO_prior, Diffuse_prior, rep(NA, times = 4000))

df_plot <- df_combine %>%
  dplyr::select(alpha_choice, alpha0) %>%
  arrange(alpha_choice) %>%
  rename(Posterior = alpha0) %>%
  cbind.data.frame(Prior) %>%
  gather(key = sample, value = value, Prior, Posterior) %>%
  mutate(sample = factor(sample, levels = c("Prior", "Posterior"))) %>%
  filter(alpha_choice %in% c("DORO", "Diffuse"))

### Check summary statistics
tbl_temp <- df_plot %>%
  group_by(alpha_choice, sample) %>%
  summarise(mean = mean(value, na.rm = TRUE), 
            var = var(value, na.rm = TRUE))

write.csv(tbl_temp, file = "tables/Prior vs Posterior_DORO and Diffuse.csv")


### Plot
p <- ggplot(data = df_plot, aes(x = value, group = sample, fill = sample)) +
  geom_density(position = "identity", size = 0.1, alpha = 0.3) + 
  geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") +
  facet_wrap(~ alpha_choice, scales = "free") + 
  labs(title = "Precision parameter: Prior vs. Posterior", 
       subtitle = NULL, 
       x = "alpha") + theme_preset + theme_minimal()

ggsave("figures/Prior vs Posterior_DORO and Diffuse.pdf", p, 
       width = 12, height = 5)



###'######################################################################
###'
###' (1) Hyperparameters $\tau$ and $\sigma_{\tau}$
###'
###' - Compare the outcomes from different alpha choices
###'
###'

### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(alpha_choice, G0_mu, G0_s2) %>%
  mutate(G0_s2 = log(G0_s2)) %>%    # Log-transformation for G0_s2
  gather(key = variable, value = value, G0_mu, G0_s2) %>%
  mutate(variable = factor(variable, levels = c("G0_mu", "G0_s2")))


### Plot!
title <- c("Posterior densities of hyperparameters")
subtitle <- c("by the choices of alpha prior")

p <- ggplot(data = df_plot, aes(x = value, group = alpha_choice, fill = alpha_choice)) +
  geom_density(position = "identity", size = 0.1, alpha = 0.3) + 
  geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
  labs(title = title, subtitle = subtitle, x = NULL) + theme_preset +
  facet_wrap( ~ variable, scales = "free") + 
  theme_minimal()

ggsave("figures/By alpha choices_Posterior densities of hyperparameters.pdf", p, 
       width = 12, height = 5)



###'######################################################################
###'
###' (2) Number of clusters
###'
###'

### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(alpha_choice, Ncluster) 


### Plot!
title <- c("Posterior densities of the (assumed) number of clusters")
p <- ggplot(data = df_plot, aes(x = Ncluster, group = alpha_choice, fill = alpha_choice)) +
  geom_histogram(position = "dodge", size = 1, alpha = 0.3) + 
  labs(title = title, subtitle = subtitle) + theme_preset + 
  # facet_wrap(~ alpha_prior, scales = "free") + 
  scale_x_continuous(breaks = seq(from = 0, to = max(df_plot$Ncluster), by = 2))

ggsave("figures/By alpha choices_Posterior densities of the number of clusters.pdf", p, 
       width = 8, height = 5)



###'######################################################################
###'
###' (3) Empirical distribution function (EDF) estimates for $\tau_{k}$s 
###'
###'

### Generate a dataframe to plot
df_plot <- df_combine %>%
  dplyr::select(alpha_choice, contains("tau_k[")) %>%
  gather(key = variable, value = value, -alpha_choice) %>%
  group_by(alpha_choice, variable) %>%
  summarise(mean = mean(value, na.rm = TRUE), 
            sd = sd(value, na.rm = TRUE)) %>%
  mutate(id = as.numeric(str_extract(variable, "\\d+"))) %>%
  arrange(alpha_choice, id) 


### Plot!
p1 <- ggplot(data = df_plot %>% rename(value = mean), 
             aes(x = value, group = alpha_choice, fill = alpha_choice)) +
  geom_density(position = "identity", size = 0.1, alpha = 0.3) + 
  geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
  labs(title = "EDF of posterior means", subtitle = subtitle, x = NULL) + theme_preset +
  theme_minimal()


p2 <- ggplot(data = df_plot %>% rename(value = sd), 
             aes(x = value, group = alpha_choice, fill = alpha_choice)) +
  geom_density(position = "identity", size = 0.1, alpha = 0.3) + 
  geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
  labs(title = "EDF of posterior SDs", subtitle = subtitle, x = NULL) + theme_preset +
  theme_minimal()

p <- plot_grid(p1, p2, ncol = 2, labels = "AUTO")

ggsave("figures/By alpha choices_EDF of posterior means and sds.pdf", p, 
       width = 15, height = 5)



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
  dplyr::select(alpha_choice, rand_var) %>%
  gather(key = variable, value = value, -alpha_choice) 


### Plot
title <- "Posterior densities of site-specific estimates"

p <- ggplot(data = df_plot, aes(x = value, group = alpha_choice, fill = alpha_choice)) +
  geom_density(position = "identity", size = 0.1, alpha = 0.3) + 
  geom_vline(aes(xintercept = 0), size = 0.3, color = "red", linetype = "dashed") + 
  labs(title = title, subtitle = subtitle, x = NULL) + theme_preset +
  facet_wrap( ~ variable, scales = "free") + 
  theme_minimal()

ggsave("figures/By alpha choices_Posterior densities of site-specific estimates.pdf", p, 
       width = 10, height = 5)




