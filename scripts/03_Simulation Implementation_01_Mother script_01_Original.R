
###'######################################################################
###'
###' Category: Simulation Implementation
###' 
###' Task: Simulation (Mother script)  
###' 
###'     - Flexible distributions for triple-goal estimates in 
###'       two-stage hierarchical models
###'       
###'     - Using Dirichlet Processes for Modeling Heterogenous 
###'       Treatment Effects Across Sites  
###'       
###' Data: Simulated data
###' 
###' Data: 2020-03-26
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
library(sn)
library(hhsim)
library(DPpackage)
library(HETOP)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "scripts/03_Simulation Implementation_04_Plot helpers.R")
source(file = "scripts/03_Simulation Implementation_05_Data management helpers.R")



###'######################################################################
###'
###' Prepare simulation variables
###' 
###' (1) N: number of sites (K = 25, 50, 100, 200)
###'
###' (2) rgt: the informativeness of the data (rgt = 0.10, 0.33, 1.00)
###' 
###' (3) rsr: the heterogeneity of the uncertainty estimates (rsr = 1^2, 5^2, 10^2)
###' 
###' (4) True population distribution of G
###'     - 1. Gaussian distribution N(0, 1)
###'     - 2. T distribution T_5(0, 1) with df = 5
###'     - 3. Gaussian mixture 0.8 N(0, 1) + 0.2 (4, 1)
###'
###'  => Change these variables different from Paddock et al. (2006)
###'  
###'  

### Set seed 
set.seed(10101)


### Generate a table of the simulation parameters
tbl_param <- expand.grid(truth = c("Gaussian", "T", "Mix", "Skew"),
                         assumed = c("Gaussian", "T", "DP"),
                         rgt = c(0.1, 0.33, 1, 10),
                         rsr = c(1, 5, 10),
                         n = c(25, 50, 100, 200), 
                         df = 5)


### Reorder the table & add a variable
tbl_param <- tbl_param %>%
  arrange(n, rgt, rsr, truth, assumed)


### Add a concatenated design indicator
sim_parm <- tbl_param %>%
  unite(fname, truth, assumed, rgt, rsr, n, sep = "_") %>%
  dplyr::select(-df) %>%
  cbind.data.frame(tbl_param) %>%
  select(truth:df, fname)

head(sim_parm, 10)


### Save as a table 
write.csv(sim_parm, file = "tables/simulation parameters.csv")


### The number of Monte Carlo iterations to perform
N_sim <- 100



###'######################################################################
###'
###' Load functions for data generation & converting quantiles to CDF
###' 
###' (1) Generate_Gaussian => Normal
###' (2) Generate_T => heavy-tailed
###' (3) Generate_Mixture => Multi- (Bi-) modal
###' (4) Generate_SkewN => Skewness
###'
###'

setwd(work_dir)
source(file = "scripts/03_Simulation Implementation_02_Data generation.R")
source(file = "scripts/03_Simulation Implementation_03_Convert quantiles to theoretical CDFs.R")



###'######################################################################
###'
###' Start loops over simulation conditions
###'
###'

### The number of simulation conditions 
N_cond <- length(sim_parm$truth)   # Total 576 conditions


### Mark start time
start_time <- Sys.time()


for (j in 1:N_cond){
  
  ###'######################################################################
  ###'
  ###' Set required objects for simulation
  ###' 
  ###'
  
  ### Extract simulation condition elements
  gen_name <- paste0("True_", sim_parm$truth[j])
  anl_name <- paste0("Model_", sim_parm$assumed[j])
  rgt <- sim_parm$rgt[j]  
  rsr <- sim_parm$rsr[j]
  nu <- sim_parm$df[j]
  K_max <- sim_parm$n[j]
  
  
  ### Set quantiles (tail probabilities)
  tail_p <- c(0.05, 0.10, 0.25, 0.75, 0.90, 0.95) 
  
  
  ### Generate an empty object to contain variances of true thetas
  var_theta <- NULL
  
  
  ### Define necessary parameters to control Bayesian model estimation
  nmcmc <- 2000    # number of MCMC iterations per each model fit
  nburn <- 2000    # burn-in iterations to discard at each round
  nDraws <- nmcmc + nburn
  
  
  ### Define the number of posterior summary methods
  post_sum_methods <- c("PM", "CB", "GR", "ML")
  J_max <- length(post_sum_methods) 
  
  
  
  ###'######################################################################
  ###'
  ###' Generate empty arrays or vectors 
  ###' to save or accumulate resulting losses
  ###' 
  ###' 
  
  ### (0) Means and variances of PM, CB, GR, and ML
  temp_names <- c(paste0("mean_", post_sum_methods), 
                  paste0("var_", post_sum_methods))
  
  array_MOM_theta <- array(NA, c(N_sim, length(temp_names)))
  colnames(array_MOM_theta) <- temp_names
  
  
  ### (1) Squared error loss (SEL) for the estimated "ranks"
  temp_names <- c("SEL_rank", paste("SEL_rank", post_sum_methods, sep = "_"))
  array_SEL_rank <- array(NA, c(N_sim, length(temp_names)))
  colnames(array_SEL_rank) <- temp_names
  
  
  ### (2) Squared error loss (SEL) for the "individual site-specific effects"
  temp_names <- c(paste0("SSEL_", post_sum_methods), 
                  paste0("WSSEL_", post_sum_methods))
  array_SSEL_theta <- array(NA, c(N_sim, length(temp_names)))
  colnames(array_SSEL_theta) <- temp_names  
  
  sim_SEL <- sim_WSEL <- rep(0, K_max*J_max) 
  sim_SSEL <- sim_WSSEL <- rep(0, J_max)
  
  
  ### (3) Integrated SEL (ISEL) for the estimated "empirical distribution function" 
  temp_names <- c(paste0("ISEL_", post_sum_methods))
  array_ISEL_EDF <- array(NA, c(N_sim, length(temp_names)))
  colnames(array_ISEL_EDF) <- temp_names
  
  
  ### (4) Tail quantiles
  tail_p <- c(0.05, 0.10, 0.25, 0.75, 0.90, 0.95) 
  sim_tail <- rep(0, length(tail_p)*J_max)
  
  
  ### (5) Histogram bars
  n_bars <- 50
  sim_hisbars <- rep(0, n_bars*J_max)
  
  
  
  ###'######################################################################
  ###'
  ###' Perform Mote Carlo iterations (N_sim)
  ###'
  ###'
  
  for (i in seq(N_sim)){
    
    ### Print progress
    cat(paste0(sprintf("%03d", j)), 
      gen_name, anl_name, paste0("N = ", K_max), 
        paste0("rgt = ", rgt), paste0("rsr = ", rsr), 
        paste0("i = ", i), "\n")
    
    
    
    ###'######################################################################
    ###'
    ###' (1) Generate data 
    ###' 
    ###' - Conforming to the "true" distribution
    ###'
    ###'
    
    if (gen_name == "True_Gaussian"){ 
      
      list_meta <- Generate_Gaussian(n = K_max, rgt = rgt, rsr = rsr)
      df_meta <- list_meta %>% list_to_df_G()
      priquan <- Q2CDF_Gaussian(tail_p)
      
      
    } else if (gen_name == "True_T"){   
      
      nu <- sim_parm$df[j]  # degrees of freedom
      
      list_meta <- Generate_T(n = K_max, rgt = rgt, rsr = rsr, nu = nu)
      df_meta <- list_meta %>% list_to_df_G()
      priquan <- Q2CDF_T(tail_p)
      
      
    } else if (gen_name == "True_Mix"){
      
      delta <- 4   # distance between two means
      eps <- 0.2   # proportion of the small component
      ups <- 1     # ratio between two variances
      
      list_meta <- Generate_Mixture(n = K_max, rgt = rgt, rsr = rsr, 
                                    delta = delta, eps = eps, ups = ups)
      
      df_meta <- list_meta %>% list_to_df_G()
      priquan <- Q2CDF_Mixture(tail_p)
      
      
    } else if (gen_name == "True_Skew"){
      
      slant <- 5   # control the degree of skewness
      
      list_meta <- Generate_SkewN(n = K_max, rgt = rgt, rsr = rsr, 
                                  mean = 0, var = 1, slant = slant)
      
      df_meta <- list_meta %>% list_to_df_G()
      priquan <- Q2CDF_Mixture(tail_p)
      
    }
    
    ### Get the variance of true theta
    var_theta <- c(var_theta, c(var(list_meta$theta)))
    
    
    
    ###'######################################################################
    ###'
    ###' (2) Fit the assumed model 
    ###' 
    ###' - To get posterior samples 
    ###'
    ###' 
    
    if (anl_name == "Model_Gaussian"){
      
      ### Fit the model
      outp <- GaussianGaussian(list_meta$Y,
                               list_meta$sigma2,
                               nDraws)
      
      ### Extract posterior samples
      df_posterior <- get_posterior_hhsim(output = outp, nburn, nDraws)
      
      
    } else if (anl_name == "Model_T"){
      
      ### Fit the model
      outp <- GaussianT(list_meta$Y,
                        list_meta$sigma2,
                        nDraws, nu)
      
      ### Extract posterior samples
      df_posterior <- get_posterior_hhsim(output = outp, nburn, nDraws)
      
      
    } else if (anl_name == "Model_DP"){
      
      ### Prepare dataset as a matrix form (with y and sigma2)
      mat_DPmeta <- df_meta %>% select(Y, sigma2) %>% as.matrix()
      
      
      ### Set MCMC parameters
      state <- NULL    # the current value of the parameters
      mcmc <- list(nburn = 4000,    # the number of burn-in scans, 
                   nsave = 4000,    # the total number of scans to be saved,
                   nskip = 20,      # the thinning interval, 
                   ndisplay = 100)  # the number of saved scans to be displayed on screen
      
      
      ### Diffuse prior with respect to \alpha (the precision parameter)
      b <- 0.1
      alpha_mean <- K_max/2
      a <- alpha_mean*b
      alpha_var <- a/b^2  
      
      prior <- list(a0 = a,      # alpha0: shape param.
                    b0 = b,      # alpha0: rate param.  
                    tau1 = 1,    # G0 variance: shape param.
                    tau2 = 1,    # G0 variance: rate param. 
                    mub = 0, 
                    Sb = 100)    # G0 mean: mean 
      
      
      ### Fit the model with the DP prior
      outp <- DPmeta(formula = mat_DPmeta ~ 1,
                     prior = prior, mcmc = mcmc, 
                     state = state, status = TRUE)
      
      
      ### Extract posterior samples
      df_posterior <- get_posterior_DPmeta(outp, nburn_DP = 4000)

    }
    
    
    
    ###'######################################################################
    ###'
    ###' (3) Compute the site-specific estimates
    ###' 
    ###' - Posterior means (PM) and standard deviations (PSD)
    ###' - Constrained Bayes estimates (CB)
    ###' - Triple goal estimates (GR)
    ###' - Posterior means of ranks (rbar)
    ###' - Integer ranks (rhat)
    ###'
    ###' and add the TRUE site-specific values etc. 
    ###' 
    ###'
    
    ### Compute PM, PSD, CB, GR, rbar, and rhat
    df_theta <- df_posterior %>%
      dplyr::select(contains("tau_k["))  # Extract only site-specific estimates
    
    df_temp <- HETOP::triple_goal(as.matrix(df_theta))
    
    
    ### Combine true theta and ML theta (maximum likelihood) with true ranks
    df_est <- cbind.data.frame(df_meta, df_temp) %>%
      dplyr::select(-index) %>%
      rename(theta_true = theta, theta_ml = Y) %>%
      mutate(rtrue = rank(theta_true), 
             rpm = rank(theta_pm),
             rcb = rank(theta_cb), 
             rgr = rank(theta_gr), 
             rml = rank(theta_ml))
    
    
    ### Calculate moments (mean and variance) of PM, CB, GR, ML
    est_vars <- paste0("theta_", str_to_lower(post_sum_methods)) 
    theta_est <- df_est %>% dplyr::select(est_vars) %>%
      as.matrix()  # Estimated theta's
    
    array_MOM_theta[i, ] <- c(apply(theta_est, 2, mean), 
                              apply(theta_est, 2, var))
    
    
    
    ###'######################################################################
    ###'
    ###' (4) Apply loss functions
    ###' 
    ###'  4-1. Squared error loss (SEL) for the estimated "ranks"
    ###' 
    ###'
    
    SEL_rank <- SEL(df_est$rhat, df_est$rtrue)
    SEL_rank_PM <- SEL(df_est$rpm, df_est$rtrue)
    SEL_rank_CB <- SEL(df_est$rcb, df_est$rtrue)
    SEL_rank_GR <- SEL(df_est$rgr, df_est$rtrue)
    SEL_rank_ML <- SEL(df_est$rml, df_est$rtrue)
    
    array_SEL_rank[i, ] <- c(SEL_rank, SEL_rank_PM, SEL_rank_CB, 
                             SEL_rank_GR, SEL_rank_ML) 
    
    
    
    ###'######################################################################
    ###'
    ###' (4) Apply loss functions
    ###' 
    ###'  4-2. SEL for the "Individual site-specific effects"
    ###'
    ###' 
    
    ###' Obtain Bayes risk (SEL and SSEL) for 
    ###' all individual site-specific effects estimates 
    ###' at this simulation iteration
    list_SEL <- get_theta_SEL(df_est = df_est, 
                              est_vec = post_sum_methods, 
                              weight = rep(1, K_max))
    
    array_SSEL_theta[i, ] <- list_SEL[[2]] %>% unlist()
    
    
    ### Accumulate SEL and SSEL over simulation iterations
    temp <- theta.simu(J_max, K_max, 
                       list_SEL[[1]]$SEL, list_SEL[[2]]$SSEL, 
                       sim_SEL, sim_SSEL, 
                       list_SEL[[1]]$WSEL, list_SEL[[2]]$WSSEL, 
                       sim_WSEL, sim_WSSEL)
    
    sim_SEL <- temp$simusel
    sim_WSEL <- temp$simuwsel
    sim_SSEL <- temp$simussel
    sim_WSSEL <- temp$simuwsse
    
    
    
    ###'######################################################################
    ###'
    ###' (4) Apply loss functions
    ###' 
    ###'  4-3. Integrated Squared Error Loss (ISEL) for  
    ###'       the Empirical Distribution Function (EDF)
    ###'
    ###' 
    
    vec_ISEL <- get_EDF_ISEL(df_est)
    
    array_ISEL_EDF[i, ] <- vec_ISEL
    
    
    
    ###'######################################################################
    ###'
    ###' (4) Apply loss functions
    ###' 
    ###'  4-4. Get tail quantile estimates
    ###'
    ###' 
    
    tail <- get_tail_quantiles(df_est, tail_p, 
                               priquan, est_vec = post_sum_methods)

    sim_tail <- sim_tail + tail
    
    
    
    ###'######################################################################
    ###'
    ###' (4) Apply loss functions
    ###' 
    ###'  4-5. the Empirical Distribution Function (EDF)
    ###'
    ###' 
  
    mat_histbars <- get_histogram_bars(df_est, 
                                       est_vec = post_sum_methods, 
                                       n_bars = n_bars, 
                                       alowhis = -6, 
                                       ahighhis = 6)
    
    sim_hisbars <- sim_hisbars + mat_histbars

    
  } # End of loops over i (Monte Carlo iteration)
  
  
  
  ###'######################################################################
  ###'
  ###' Tidy up and save the results from one condition (j)
  ###'
  ###'
  
  ### Create a directory to save the results
  folder_name <- paste0(sprintf("%03d", j), "_", sim_parm$fname[j])
  
  folder_dir <- file.path(work_dir, "datasets", 
                          "04_Simulation results", 
                          folder_name)
  
  dir.create(folder_dir, showWarnings = FALSE)
  
  setwd(folder_dir)
  
  
  
  ###'######################################################################
  ###'
  ###' (1) Means and variances of PM, CB, GR, and ML
  ###'
  ###'
  
  ### Save the dataframe 
  array_MOM_theta %>%
    data.frame() %>%
    write.csv(file = "01_Means and Variances.csv")
  
  
  ### Compute the means of the means and variances of estimators 
  avg_MOM <- apply(array_MOM_theta, 2, mean)
  
  
  
  ###'######################################################################
  ###'
  ###' (2) Sum of Squared Error Loss (SSEL) for 
  ###'     the individual site-specific effects
  ###'
  ###'
  
  ### Save the data.frame
  array_SSEL_theta %>%
    data.frame() %>%
    write.csv(file = "02_SSEL for the individual effects.csv")
  
  
  ### Compute the means of SSEL
  avg_SSEL_theta <- apply(array_SSEL_theta, 2, mean)
  
  
  
  ###'######################################################################
  ###'
  ###' (3) Squared Error Loss (SEL) for the
  ###'     the individual site-specific effects
  ###'
  ###'
  
  ### Tidy up the results
  avg_SEL_theta <- matrix(sim_SEL/N_sim, nrow = K_max, ncol = J_max) %>%
    data.frame() 
  
  avg_WSEL_theta <- matrix(sim_WSEL/N_sim, nrow = K_max, ncol = J_max) %>%
    data.frame()
  
  names(avg_SEL_theta) <- paste0("SEL_", post_sum_methods)
  names(avg_WSEL_theta) <- paste0("WSEL_", post_sum_methods)
  
  
  ### Save the dataframes
  avg_SEL_theta %>%
    write.csv(file = "03_SEL for the individual effects.csv")
  
  avg_WSEL_theta %>%
    write.csv(file = "03_WSEL for the indificual effects.csv")
  
  
  
  ###'######################################################################
  ###'
  ###' (4) Integrated SEL (ISEL) for 
  ###'     the estimated "empirical distribution function (EDF)"
  ###'
  ###'
  
  ### Save the dataframe
  array_ISEL_EDF %>%
    data.frame() %>%
    write.csv(file = "04_ISEL for the EDF.csv")
  
  
  ### Compute the mean of the ISEL
  avg_ISEL_EDF <- apply(array_ISEL_EDF, 2, mean)
  
  
  
  ###'######################################################################
  ###'
  ###' (5) Squared Error Loss (SEL) for 
  ###'     the estimated ranks
  ###'
  ###'
  
  ### Save the dataframe
  array_SEL_rank %>%
    data.frame() %>%
    write.csv(file = "05_SEL for the ranks.csv")
  
  
  ### Compute the mean
  avg_SEL_rank <- apply(array_SEL_rank, 2, mean)
  
  
  
  ###'######################################################################
  ###'
  ###' (6) Tail quantiles
  ###'
  ###'
  
  ### Tidy up and save the result
  avg_tail <- data.frame(rep(post_sum_methods, times = length(tail_p)), 
                         rep(tail_p, each = J_max), 
                         sim_tail/N_sim)
  
  names(avg_tail) <- c("estimator", "tail_p", "quantile")
  
  avg_PCTs <- avg_tail %>% 
    mutate(estimator = factor(estimator, levels = post_sum_methods), 
           diff = round(quantile - tail_p, 10)) %>%
    arrange(estimator) 
  
  write.csv(avg_PCTs, file = "06_Percentile estimates of G.csv")
  
  
  
  ###'######################################################################
  ###'
  ###' (7) Histogram bars
  ###'
  ###'
  
  ### Tidy up and save the result
  df_hist_bars <- sim_hisbars/(N_sim*K_max)
  
  write.csv(df_hist_bars, file = "07_Histogram bars.csv")
  
  
  
  ###'######################################################################
  ###'
  ###' (8) Variances of true theta
  ###' 
  ###'
  
  df_var_theta <- var_theta %>% data.frame()
  
  names(df_var_theta) <- c("var_true_theta")
  
  df_var_theta %>%
    write.csv(file = "08_Variance of true theta over iterations.csv")
  
  var_true <- mean(var_theta)
  
  names(var_true) <- c("var_true")
  
  
  
  ###'######################################################################
  ###'
  ###' (9) Collect average losses
  ###'
  ###'

  ### Collect loss estimates as a vector 
  vec_collect <- c(avg_MOM, var_true, 
                   avg_SSEL_theta, avg_ISEL_EDF, 
                   avg_SEL_rank[-1])
  
  
  ### Tidy up as a dataframe format and save it
  df_collect <- data.frame(vec_collect) %>% 
    rownames_to_column() %>% 
    mutate(rowname = str_replace(rowname, "SEL_rank", "SELrank")) %>%
    separate(rowname, into = c("quantity", "estimator"), sep = "_") %>%
    rename(value = vec_collect)
  
  write.csv(df_collect, file = "09_Collection of loss estimates.csv")
  
  
}  # End of loops over j (simulation condition)



###'######################################################################
###'
###' Measure computation time
###'
###'

### Mark end time
end_time <- Sys.time()

comp_time <- end_time - start_time

comp_time


