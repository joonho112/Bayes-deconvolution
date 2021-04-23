
###'######################################################################
###'
###' Category: Real Data Analysis
###' 
###' Task: Summarize effect size estimates and their uncertainty 
###' 
###' Data: Simulated data
###' 
###' Date: 2020-05-02 
###'       2021-04-20  updated
###'
###' 
###' Author: JoonHo Lee (`joonho@berkeley.edu`)
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
data_dir <- file.path(work_dir, "datasets", 
                      "13_Real data analysis")


### Call libraries
library(tidyverse)
library(cowplot)
library(metafor)
library(collateral)
library(DPpackage)
library(hhsim)
library(HETOP)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)
source(file = "functions/03_Simulation Implementation_04_Plot helpers.R")
source(file = "functions/03_Simulation Implementation_05_Data management helpers.R")



###'######################################################################
###'
###' Load the cleaned dataset
###'
###'

setwd(data_dir)

load(file = "Conditional_Cash_Transfer_Data_cleaned.rda") 

df <- df_to_save; rm(df_to_save)



###'######################################################################
###'
###' Prepare loops or subsetting the dataset - Only for the Suba District 
###'
###'

### Loop over all, upper secondary, and lower secondary grade students (i)
vec_Grades_label <- c("All Grades", 
                      "Upper Secondary (9-11)")

vec_all_G <- paste0(seq(6, 11, 1), "th")
vec_upper_G <- paste0(seq(9, 11, 1), "th")

list_Grades <- list(vec_all_G, 
                    vec_upper_G)


### Loop over outcome variables (j) 
vec_ynames <- c("enrolled_ontime", "Heldback", "Dropout", "Graduate_ICFES", 
                "MED_Enroll_SPADIES", "MED_University", "MED_Vocational", "MED_Unclassified", 
                "LON_Enroll_SPADIES", "LON_Enroll_ontime", "LON_Enroll_ontime_uncond", 
                "LON_Grad", "LON_Grad_uncond")

vec_ylabs <- c("On-time secondary enrollment", 
               "Held back", 
               "Dropout", 
               "Secondary graduation", 
               "Tertiary enrollment (Medium-term), any time", 
               "Tertiary enrollment (Medium-term): University", 
               "Tertiary enrollment (Medium-term): Vocational", 
               "Tertiary enrollment (Medium-term): Unclassified", 
               "Tertiary enrollment (Long-term), any time", 
               "Tertiary enrollment (Long-term), on time", 
               "Tertiary enrollment (Long-term), on time - unconditional", 
               "Tertiary graduation (Long-term)", 
               "Tertiary graduation (Long-term), unconditional")

vec_treat <- c("Tertiary Treatment")


###'######################################################################
###'
###' Start loops
###'
###'

for (i in seq_along(list_Grades)){  # (1) Loop over grade levels (i)
  
  for (j in seq_along(vec_ynames)){ # (2) Loop over outcome variables (j)
    
    
    ### Extract vector/list elements
    Grades_lab <- vec_Grades_label[i]
    vec_Grades <- list_Grades[[i]]
    y_name <- vec_ynames[j]
    y_lab <- vec_ylabs[j]
    
    treat_name <- c("Suba_TRT")
    treat_sort <- vec_treat
    District_name <- c("Suba")
    
    
    
    ###'######################################################################
    ###'
    ###' Create a folder to save results
    ###' 
    ###'       
    
    ### Create a folder
    setwd(work_dir)
    
    folder_dir <- file.path("datasets",
                            "13_Real data analysis", 
                            paste(sprintf("%02d", 1), District_name,
                                  sprintf("%02d", i), Grades_lab,  
                                  sprintf("%02d", j), y_name, 
                                  sep = "_"))
    
    dir.create(folder_dir, showWarnings = FALSE)
    
    
    ### Reset working directory
    setwd(file.path(work_dir, folder_dir))
    
    
    
    ###'######################################################################
    ###'
    ###' Prepare a dataframe for model fitting
    ###' 
    ###' df_analytic sample
    ###'
    ###'
    
    ### Filter only relevant district and grades
    df_temp <- df %>%
      filter(District == District_name) %>%
      filter(Grade %in% vec_Grades)
    
    
    ### Rename the outcome variable
    idxl <- names(df_temp) == y_name
    names(df_temp)[idxl] <- "y"
    
    
    ### Rename the treatment variable
    idxl <- names(df_temp) == treat_name
    names(df_temp)[idxl] <- "treat"
    
    
    ### Remove missing outcomes
    idx_miss <- is.na(df_temp$y)
    
    df_ymiss <- data.frame(table(idx_miss)) %>%
      mutate(Percent = 100*(Freq/nrow(df_temp)))
    
    write.csv(df_ymiss, "Table_01_Percent of missing outcomes.csv")
    
    df_temp <- df_temp[!idx_miss, ]
    
    
    ### Save the subsetted dataset
    dualsave(df_temp, "df_analytic sample")
    
    
    
    ###'######################################################################
    ###'
    ###' [Level-1 model] Fit the logistic regression models per group
    ###'
    ###'
    
    ### Convert the "treatment" variable into factor type
    classmode(df_temp, treat)
    tabdf(df_temp, treat)
    df_temp$treat <- relevel(df_temp$treat, ref = "Control")
    
    
    ### Filter only schools with 15 or more cluster size
    df_valid <- df_temp %>%
      group_by(School_Code) %>%
      mutate(N = n()) %>%
      filter(N >= 10) %>%
      filter(!is.na(School_Code))
    
    vec_ID <- unique(df_valid$School_Code)
    
    
    ### Fit logistic regression model and collect estimates
    
    list_temp <- list()
    
    
    for (k in seq_along(vec_ID)){
      
      df_fit <- df_valid %>%
        filter(School_Code == vec_ID[k])
      
      mod_fit <- glm(y ~ treat, data = df_fit, family = "binomial")
      
      df_coef <- summary(mod_fit) %>% 
        {.$coefficients} %>%
        data.frame() %>% 
        rownames_to_column() %>%
        mutate(School_Code =  vec_ID[k]) %>%
        dplyr::select(School_Code, everything())
      
      list_temp[[k]] <- df_coef
      
    }
    
    df_collect <- bind_rows(list_temp) %>%
      set_names(c("School_Code", "Coef", "Estimate", "SE", "z", "p"))
    
    df_es <- df_collect %>%
      filter(Coef != "(Intercept)")
    
    ### Save the estimated effect sizes
    dualsave(df_es, "effect size estimates")
    
    
    
    ###'######################################################################
    ###'
    ###' [Level-2 model] Fit Bayesian deconvolution models
    ###'
    ###'
    
    ### Define necessary parameters to control Bayesian model estimation
    nmcmc <- 2000    # number of MCMC iterations per each model fit
    nburn <- 2000    # burn-in iterations to discard at each round
    nDraws <- nmcmc + nburn
    
    
    
    ###'######################################################################
    ###'
    ###' (1) Gaussian model
    ###' 
    ###' 
    
    outp <- GaussianGaussian(df_es$Estimate, 
                             df_es$SE, 
                             nDraws)
    
    df_posterior1 <- get_posterior_hhsim(output = outp, nburn, nDraws)
    
    
    
    ###'######################################################################
    ###'
    ###' (2) T model
    ###' 
    ###' 
    
    outp <- GaussianT(df_es$Estimate, 
                      df_es$SE, 
                      nDraws)

    df_posterior2 <- get_posterior_hhsim(output = outp, nburn, nDraws)

    
    
    ###'######################################################################
    ###'
    ###' (3) DP mixture model with diffuse prior
    ###' 
    ###' 
    
    ### Prepare dataset as a matrix form (with y and sigma2)
    mat_DPmeta <- df_es %>% 
      dplyr::select(Estimate, SE) %>% 
      set_names(c("Y", "sigma2")) %>%
      as.matrix()
    
    
    ### Set MCMC parameters
    state <- NULL    # the current value of the parameters
    mcmc <- list(nburn = 4000,    # the number of burn-in scans, 
                 nsave = 4000,    # the total number of scans to be saved,
                 nskip = 20,      # the thinning interval, 
                 ndisplay = 100)  # the number of saved scans to be displayed on screen
    
    
    ### Diffuse prior with respect to \alpha (the precision parameter)
    b <- 0.1
    alpha_mean <-  nrow(mat_DPmeta)/2
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
    
    df_posterior3 <- get_posterior_DPmeta(outp, nburn_DP = 4000)
    
    
    
    ###'######################################################################
    ###'
    ###' (4) DP mixture model with informative prior
    ###'
    ###'

    ### Prepare dataset as a matrix form (with y and sigma2)
    mat_DPmeta <- df_es %>% 
      dplyr::select(Estimate, SE) %>% 
      set_names(c("Y", "sigma2")) %>%
      as.matrix()
    
    
    ### Set MCMC parameters
    state <- NULL    # the current value of the parameters
    mcmc <- list(nburn = 4000,    # the number of burn-in scans, 
                 nsave = 4000,    # the total number of scans to be saved,
                 nskip = 20,      # the thinning interval, 
                 ndisplay = 100)  # the number of saved scans to be displayed on screen
    
    
    ### Informative prior with respect to \alpha (the precision parameter)
    info_priors <- list(c(1.24, 0.64), 
                        c(1.60, 1.22), 
                        c(1.84, 1.82),
                        c(1.96, 2.38))
    
    names(info_priors) <- c("25", "50", "100", "200")
    
    K_max <- names(info_priors) %>%
      as.numeric() - nrow(mat_DPmeta) 
    
    idx <- K_max == abs(K_max) %>% min()
    
    K_max <- names(info_priors)[idx]
    
    ab_info <- info_priors[[as.character(K_max)]]
    
    a <- ab_info[1]; b <- ab_info[2]
    
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
    df_posterior4 <- get_posterior_DPmeta(outp, nburn_DP = 4000)
    
    
    
    ###'######################################################################
    ###'
    ###' Collect and save posterior samples
    ###'
    ###'
    
    list_collect <- list(df_posterior1, 
                         df_posterior2, 
                         df_posterior3, 
                         df_posterior4)
    
    names(list_collect) <- c("Gaussian", "T", "DP-diffuse", "DP-inform")
    
    saveRDS(list_collect, file = "list_collection of posterior samples.rds")
    
    
    
    ###'######################################################################
    ###'
    ###' Compute the site-specific estimates 
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
    
    list_est <- list()
    
    for (l in seq(length(list_collect))){
      
      ### Compute PM, PSD, CB, GR, rbar, and rhat
      df_posterior <- list_collect[[l]]
      
      df_theta <- df_posterior %>%
        dplyr::select(contains("tau_k["))  # Extract only site-specific estimates
      
      df_estimates <- HETOP::triple_goal(as.matrix(df_theta))
      
      
      ### Combine true theta and ML theta (maximum likelihood) with true ranks
      df_meta <- df_es %>%
        dplyr::select(School_Code, Estimate, SE)
        
      df_est <- cbind.data.frame(df_meta, df_estimates) %>%
        dplyr::select(-index) %>%
        rename(theta_ml = Estimate) %>%
        mutate(rpm = rank(theta_pm),
               rcb = rank(theta_cb), 
               rgr = rank(theta_gr), 
               rml = rank(theta_ml))
      
      
      ### Save as a list component 
      list_est[[l]] <- df_est
    }
    
    
    ### Save as .rds file
    names(list_est) <- c("Gaussian", "T", "DP-diffuse", "DP-inform")
    saveRDS(list_est, file = "list_collection of posterior summary estimates.rds")
    
    
    
    ###'######################################################################
    ###'
    ###' Tidy up the posterior summary estimates
    ###'
    ###'
    
    for (m in seq(length(list_est))){
      
      df_temp <- list_est[[m]]
      
      df_temp <- df_temp %>%
        mutate(G_prior = names(list_est)[m]) %>%
        dplyr::select(G_prior, everything())
      
      list_est[[m]] <- df_temp
      
    }
    
    df_final <- bind_rows(list_est) %>%
      mutate(G_prior = factor(G_prior, 
                              levels = c("Gaussian", "T", "DP-diffuse", "DP-inform")))
    
    
    dualsave(df_final, "df_site-level effect distributions")
    
    
    ###'######################################################################
    ###'
    ###' Compare distributions
    ###'
    ###'

    ### Generate dataframe to plot
    vec_est <- c("ML", "PM", "CB", "GR")
    vec_lab <- c("Raw estimates", 
                 "Posterior means", 
                 "Constrained Bayes estimates", 
                 "Triple-goal estimates")
    
    vec_est_vars <- paste0("theta_", str_to_lower(vec_est))
    
    df_plot <- df_final %>%
      dplyr::select(G_prior, theta_ml, theta_pm, theta_cb, theta_gr) %>%
      gather(key = estimator, value = value, all_of(vec_est_vars)) %>%
      mutate(estimator = factor(estimator, 
                                levels = vec_est_vars, 
                                labels = vec_lab)) %>%
      filter(estimator != "Raw estimates") %>%
      filter(G_prior != "T")
      
    
    ### Plot!
    color_manu <- color_palette[seq(unique(df_plot$estimator))]
    shape_manu <- shape_palette[seq(unique(df_plot$estimator))]
    
    p <- ggplot(data = df_plot, aes(x = value, 
                                    group = G_prior, 
                                    color = G_prior)) + 
      geom_density(adjust = 1, size = 0.7) + 
      geom_vline(xintercept = 0, 
                 color = "black", linetype = "dashed", size = 0.5) + 
      facet_wrap(.~ estimator, ncol = 1, scales = "free_y") + 
      scale_color_manual(values = color_manu) +
      scale_shape_manual(values = shape_manu) + 
      theme_classic() +
      theme(strip.background = element_blank()) + 
      theme(legend.position = "right", legend.direction = "vertical") + 
      labs(title = paste0("Outcome: ", y_lab),
           subtitle = paste0(District_name, " district, ", treat_sort), 
           caption = NULL, 
           y = "Density",  
           x = "Site-level effects (logit scale)") 

    
    ### Save the resulting graph
    ggsave("Figure_01_Compare distributions of estimated G_.pdf", 
           plot = p, 
           width = 7, height = 8)
  
    
    
    ###'######################################################################
    ###' 
    ###' Print the progress  
    ###' 
    
    cat(paste(sprintf("%02d", i), Grades_lab,  
              sprintf("%02d", j), y_lab, 
              sep = "_"), 
        " => completed", "\n")
    
  }    # End of loop over grade levels (j) 
}      # End of loop over districts (i)

