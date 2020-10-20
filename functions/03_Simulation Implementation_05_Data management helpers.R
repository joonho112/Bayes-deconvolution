
###'######################################################################
###'
###' Category: Define functions
###' 
###' Task: Define Data management helper functions
###'
###' Date: 2020-03-29
###' 
###' Author: JoonHo Lee (joonho@berkeley.edu)
###' 
###'

###'######################################################################
###'
###' get_posterior_hhsim()
###' 
###' - Tidy up the output object from hhsim package
###'
###'

get_posterior_hhsim <- function(output = NULL, nburn = 2000, nDraws = 4000){
  
  # (1) Tidy up the "theta" outpututut object
  df_theta <- output[["theta"]] %>%
    unlist() %>%
    matrix(nrow = length(output[["theta"]]), byrow = TRUE) %>%
    data.frame() %>%
    slice(c((nburn + 1):nDraws))  # throw away burn-ins
  colnames(df_theta) <- paste0("tau_k[", seq(ncol(df_theta)), "]")
  
  # (2) Tidy up the "tau" and "sigma_tau" outputut objects
  df_hyperparm <- output[["m"]] %>%
    cbind.data.frame(output[["tau2"]]) %>%
    slice(c((nburn + 1):nDraws))  # throw away burn-ins
  colnames(df_hyperparm) <- c("tau", "sigma_tau")
  
  ### Return the resulting object
  cbind.data.frame(df_hyperparm, df_theta)
}



###'######################################################################
###'
###' get_posterior_DPmeta()
###' 
###' - Tidy up the output object from DPmeta package
###'
###'

get_posterior_DPmeta <- function(output = outp, nburn_DP = 4000){
  
  # (1) Tidy up the "theta" outpututut object
  df_theta <- output$save.state$randsave %>%
    data.frame() %>% 
    dplyr::select(-Prediction) %>%
    slice(c((nburn_DP/2 + 1):nburn_DP))  # throw away burn-ins
  colnames(df_theta) <- paste0("tau_k[", seq(ncol(df_theta)), "]")
  
  # (2) Tidy up the G0 parameters (m, s2), alpha0, and N of clusters outputs
  df_hyperparm <- output$save.state$thetasave %>%
    data.frame() %>%
    dplyr::select(-Y) %>%
    rename(G0_mu = mu, G0_s2 = sigma2, Ncluster = ncluster, alpha0 = alpha) %>%
    dplyr::select(G0_mu, G0_s2, alpha0, Ncluster) %>%
    slice(c((nburn_DP/2 + 1):nburn_DP))  # throw away burn-ins
  
  # Return the resulting object
  cbind.data.frame(df_hyperparm, df_theta)
}



###'######################################################################
###'
###' SEL()
###' 
###' - Calculate squared-error loss (SEL)
###'
###'

SEL <- function(est_vec, true_vec){
  
  K_max <- length(true_vec)
  
  true_k <- est_vec/(K_max + 1)
  est_k <- true_vec/(K_max + 1)
  
  squared_error <- (est_k - true_k)^2
  SEL <- sum(squared_error)*(1/K_max)
  return(SEL)
}



###'######################################################################
###'
###' get_theta_SEL()
###'
###' - Obtain Bayes risk (SEL and SSEL) for 
###'   all individual site-specific effects estimates 
###'   at each simulation 
###' 
###' 

get_theta_SEL <- function(df_est, 
                          est_vec = c("PM", "CB", "GR", "ML"), 
                          weight = rep(1, K_max)){
  
  ### Extract the numbers of sites and posterior summary methods
  K_max <- length(df_est$K)
  J_max <- length(est_vec)
  
  
  ### Get true and estimated theta's
  theta_true <- df_est$theta_true 
  est_vars <- paste0("theta_", str_to_lower(est_vec)) 
  theta_est <- df_est %>% dplyr::select(est_vars) %>%
    as.matrix()  # Estimated theta's
  
  
  ### Generate zero vectors to save (W)SEL and (W)SSEL 
  SEL <- WSEL <- rep(0, K_max*J_max)  # (Weighted) Squared Error Loos
  SSEL <- WSSEL <- rep(0, J_max)      # (Weighted) Sum of SEL
  
  
  ###' Obtain Bayes risk (SEL and SSEL) 
  temp <- theta.ssel(K_max, J_max, 
                     theta_true, theta_est, 
                     weight, 
                     SEL, SSEL, WSEL, WSSEL)
  
  
  ### Return a list of (W)SEL and (W)SSEL
  df_SEL <- cbind.data.frame(temp$sel, temp$wsel)
  df_SSEL <- cbind.data.frame(temp$ssel, temp$wssel)
  
  names(df_SEL) <- c("SEL", "WSEL")
  names(df_SSEL) <- c("SSEL", "WSSEL")
  
  return(list(df_SEL, df_SSEL))
  
}



###'######################################################################
###'
###' get_EDF_ISEL()
###'
###' - Obtain Integrated Squared Error Loss (ISEL) for 
###'   the empirical distribution function (EDF)
###' 
###' 

get_EDF_ISEL <- function(df_est, 
                         est_vec = c("PM", "CB", "GR", "ML")){
  
  ### Extract the numbers of sites and posterior summary methods
  K_max <- length(df_est$K)
  J_max <- length(est_vec)
  
  
  ### Get true and estimated theta's
  theta_true <- df_est$theta_true 
  est_vars <- paste0("theta_", str_to_lower(est_vec)) 
  theta_est <- df_est %>% dplyr::select(est_vars) %>%
    as.matrix()  # Estimated theta's
  
  
  ### Compute ISEL for the selected estimates
  list_temp <- list()
  
  for (m in seq_along(est_vars)){
    
    theta_est_select <- theta_est[, est_vars[m]]
    est_ISEL <- ensrisk(K_max, theta_true, theta_est_select, 0)$isel
    list_temp[[m]] <- est_ISEL
    
  }

  names(list_temp) <- paste0("ISEL_", est_vec)
  
  return(unlist(list_temp))
}



###'######################################################################
###'
###' get_tail_quantiles
###'
###' - Obtain tail quantiles
###'
###'

get_tail_quantiles <- function(df_est, 
                               tail_p, 
                               priquan, 
                               est_vec = c("PM", "CB", "GR", "ML")){
  
  ###' Extract the number of quantiles, and the number of estimators
  ###' Then, initialize an empty vector
  nquan <- length(tail_p)
  J_max <- length(est_vec)
  tail <- rep(0, nquan*J_max)
  K_max <- length(df_est$K)
  
  
  ### Extract the theoretical CDF 
  priquan_vec <- priquan$CDF
  
  
  ### Get true and estimated theta's
  theta_true <- df_est$theta_true 
  est_vars <- paste0("theta_", str_to_lower(est_vec)) 
  theta_est <- df_est %>% dplyr::select(est_vars) %>%
    as.matrix()  # Estimated theta's
  
  
  ### Obtain tail quantiles
  for (p in 1:J_max) {
    for (q in 1: nquan) {
      for (r in 1:K_max) {
        
        if (theta_est[(p - 1)*K_max + r] <= priquan_vec[q])
          tail[(q - 1)*J_max + p] = tail[(q - 1)*J_max + p] + 1
        
      }
    }
  }
  
  for (p in 1:J_max){
    for (q in 1:nquan){
      
      tail[(q - 1)*J_max + p] =  tail[(q - 1)*J_max + p] / K_max
      
    }
  }
  
  ### Return the resulting vector
  tail
}



###'######################################################################
###'
###' get_histogram_bars()
###'
###' - Obtain histogram bars across possible values
###' 
###' 

get_histogram_bars <- function(df_est, 
                               est_vec = c("PM", "CB", "GR", "ML"), 
                               n_bars = 50,  
                               alowhis = -6, ahighhis = 6){
  
  ### Extract the numbers of sites
  K_max <- length(df_est$K)
  
  
  ### Generate a zero-filled vector
  vec_his <- rep(0, n_bars)
  
  
  ### Get estimated theta's
  theta_true <- df_est$theta_true 
  est_vars <- paste0("theta_", str_to_lower(est_vec)) 
  theta_est <- df_est %>% dplyr::select(est_vars) %>%
    as.matrix()  # Estimated theta's
  
  
  ### Compute histogram bars for the selected estimates
  list_temp <- list()
  
  for (m in seq_along(est_vars)){
    
    theta_est_select <- theta_est[, est_vars[m]]
    est_his <- enshis(n_bars, K_max, alowhis, ahighhis, 
                      vec_his, theta_est_select)$his
    list_temp[[m]] <- est_his %>% data.frame()
    
  }
  
  names(list_temp) <- paste0("hisbars_", est_vec)
  

  ### Generate the resulting matrix
  df_temp <- bind_cols(list_temp) %>%
    dplyr::select(Var1, contains("Freq"))
  
  rownames(df_temp) <- df_temp$Var1
  df_temp <- df_temp %>% dplyr::select(-Var1)
  names(df_temp) <- est_vars
  mat_temp <- df_temp %>% as.matrix()
  
  return(mat_temp)
}



###'######################################################################
###'
###' folder_select()
###' 
###' - Conditionally select a subgroup of folders
###'
###'

folder_select <- function(df_select = NULL, 
                          folder_dir = NULL, 
                          leading_zeros = FALSE){
  
  ### Set a folder containing directory
  setwd(folder_dir)
  
  
  ### List all folders within the directory
  vec_folders_all <- list.files()
  
  
  ### Remove leading numbers (if any)
  vec_folders_sub <- str_sub(vec_folders_all, 
                             start = if_else(leading_zeros == TRUE, 5, 0))
  
  
  ### Extract folder names (fname) from the subsetted sim_parm dataframe
  vec_selected <- df_select$fname
  
  
  ### Get a logic vector (T/F) 
  indexl <- vec_folders_sub %in% vec_selected 
  
  
  ### Return a vector of selected folders
  vec_folders <- vec_folders_all[indexl]
  
  return(vec_folders)
}



###'######################################################################
###'
###' collect_selected_files()
###' 
###' - (1) Read in the selected .csv result file 
###'   (2) Attach simulation parameters
###'   (3) Loop over the selected folder list
###'   (4) Collect as a data frame format
###'
###'

collect_selected_files <- function(folder_dir, 
                                   folder_list, 
                                   file_name, 
                                   df_select, 
                                   leading_zeros = FALSE){
  
  
  ### Generate an empty list to collect the files
  list_temp <- list()
  
  
  ### Implement a loop to collect files over folder lists
  
  for (i in seq_along(folder_list)){
    
    ### Extract a folder name in the list
    folder_name <- folder_list[i] 
    folder_name_sub <- str_sub(folder_name, 
                               start = if_else(leading_zeros == TRUE, 5, 0))
    
    
    ### Extract simulation parameters
    sim_parm_sub <- df_select %>%
      filter(fname == folder_name_sub)
    
    
    ### Load the selected .csv file
    file_path <- file.path(folder_dir,
                           folder_name, 
                           file_name)
    
    df_temp <- read.csv(file = file_path)
    
    
    ### Attach simulation parameters
    df_result <- df_temp %>%
      mutate(truth = as.character(sim_parm_sub$truth), 
             assumed = as.character(sim_parm_sub$assumed), 
             ICC = as.numeric(sim_parm_sub$ICC), 
             rsr = as.numeric(sim_parm_sub$rsr), 
             N = as.numeric(sim_parm_sub$N)) %>%
      dplyr::select(truth, assumed, N, ICC, rsr, everything())
    
    ### Embed in the list
    list_temp[[i]] <- df_result
    
  }
  
  ### Convert the list into dataframe format
  df_collect <- bind_rows(list_temp)
  
  return(df_collect)
}



###'######################################################################
###'
###' Gen_True_his()
###'
###' - Define a function to generate the true EDF
###'   conforming to the true distribution
###'
###'

Gen_True_his <- function(gen_name,
                         N = 2000000, 
                         N_bins = 50, 
                         alow = -6, ahigh = 6){
  
  if (gen_name == "Gaussian"){ 
    
    mn <- 0; tau <- sqrt(1) 
    true_theta <- rnorm(N, mn, tau) 
    
  } else if (gen_name == "T"){   
    
    mn <- 0; nu <- 5
    true_theta <- rt(N, nu)*sqrt((nu - 2)/nu)
    
  } else if (gen_name == "ALD"){
    
    mean <- 0; var <- 1; p <- 0.1   
    scale <- sqrt((2*p^2*var)/(1 + p^4))
    location <- mean - ((scale*(1/p - p))/sqrt(2))
    true_theta <- LaplacesDemon::ralaplace(N, location, scale, p)
    
  } else if (gen_name == "Bimodal"){
    
    delta <- 4; eps <- 0.3; ups <- 1
    a <- sqrt((1 - eps) + eps*ups^2 + eps*(1 - eps)*delta^2)
    ind <- runif(N) < (1 - eps)
    true_theta <- ind*rnorm(N, -eps*delta/a, sqrt(1/a^2)) + 
      (1 - ind)*rnorm(N, (1 - eps)*delta/a, sqrt(ups^2/a^2))
    
  } else if (gen_name == "Skew"){
    
    mean <- 0; var <- 1; slant <- 10
    delta <- slant/sqrt(1 + slant^2)
    scale <- sqrt(var/(1-(2*delta^2/pi)))
    location <- mean - scale*sqrt(2/pi)*delta
    true_theta <- sn::rsn(n = N, xi = location, omega = scale, 
                          alpha = slant)[1:N]
    
  } else if (gen_name == "Mixed"){
    
    delta <- 5; eps <- 0.3; ups <- 2
    a <- sqrt((1 - eps) + eps*ups^2 + eps*(1 - eps)*delta^2)
    ind <- runif(N) < (1 - eps)
    true_theta <- ind*rnorm(N, -eps*delta/a, sqrt(1/a^2)) + 
      (1 - ind)*rnorm(N, (1 - eps)*delta/a, sqrt(ups^2/a^2))
    
  }
  
  ### Calculate EDF estimates for each bin 
  vec_his <- rep(0, N_bins)
  est_his <- enshis(N_bins, N, alow, ahigh, vec_his, true_theta)$his
  true_his <- est_his/N 
  
  ### Returen the resulting dataframe 
  true_his %>% as.data.frame() %>%
    rename(Bin = Var1, true_his = Freq) %>%
    mutate(truth = gen_name) %>%
    select(truth, Bin, true_his)
}



###'######################################################################
###'
###' get_lincom_results()
###' 
###' - Define a function to return lincom results
###' 
###'

get_lincom_results <- function(fit, df_sub, 
                               var1, var2, 
                               ref1, ref2, 
                               vec_ref){
  
  ###' Generate a dataframe for regression terms
  ###' Tag reference categories
  ###' Generate hypotheses
  var1_levels <- unique(df_sub[var1])[[var1]]
  var2_levels <- unique(df_sub[var2])[[var2]]
  
  df_vars <- expand.grid(var1_levels, var2_levels) %>% 
    data.frame() %>% set_names(c("ovar1", "ovar2")) %>%
    mutate(evar1 = paste0(var1, ovar1), 
           evar2 = paste0(var2, ovar2)) %>%
    unite(int_terms, evar1, evar2, sep = ":", remove = FALSE) %>%
    dplyr::select(ovar1, ovar2, evar1, evar2, int_terms) %>%
    mutate(tag1 = (ovar1 == ref1) %>% as.numeric(), 
           tag2 = (ovar2 == ref2) %>% as.numeric(), 
           tag_mod = if_else(tag1 == 0 & tag2 == 0, 1, 0), 
           elem1 = if_else(tag1 == 1, "", as.character(evar1)), 
           elem2 = if_else(tag2 == 1, "", as.character(evar2)), 
           elem3 = if_else(tag_mod == 0, "", as.character(int_terms)), 
           hypo = case_when(
             tag_mod == 1 ~ paste0(elem1, " + ", elem3, " = 0"), 
             tag1 == 0 & tag2 == 1 ~ paste0(elem1, " = 0"),
             tag1 == 1 & tag2 == 0 ~ paste0(elem2, " = 0"), 
             tag1 == 1 & tag2 == 1 ~ NA_character_, 
             TRUE ~ NA_character_)) %>%
    drop_na()
  
  
  ### Get lincom test results
  obj_ht <- glht(fit, linfct = df_vars$hypo)
  
  df_glht <- tidy(summary(obj_ht)) %>%
    full_join_track(tidy(confint(obj_ht)), by = c("lhs", "rhs", "estimate")) 
  
  df_lincom <- cbind.data.frame(df_vars, df_glht) %>%
    dplyr::select(-lhs, -rhs)
  
  
  ### Plot!
  df_plot <- df_lincom %>%
    dplyr::select(ovar1, ovar2, hypo, estimate, conf.low, conf.high) %>%
    filter(ovar1 != ref1)
  
  p <- ggplot(data = df_plot, aes(x = ovar2, y = estimate, 
                                  group = ovar1, color = ovar1, shape = ovar1)) +
    geom_point(aes(y = estimate), size = 3, position = position_dodge(width = 0)) + 
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), 
                  position = position_dodge(width = 0), width = 0.1) + 
    geom_line(aes(y = estimate)) + 
    geom_hline(yintercept = 0, color = "gray30", linetype = "dashed") + 
    # geom_text(aes(label = round(estimate, 0)), color = "black", 
    #           position = position_stack(vjust = 0.1, reverse = TRUE), size = 3) + 
    scale_y_continuous(labels = comma) + 
    scale_color_manual(values = color_palette[seq(unique(df_plot$ovar1))]) + 
    scale_shape_manual(values = shape_palette[seq(unique(df_plot$ovar1))]) + 
    theme_trend +
    labs(subtitle = vec_ref)
  
  ### Return a list object
  list(df_lincom, df_plot, p)
}