
###'######################################################################
###'
###' Category: Define functions
###' 
###' Task: Define functions for data generation
###'
###' Date: 2020-03-27
###' 
###' Author: JoonHo Lee (joonho@berkeley.edu)
###' 
###'

###'######################################################################
###'
###' Call necessary libraries
###'
###'

library(LaplacesDemon)
library(sn)


###'######################################################################
###'
###' Gen_Gaussian()
###' 
###' - Generate data assuming G follows Gaussian
###'
###'

Gen_Gaussian <- function(n, ICC, rsr) {
  
  ### Calculate the geometric mean of level-1 error variances
  rgt <- (1 - ICC)/ICC
  
  ### Generate within-site variances (sds) for K sites
  sigma2max <- rgt*rsr
  sigma2min <- sigma2max/(rsr^2)
  sigma2 <- exp(seq(from = log(sigma2min), to = log(sigma2max), length = n))
  sd <- sqrt(sigma2)
  
  ### Set mean and SD of standard normal distribution   
  mn <- 0; tau <- sqrt(1) 
  
  ### Generate data (Note the difference between theta and Y)
  data_obj <- list()
  data_obj$theta <- rnorm(n, mn, tau)        # site means
  data_obj$Y <- rnorm(n, data_obj$theta, sd) # within-site outcome (only one per site) 
  data_obj$sigma2 <- sigma2                  # site SDs
  return(data_obj)
}



###'######################################################################
###'
###' Gen_T()
###' 
###' - Generate data assuming G follows Student T distribution with df = nu
###'
###'

Gen_T <- function(n, ICC, rsr, nu) {
  
  ### Calculate the geometric mean of level-1 error variances
  rgt <- (1 - ICC)/ICC
  
  ### Generate within-site variances (sds) for K sites
  sigma2max <- rgt*rsr
  sigma2min <- sigma2max/(rsr^2)
  sigma2 <- exp(seq(from = log(sigma2min), to = log(sigma2max), length = n))
  sd <- sqrt(sigma2)
  
  ## Mean of T-distribution
  mn <- 0
  
  ## Generate data object
  data_obj <- list()
  data_obj$theta <- rt(n, nu)*sqrt((nu - 2)/nu)
  data_obj$Y <- rnorm(n, data_obj$theta, sd)
  data_obj$sigma2 <- sigma2
  return(data_obj)
}



###'######################################################################
###'
###' Gen_ALD()
###' 
###' - Generate data assuming G follows asymmetric Laplace distribution
###'
###'

Gen_ALD <- function(n, ICC, rsr, 
                    mean = 0, var = 1, p = 0.1) {
  
  ### Calculate the geometric mean of level-1 error variances
  rgt <- (1 - ICC)/ICC
  
  ### Generate within-site variances (sds) for K sites
  sigma2max <- rgt*rsr
  sigma2min <- sigma2max/(rsr^2)
  sigma2 <- exp(seq(from = log(sigma2min), to = log(sigma2max), length = n))
  sd <- sqrt(sigma2)
  
  ###' Generate location, scale, and skewness parameters of ALD distribution
  ###' to get E(theta) = 0 and Var(theta) = 1
  scale <- sqrt((2*p^2*var)/(1 + p^4))
  location <- mean - ((scale*(1/p - p))/sqrt(2))
  
  ## Generate data object
  data_obj <- list()
  data_obj$theta <- LaplacesDemon::ralaplace(n, location, scale, p)
  data_obj$Y <- rnorm(n, data_obj$theta, sd)
  data_obj$sigma2 <- sigma2
  return(data_obj)
}



###'######################################################################
###'
###' Gen_Mixture()
###' 
###' - Generate data assuming G follows a mixture of two Gaussian distributions
###'
###'  [Additional arguments]
###'  
###'  (1) delta: distance between two mixtures
###'  (2) eps: proportion of the small portion of two mixtures
###'  (3) ups: SD for components => want components to have equal variance
###'  
###'

Gen_Mixture <- function(n, ICC, rsr, 
                        delta = 4, eps = 0.30, ups = 1){
  
  ### Initiate empty list 
  data_obj <- list()
  
  
  ### Calculate the geometric mean of level-1 error variances
  rgt <- (1 - ICC)/ICC
  
  
  ### Define a normalizing factor `a`
  a <- sqrt((1 - eps) + eps*ups^2 + eps*(1 - eps)*delta^2)
  

  ###' Simulate a mixture of 2 normals with mean 0 and var 1  
  ###' with this parameterization
  ind <- runif(n) < (1 - eps)
  
  data_obj$theta <- ind*rnorm(n, -eps*delta/a, sqrt(1/a^2)) + 
    (1 - ind)*rnorm(n, (1 - eps)*delta/a, sqrt(ups^2/a^2))
  
  
  ### Generate uncertainty estimates for K sites
  sigma2max <- rgt*rsr
  sigma2min <- sigma2max/(rsr^2)
  sigma2 <- exp(seq(from = log(sigma2min), to = log(sigma2max), length = n))
  sd <- sqrt(sigma2)
  
  
  ### Generate the data objects (observed Y)
  data_obj$Y <- rnorm(n, data_obj$theta, sd)
  data_obj$sigma2 <- sigma2
  return(data_obj)
  
}



###'######################################################################
###'
###' Gen_SkewN()
###' 
###' - Generate data assuming G follows skewed normal distribution with 
###'   marginal mean and variance of SkewN(0, 1^2)
###'
###'

Gen_SkewN <- function(n, ICC, rsr, 
                      mean = 0, var = 1, slant = 5) {
  
  ### Calculate the geometric mean of level-1 error variances
  rgt <- (1 - ICC)/ICC
  
  ### Generate within-site variances (sds) for K sites
  sigma2max <- rgt*rsr
  sigma2min <- sigma2max/(rsr^2)
  sigma2 <- exp(seq(from = log(sigma2min), to = log(sigma2max), length = n))
  sd <- sqrt(sigma2)
  
  ###' Generate location and scale parameters of skewed normal
  ###' to get E(theta) = 0 and Var(theta) = 1
  delta <- slant/sqrt(1 + slant^2)
  scale <- sqrt(var/(1-(2*delta^2/pi)))
  location <- mean - scale*sqrt(2/pi)*delta
  
  ### Generate data object
  data_obj <- list()
  data_obj$theta <- sn::rsn(n = n, xi = location, omega = scale, 
                            alpha = slant)[1:n]
  data_obj$Y <- rnorm(n, data_obj$theta, sd)
  data_obj$sigma2 <- sigma2
  return(data_obj)
}



###'######################################################################
###'
###' list_to_df_G()
###' 
###' - A helper function to convert the resulting list to dataframe
###'
###'

list_to_df_G <- function(list_G){
  
  list_G %>%
  data.frame() %>%
    rownames_to_column("K") %>%
    mutate_at(.vars = c("K"), .funs = as.numeric)
  
}
