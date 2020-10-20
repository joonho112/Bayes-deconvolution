
###'######################################################################
###'
###' Category: Define functions
###' 
###' Task: Define functions for converting quantiles to the CDF 
###'       (for the theoretical quantiles)
###'
###' Date: 2020-03-27
###' 
###' Author: JoonHo Lee (joonho@berkeley.edu)
###' 
###'

###'######################################################################
###' 
###'  Set quantiles (tail probabilities)
###' 
###'

# tail_p <- c(0.05, 0.10, 0.25, 0.75, 0.90, 0.95)



###'######################################################################
###' 
###' Q2CDF_Gaussian()
###' 
###' - Convert quantiles to the "Gaussian" CDF
###' 
###' 

Q2CDF_Gaussian <- function(tail_p){
  
  Quantile <- tail_p
  CDF <- qnorm(tail_p)
  
  df_temp <- tibble(Quantile, CDF)
  return(df_temp)
}



###'######################################################################
###' 
###' Q2CDF_T()
###' 
###' - Convert quantiles to the "Student T" CDF
###' 
###' 

Q2CDF_T <- function(tail_p, nu = 5){
  
  Quantile <- tail_p
  sd <- sqrt(nu/(nu - 2))
  CDF <- qt(Quantile, nu)*(1/sd) 
  
  df_temp <- tibble(Quantile, CDF)
  return(df_temp)
}



###'######################################################################
###' 
###' Q2CDF_Mixure()
###' 
###' - Convert quantiles to the CDF of "a mixture of two Gaussian distributions"
###' - Why theoretical? Because we simulate the large number (20,000)
###' 

Q2CDF_Mixture <- function(tail_p, n_size = 20000, 
                          delta = 4, eps = 0.3, ups = 1){
  
  ### Get quantiles
  Quantile <- tail_p 
  
  
  ### Calculate PDF (not CDF)
  ind <- runif(n_size) < (1 - eps)
  
  a <- sqrt((1 - eps) + eps*ups^2 + eps*(1 - eps)*delta^2)
  
  PDF <- ind*rnorm(n_size, -eps*delta/a, sqrt(1/a^2)) + 
    (1 - ind)*rnorm(n_size, (1 - eps)*delta/a, sqrt(ups^2/a^2))
  
  
  ### Calculate CDF (quantile() for inverse CDF)
  CDF <- as.numeric(quantile(PDF, tail_p))
  
  df_temp <- tibble(Quantile, CDF)
  return(df_temp)
}



###'######################################################################
###' 
###' Q2CDF_ALD()
###' 
###' - Convert quantiles to the CDF of asymmetric Laplace distribution
###' 
###' 

Q2CDF_ALD <- function(tail_p, n_size = 20000, 
                      mean = 0, var = 1, p = 0.1){
  
  ### Get quantiles
  Quantile <- tail_p 
  
  
  ###' Generate location, scale, and skewness parameters of ALD distribution
  ###' to get E(theta) = 0 and Var(theta) = 1
  scale <- sqrt((2*p^2*var)/(1 + p^4))
  location <- mean - ((scale*(1/p - p))/sqrt(2))
  
  
  ### Calculate PDF (not CDF)
  PDF <- LaplacesDemon::ralaplace(n_size, location, scale, p)
  
  
  ### Calculate CDF (quantile() for inverse CDF)
  CDF <- as.numeric(quantile(PDF, tail_p))
  
  df_temp <- tibble(Quantile, CDF)
  return(df_temp)
}


###'######################################################################
###' 
###' Q2CDF_SkewN()
###' 
###' - Convert quantiles to the CDF of "Skewed Normal" distribution
###' 
###' 

Q2CDF_SkewN <- function(tail_p, n_size = 20000, 
                        mean = 0, var = 1, slant = 5){
  
  ### Get quantiles
  Quantile <- tail_p 
  
  
  ###' Generate location and scale parameters of skewed normal
  ###' to get E(Y) = 0 and Var(Y) = 1
  delta <- slant/sqrt(1 + slant^2)
  scale <- sqrt(var/(1-(2*delta^2/pi)))
  location <- mean - scale*sqrt(2/pi)*delta
  
  
  ### Calculate PDF (not CDF)
  PDF <- sn::rsn(n = n_size, xi = location, omega = scale, 
                 alpha = slant)[1:n_size]
  
  
  ### Calculate CDF (quantile() for inverse CDF)
  CDF <- as.numeric(quantile(PDF, tail_p))
  
  df_temp <- tibble(Quantile, CDF)
  return(df_temp)
}

