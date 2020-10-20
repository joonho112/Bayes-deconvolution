
###'######################################################################
###'
###' Category: Code Experiment
###' 
###' Task: Hierarchical model simulations  
###' 
###'       Paddock, S. M., Ridgeway, G., Lin, R., & Louis, T. A. (2006). 
###'       Flexible distributions for triple-goal estimates in 
###'       two-stage hierarchical models. 
###'       Computational statistics & data analysis, 50(11), 3243-3262.
###'       
###' Data: Simulated data
###' 
###' Data: 2020-02-20
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
work_dir <- c("~/Treatment-effect-heterogeneity/Multisite Trials")
setwd(work_dir)


### Set a data directory
data_dir <- file.path(work_dir, "datasets")


### Call libraries
library(tidyverse)
library(readxl)
library(Hmisc)
library(foreign)
library(haven)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)



###'######################################################################
###'
###' Install and load hhsim R package => Working!
###'
###'

# install.packages("remotes")
# 
# remotes::install_github("gregridgeway/hhsim")

library(hhsim)



###'######################################################################
###'
###' Prepare simulation
###'
###' (1) rgt
###' - Geometric mean of the within-site variance
###' - The informativeness of the data
###' - Large within-site variance => small information for the site-specific effect
###' 
###' (2) rsr
###' - the heterogeneity of within-site variances across sites
###'   
###' (3) n
###' - Number of sites (K)
###' 
###' (4) df
###' - df of T-distributions 
###'
###'

### Set seed 
set.seed(11172)


### Generate a table of the simulation parameters
simulparam <- expand.grid(truth = c("Gaussian", "Mdp", "Tdist"),
                          assumed = c("Gaussian", "Mdp", "Tdist", "NPML", "SBR"),
                          rgt = c(0.1, 0.33, 1, 10),
                          rsr = c(1, 5, 10),
                          n = c(25, 50, 100, 200),
                          df = 5)


### Reorder the table
i <- with(simulparam, order(n, rgt, rsr, truth, assumed))

simulparam <- simulparam[i, ]


### Add a variable 
simulparam$fname <- rep("", nrow(simulparam))

for(j in 1:nrow(simulparam)){
  
  simulparam$fname[j] <- with(simulparam[j,], 
                              paste(paste(truth, assumed, "_", 
                                          rgt, "_", rsr, "_", 
                                          n, sep = "")))
  
}



###'######################################################################
###'
###' Define functions to generate simulated data
###' 
###' => Specifying the population distribution of G
###'
###'

### (1) Gaussian distribution

GenerateGaussian <- function(n, rgt, rsr) {
  
  ##' Maximum within-site variance
  ##' geometric mean of sigma * sigma^2_largest/sigma^2_smallest
  sigma2max = rgt*rsr
  
  ## Minimum within-site variance
  sigma2min = sigma2max/(rsr^2)
  
  ## Generate within-site variances (sds) for each n = 100 site
  sigma2 = exp(seq(from = log(sigma2min), to = log(sigma2max), length = n))
  sd = sqrt(sigma2)
  
  
  ## Mean and SD of standard normal distribution   
  tau = sqrt(1)      # sd
  mn = 0             # mean
  
  ##' Generate data
  ##' Note the difference between theta and Y
  data.obj = list()
  data.obj$theta = rnorm(n, mn, tau)        # site means
  data.obj$Y = rnorm(n, data.obj$theta, sd) # within-site outcome (only one) 
  data.obj$sigma2 = sigma2                  # site SDs
  return(data.obj)

}


### (2) T-distribution with df = 5 (nu)

GenerateT <- function(n, rgt, rsr, nu) {
  
  ##' Maximum within-site variance
  ##' geometric mean of sigma * sigma^2_largest/sigma^2_smallest
  sigma2max = rgt*rsr
  
  ## Minimum within-site variance
  sigma2min = sigma2max/(rsr^2)
  
  ## Generate within-site variances (sds) for each n = 100 site
  sigma2 = exp(seq(from = log(sigma2min), to = log(sigma2max), length = n))
  sd = sqrt(sigma2)
  
  ## Mean of T-distribution
  mn = 0
  
  ## Generate data object
  data.obj = list()
  data.obj$theta = rt(n, nu)*sqrt((nu - 2)/nu)
  data.obj$Y = rnorm(n, data.obj$theta, sd)
  data.obj$sigma2 = sigma2
  return(data.obj)
  
}


###' (3) Dirichlet process 1 (DP-1): A mixture of two T-distributions
###'     Generate 75% from a t(8) and 25% from a t(50)

GenerateSmix <- function(n, rgt, rsr) {
  
  ##' Maximum within-site variance
  ##' geometric mean of sigma * sigma^2_largest/sigma^2_smallest
  sigma2max = rgt*rsr
  
  ## Minimum within-site variance
  sigma2min = sigma2max/(rsr^2)
  
  ## Generate within-site variances (sds) for each n = 100 site
  sigma2 = exp(seq(from = log(sigma2min), to = log(sigma2max), length = n))
  sd = sqrt(sigma2)
  
  ##' Generate data
  ##' 75% from a t(8) and 25% from a t(50)
  data.obj = list()
  nnorm = floor(0.25*n)
  data.obj$theta = rt(n, 8)
  data.obj$theta[1:nnorm] = rt(nnorm, 50)
  data.obj$Y = rnorm(n, data.obj$theta, sd)
  data.obj$sigma2 =  sigma2
  return(data.obj)

}


### (4) Dirichlet process 2: smoother nonparametric G

GenerateMdp = function(n, rgt, rsr){
  
  ## Set parameters 
  data.obj = list()
  data.obj = NULL
  delta = 4   # distance between two mixtures
  ups = 1     # want components to have equal variance
  eps = 0.2   # proportion of the small portion of two mixtures
  
  
  ##' Simulate a mixture of 2 normals with mean 0 and var 1 
  ##' with this parameterization
  a = sqrt((1 - eps) + eps*ups^2 + eps*(1 - eps)*delta^2)
  
  ind = runif(n) < (1 - eps)
  
  data.obj$theta = ind*rnorm(n, -eps*delta/a, sqrt(1/a^2)) + 
                   (1 - ind)*rnorm(n, (1 - eps)*delta/a, sqrt(ups^2/a^2))
  
  
  ## Simulate true distribution for quantiles etc.
  ind = runif(n) < (1 - eps)
  
  data.obj$theta = ind*rnorm(n, -eps*delta/a, sqrt(1/a^2)) + 
                   (1 - ind)*rnorm(n, (1 - eps)*delta/a, sqrt(ups^2/a^2))
  
  tailp = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95)
  
  sigma2max = rgt*rsr
  sigma2min = sigma2max/(rsr^2)
  sigma2 = exp(seq(from = log(sigma2min), to = log(sigma2max), length = n))
  sd = sqrt(sigma2)
  
  data.obj$Y = rnorm(n, data.obj$theta, sd)
  data.obj$sigma2 = sigma2
  #   ind = runif(1000) < (1-eps)
  #   quant=ind*rnorm(1000,-eps*delta/a,sqrt(1/a^2)) + (1-ind)*rnorm(1000,(1-eps)*delta/a,sqrt(ups^2/a^2))
  data.obj$priquan = mdppriquan  # "True" CDF
  
  return(data.obj)

}



###'######################################################################
###'
###' Get true CDF of thetas from mixture, via simulation 
###' 
###' These pararmeters will change depending on 
###' 
###' what changes in GenerateMdp function
###' 
###' => Why TRUE? 
###'    because we simulate the very large number of cases (20,000)
###'
###'

### Set parameters
tailp = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95)
delta = 4   # distance between two mixtures
ups = 1     # want components to have equal variance
eps = 0.2   # proportion of the small portion of two mixtures


### Get true CDF of thetas from mixture
a = sqrt((1 - eps) + eps*ups^2 + eps*(1 - eps)*delta^2)

#ind = runif(2000000) < (1-eps)
#quant=ind*rnorm(2000000,-eps*delta/a,sqrt(1/a^2)) + (1-ind)*rnorm(2000000,(1-eps)*delta/a,sqrt(ups^2/a^2))

ind = runif(20000) < (1 - eps)

quant = ind*rnorm(20000, -eps*delta/a, sqrt(1/a^2)) + 
        (1 - ind)*rnorm(20000, (1 - eps)*delta/a, sqrt(ups^2/a^2))

mdppriquan = quantile(quant, tailp)



###'######################################################################
###'
###' Start simulation loops
###'
###'

LL = length(simulparam$truth)  # total 180 cells


for(LLind in 1:LL) {  ### Start of loop over LLind
  
  
  ###'######################################################################
  ###'
  ###' Initial values and required parameters for simulation
  ###'
  ###'

  ### Generate an empty vector to contain variance of true thetas
  evartheta = NULL
  
  
  ### Set numbers of iterations
  nsimu = 500    # number of Monte Carlo iterations to perform
  nmcmc = 500    # number of MCMC iterations per each model fit
  nburn = 100    # burn-in iterations to discard at each round
  
  nmcmcall = nmcmc + nburn
  
  nmcm = as.vector(as.integer(nmcmc)) # need for R to C conversion
  
  
  ### Sample size to generate
  kmax = simulparam$n[LLind] 
  
  
  ### Set DP quantile parameters
  alow = -6
  ahigh = 6
  alowhis = -6
  ahighhis = 6
  nquan = 6
  tailp = c(0.05, 0.1, 0.25, 0.75, 0.9, 0.95)
  
  
  ### Number of estimators to consider (bayes, gr, mle)
  jmax = 3   
  
  
  ### Weight (not used yet)
  weight = rep(1, kmax)  
  
  
  ### Number of bars in average histogram
  nhis = 50 
  
  
  ### Rank estimates
  simurankest = rep(0, jmax)
  
  
  
  ###'######################################################################
  ###'
  ###' Initialize key vectors before simulation study
  ###'
  ###'
  
  ### Site index
  kindex = c(0:(kmax - 1))
  
  
  ### Site index quantiles
  p = (2*(kindex + 1) - 1)/(2*kmax)
  
  
  ### A vector for rank??
  irstar = rep(0, kmax)
  
  
  ### Zero vectors for the simulated parameters for each estimator
  simumean = simuvar = simuisel = simussel = simuwsse = rep(0, jmax)
  
  simusel = rep(0, jmax*kmax)
  
  simuwsel = rep(0, jmax*kmax)
  
  simuselrank = 0
  
  simutail = rep(0, nquan*jmax)
  
  simuhis = rep(0, nhis*jmax)
  
  avgwsel = avgsel = rep(0, jmax*kmax)
  
  avgmean = avgvar = avgisel = weff = eff = avgwssel = avgssel = rep(0, jmax)
  
  avghis = rep(0, nhis*jmax)
  
  avgtail = rep(0, nquan*jmax)
  
  niter = 50
  
  
  
  ###'######################################################################
  ###'
  ###' Start Monte Carlo iteration
  ###'
  ###'
  
  ### Extract True data generating model name and assumed model name
  genname = paste("Generate", simulparam$truth[LLind], sep = "")
  anlname = paste("Gaussian", simulparam$assumed[LLind], sep="")
  
  
  ### Start Monte Carlo loop
  for (nmccount in 1:nsimu){
    
    ### Print progress
    cat(genname, anlname, nmccount, "\n")
    
    
    
    ###'######################################################################
    ###'
    ###' Generate data: True data generating model
    ###'
    ###'
    
    if(genname == "GenerateGaussian") {
      
      ## Convert quantiles to CDF
      priquan = qnorm(tailp)
      
      ## Generate data based on normal distribution 
      data.obj = GenerateGaussian(kmax,
                                  simulparam$rgt[LLind],
                                  simulparam$rsr[LLind])
      
    } else if (genname == "GenerateTdist") {
      
      ## Convert quantiles to CDF
      priquan =  qt(tailp, simulparam$df[LLind])*
        sqrt((simulparam$df[LLind] - 2)/simulparam$df[LLind])
      
      ## Generate data based on T distribution
      data.obj = GenerateT(kmax,
                           simulparam$rgt[LLind],
                           simulparam$rsr[LLind],
                           simulparam$df[LLind])
    
    } else if (genname == "GenerateMdp") {
      
      ## Generate data based on two normal mixtures
      data.obj = GenerateMdp(kmax,
                             simulparam$rgt[LLind],
                             simulparam$rsr[LLind])
      
      ## Get CDF for each quantile
      priquan = data.obj$priquan
      
    }
    
    
    ### Get the variance of true theta
    evartheta = c(evartheta, c(var(data.obj$theta)))
    
    
    
    ###'######################################################################
    ###'
    ###' Fit the assumed model & get posterior sample 
    ###'
    ###' 
    
    if(anlname == "GaussianGaussian"){
      
      outp = GaussianGaussian(data.obj$Y,
                              data.obj$sigma2,
                              nmcmcall)
      
    } else if (anlname == "GaussianTdist"){
      
      outp = GaussianT(data.obj$Y,
                       data.obj$sigma2,
                       nDraws = nmcmcall,
                       nu = simulparam$df[LLind])
      
    } else if (anlname == "GaussianMdp") {
      
      outp = GaussianMDP(data.obj$Y,
                         data.obj$sigma2,
                         nDraws = nmcmcall,
                         numconfig = 3,
                         m = 0,      
                         # initial value of m
                         #  w = 1,
                         #  W = 10,
                         #  par1 = 10,
                         #  par2 = 0.1)
                         w = 1,
                         W = 1,
                         par1 = 4,
                         par2 = 4)
      
    } else if(anlname =="GaussianNPML") {
      
      EBestimate = GaussianNPML(data.obj$Y,
                                data.obj$sigma2)
      
      outp = GaussianEBSimulate(data.obj$Y,
                                data.obj$sigma2,
                                mu = EBestimate$mu,
                                alpha = EBestimate$alpha,
                                size = nmcmcall)
      
    } else if(anlname == "GaussianSBR") {
      
      EBestimate = GaussianSBR(data.obj$Y,
                               data.obj$sigma2,
                               n.iters = round(3*log(kmax)))
      
      outp = GaussianEBSimulate(data.obj$Y,
                                data.obj$sigma2,
                                mu = EBestimate$mu,
                                alpha = EBestimate$alpha,
                                size = nmcmcall)
    }
    
    # for preferring the Gaussian (Kmax=100, sigma2=1): w=1,W=10,par1=10,par2=0.1
    # to prefer the MDP / mxiture: w=1,W=1,par1=4,par2=4
    
    
    
    ###'######################################################################
    ###'
    ###' Calculate posterior means from the estimated "outp"
    ###'
    ###'
    
    postmean = data.matrix(data.frame(outp$theta))
    
    if(is.element(anlname, c("GaussianGaussian", "GaussianTdist", "GaussianMdp"))){
      
      postmean = postmean[ , c((nburn + 1):nmcmcall)]
      pmmean = rowMeans(postmean)
      
    } else {
      
      pmmean = EBestimate$theta
      postmean = t(postmean[c((nburn + 1):nmcmcall), ])
    
    }
    
    
    ###'######################################################################
    ###'
    ###' Compute GR esimate (Triple-goal estimator)
    ###' 
    ###' => Get rank estimates and their squared error losses as well
    ###'
    ###'
    
    ### Generate zero vectors to save estimates
    grest = ensemble = probens = rep(0, kmax)
    
    
    ###' Get triple-goal estimates (using the gr.quantile() function)
    ###' Get quantiles - for ensemble estimation  
    
    tmp = gr.quantile(kmax, 
                      alow, ahigh, 
                      p, 
                      postmean, ensemble, probens,
                      nmcmc, niter)
    
    ensemble = tmp$ensemble
    
    
    ###' Get rank estimates
    ###' Rank estimates are not necessarily integers
    
    rhat = rep(0, kmax)
    
    for (k in 1:kmax) { 
      
      for(j in 1:kmax) {
        
        tmp = ifelse(j != k, 
                     sum(postmean[k, ] >= postmean[j, ])/length(postmean[k, ]), 
                     1)
        
        rhat[k] = rhat[k] + tmp
        
      }
    }
    
    irstar <- rank(rhat)
    
    
    ### Order GR estimates based on the estimated ranks
    grest <- ensemble[irstar]
    
    
    ### Extract true theta rank
    thetatrue = data.obj$theta
    ranktrue  = rank(thetatrue)
    
    
    ### Calculate squared error losses for the estimated ranks
    selrank    = sum((irstar/(kmax + 1) - ranktrue/(kmax + 1))^2)/kmax
    selpmrank  = sum((rank(pmmean)/(kmax + 1) - ranktrue/(kmax + 1))^2)/kmax
    selgrrank  = sum((rank(grest)/(kmax + 1) - ranktrue/(kmax + 1))^2)/kmax
    selmlrank  = sum((rank(data.obj$Y)/(kmax + 1) - ranktrue/(kmax + 1))^2)/kmax
    
    
    ### Collect estimates (1. posterior means, 2. GR estimates, 3. raw estimates)
    est = cbind(pmmean, grest, data.obj$Y)
    
    rankest = c(selpmrank, selgrrank, selmlrank)
    
    simurankest = simurankest + rankest
    
    
    
    ###'######################################################################
    ###'
    ###' Evaluate estimates with loss functions (squared error loss, SEL)
    ###' 
    ###' at this iteration of sim study
    ###' 
    ###' => Bayes risk (SEL and SSEL) for all estimates at each simulation 
    ###'
    ###'
    
    ### Prepare zero vectors to save SEL and SSEL estimates
    sel = wsel = rep(0, kmax*jmax)
    ssel = wssel = rep(0, jmax)
    
    
    
    ### SSEL for individual estimation (using the "theta.ssel" function)
    
    weight = rep(1, kmax)

    
    tmp = theta.ssel(kmax, jmax, 
                     thetatrue, est, 
                     weight, 
                     sel, ssel, 
                     wsel, wssel)
    
    sel = tmp$sel     # SEL
    ssel = tmp$ssel   # SSEL
    wsel = tmp$wsel   # Weighted SEL
    wssel = tmp$wssel # Weighted SSEL
    
    
    
    ### Get SEL and SSEL over simulations 
    tmp =  theta.simu(jmax, kmax, 
                      sel, ssel, 
                      simusel, simussel,
                      wsel, wssel, 
                      simuwsel, simuwsse)
    
    simusel = tmp$simusel
    simuwsel = tmp$simuwsel
    simussel = tmp$simussel
    simuwsse = tmp$simuwsse
    
    
    ###' Get estimated means and variances for each three estimator
    ###' 1. posterior means, 2. GR estimates, 3. raw estimates
    estvar = estmean = rep(0, jmax)
    
    tmp = ensmom(kmax, jmax, 
                 est, 
                 estmean, estvar)
    
    estmean = tmp$estmean
    
    estvar  = tmp$estvar
    
    
    ###' Get ISEL (Ensemble estimation) for each three estimator
    ###' 1. posterior means, 2. GR estimates, 3. raw estimates
    
    grisel = mlisel = pmisel = 0
    
    pmisel  = ensrisk(kmax, 
                      thetatrue,
                      pmmean, pmisel)$isel
    
    grisel  = ensrisk(kmax,
                      thetatrue,
                      grest, grisel)$isel
    
    mlisel  = ensrisk(kmax,
                      thetatrue,
                      data.obj$Y, mlisel)$isel
    
    
    ### Get histogram bars    
    
    grhis = mlhis = pmhis = rep(0, nhis)
    
    pmhis  = enshis(nhis, 
                    kmax, 
                    alowhis, ahighhis,
                    pmhis, pmmean)$his
    
    mlhis  = enshis(nhis, 
                    kmax, 
                    alowhis, ahighhis, 
                    mlhis, data.obj$Y)$his
    
    grhis  = enshis(nhis, 
                    kmax, 
                    alowhis, ahighhis, 
                    grhis, grest)$his
    
    
    ### Get tail quantile estimates
    
    tail = rep(0, nquan*jmax)
    
    for (j in 1:jmax) {
      for (i in 1: nquan) {
        for (k in 1:kmax) {
          
          if (est[(j - 1)*kmax+k] <= priquan[i])
            tail[(i - 1)*jmax + j] = tail[(i - 1)*jmax + j] + 1
          
        }
      }
    }
    
    for (j in 1:jmax){
      for (i in 1:nquan){
        
        tail[(i - 1)*jmax + j] =  tail[(i - 1)*jmax + j] / kmax
      
      }
    }

    
    ### Collect all ISEL estimates
    
    allisel = c(pmisel, grisel, mlisel)
    
    
    ### Get Bayes risk (SEL and SSEL) for all estimates at this simulation 
    
    simuisel = simuisel + allisel
    
    simumean = simumean + estmean
    
    simuvar = simuvar + estvar
    
    simuhis = simuhis + cbind(pmhis, grhis, mlhis)
    
    simutail = simutail + tail
    
    simuselrank = simuselrank + selrank
    
    
    
    ###'######################################################################
    ###'
    ###' End of loop over Monte Carlo iteration (nmccount)
    ###'
    ###'
    
  }  ### End of loop over Monte Carlo iteration (nmccount)
  
  
  
  ###'######################################################################
  ###'
  ###' Average over all simulations
  ###'
  ###'
  
  ### SEL and WSEL
  avgsel = simusel/nsimu
  avgwsel = simuwsel/nsimu
  
  
  ### SSEL and WSSEL
  avgssel = simussel / nsimu
  avgwssel = simuwsse / nsimu
  
  
  ### ISEL
  avgisel = simuisel / nsimu
  
  
  ### Posterior means and variances
  avgmean = simumean / nsimu
  avgvar = simuvar / nsimu
  
  
  ### Histogram bars and tail quantiles
  avghis = simuhis/(nsimu*kmax)
  avgtail = simutail/nsimu
  
  
  ### Rank and Rank estimate
  avgselrank = simuselrank / nsimu
  avgrankest = simurankest / nsimu
  
  
  
  ###'######################################################################
  ###'
  ###' Generate output list to save
  ###' 
  ###' 
  
  outpt = NULL
  outpt$nsimulations = nsimu
  outpt$grboundary = c(alow, ahigh)
  outpt$nmcmc = nmcmc
  outpt$nburn = nburn
  outpt$numintervalhis = nhis
  
  avgsel = avgssel  ## important: want this in the output
  
  tt = c(1, 4, 7, 10, 13, 16)
  
  tailprob = cbind(c(0.05, 0.10, 0.25, 0.75, 0.90, 0.95),
                   avgtail[tt], avgtail[c(tt + 1)], avgtail[c(tt + 2)])
  
  output = rbind(avgsel, avgisel, 
                 avgmean, avgvar, avgrankest)*10000 # to better read output
  
  fname <- simulparam$fname[LLind]
  fname1 = paste(fname, ".txt", sep = "")
  fname2 = paste(fname, "_hist.txt", sep = "")
  fname3 = paste(fname, "selrnk.txt", sep = "")
  fname4 = paste(fname, "trueout.txt", sep = "")
  
  avgselrank = 10000*avgselrank # just for output to be clearer
  
  
  ###'######################################################################
  ###'
  ###' Save as .txt files
  ###'
  ###'
  
  ### Set working directory
  setwd("~/Bayes-deconvolution/0_Base Camp")
  
  
  ### Outputs: SEL, ISEL, posterior means and variances, rank estimates
  
  dput(outpt,fname1)
  
  write.table(round(output, 2), 
              fname1, 
              append = TRUE, quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  
  
  write.table(round(tailprob, 2), 
              fname1, 
              append = TRUE, quote=FALSE,
              row.names = FALSE, col.names = FALSE)
  
  
  write.table(avghis, 
              fname2, 
              quote =FALSE,
              row.names = FALSE, col.names = FALSE)
  
  write.table(round(avgselrank, 2),
              fname3, 
              quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  
  
  ### Posterior variances
  tmp = c(LLind, 
          mean(evartheta), 
          c(range(evartheta)), 
          median(evartheta))
  
  tmp = matrix(tmp, nrow = 1)
  
  write.table(round(tmp,4), 
              fname4, 
              append=TRUE, quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  
  
  ###'######################################################################
  ###'
  ###' End of loop over LLind
  ###' 
  ###'
  
}    ### End of loop over LLind



###'######################################################################
###'
###' Create tables from saved results
###'
###'

corefc=function(nn)
{
  g1a=matrix(scan(nn,skip=8),ncol=4,byrow=T,quiet = TRUE)
  g1a=g1a[c(2:5),-1]
  g1b=matrix(scan(nn,skip=3),ncol=3,byrow=T,quiet = TRUE)[c(1,2,5),]
  g1=rbind(g1b,g1a)
  h1=g1
  h1[1:3,c(1:2)]=g1[1:3,c(1:2)]/g1[1:3,3]*100
  ind=c(3,1,2)
  return(h1[,ind])
}

i <- with(simulparam, order(n,rgt,rsr,truth,assumed))
simulparam <- simulparam[i,]

temp <- apply(simulparam[,c("rgt","rsr","n")],1,paste,collapse="+")
temp <- as.numeric(factor(temp,levels=unique(temp),ordered=TRUE))
simulparam$group <- temp

nind = 1
for(i in unique(simulparam$group))
{
  A <- with(subset(simulparam,group==i),
            split(fname,truth))
  A <- lapply(A,function(x){paste(x,".txt",sep="")})
  
  all <- NULL
  for(j in 1:length(A))
  {
    all1 <- NULL
    for(k in 1:length(A[[j]]))
    {
      temp <- corefc(A[[j]][k])
      if(k>1) temp <- temp[,-1]
      all1 <- cbind(all1,temp[-3,]) # remove selr
    }
    all <- rbind(all,all1)
  }
  
  print(round(all))
}








