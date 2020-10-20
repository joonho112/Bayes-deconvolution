
###'######################################################################
###'
###' Category: Code Experiment
###' 
###' Task: HETOP
###'       
###'       (2) Look into Efron prior components within the HETOP function
###'       
###' Data: Simulated data
###' 
###' Data: 2020-03-07
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
gc(); rm(list=ls())


### Set working directory 
work_dir <- c("~/Bayes-deconvolution")
setwd(work_dir)


### Set a data directory
data_dir <- file.path(work_dir, "datasets")


### Call libraries
library(tidyverse)
library(cowplot)
library(HETOP)


### Call functions
list.files("functions", full.names = TRUE) %>% walk(source)



###'######################################################################
###' 
###' Generate data
###'
###'

set.seed(1001)

### Define mean-centered covariates
G  <- 12
z1 <- sample(c(0, 1), size = G, replace = TRUE)
z2 <- 0.5*z1 + rnorm(G)
Z  <- cbind(z1 - mean(z1), z2 = z2 - mean(z2))


### Define true parameters dependent on covariates
beta_m    <- c(0.3,  0.8)
beta_s    <- c(0.1, -0.1)
mug       <- Z[, 1]*beta_m[1] + Z[, 2]*beta_m[2] + rnorm(G, sd = 0.3)
sigmag    <- exp(0.3 + Z[,1]*beta_s[1] + Z[,2]*beta_s[2] + 0.2*rt(G, df = 7))
cutpoints <- c(-1.0, 0.0, 1.2)


### Generate data
ng   <- rep(200, G)
ngk  <- gendata_hetop(G, K = 4, ng, mug, sigmag, cutpoints)
print(ngk)



###'######################################################################
###' 
###'  Assign values for arguments
###' 
###' 

### Data
ngk


### The first two cutpoints
fixedcuts = c(-1.0, 0.0)


### Arguments needed to prameterize Efron prior
p = c(10,10)     # degrees of freedom for cubic spline basis  
m = c(100, 100)  # number of grid points  
gridL = c(-5.0, log(0.10))  # lower bounds for grids
gridU = c(5.0, log(5.0))    # upper bounds for grids


### Covariates
Xm = Z  
Xs = Z


### MCMC & JAGS options
n.iter = 100 
n.burnin = 50
seed = 12345
modelfileonly = FALSE
modloc = NULL



###'######################################################################
###' 
###' (1) Prepare values for data
###' 
###' 

set.seed(seed)


### Set temporary directory and model file name
tmpdir <- tempdir()

if (is.null(modloc)) {
  modloc <- paste0(tmpdir, "/model.txt")
}


### Extract the number of groups (G) and the number of categories (K)
G <- nrow(ngk)
K <- ncol(ngk)


### Get rowwise sum & group probabilities
ng <- apply(ngk, 1, sum)
pg <- ng/sum(ng)


### Get cutpoints
cuts12 <- sort(fixedcuts)

if (K == 3) {
  
  cuts <- cuts12

} else {
  
  cuts <- c(cuts12, rep(NA, K - 3))

}


### Extract dimensions for covariates
meanX <- !is.null(Xm)
sdX <- !is.null(Xs)

dimXm <- ncol(Xm)
dimXs <- ncol(Xs)



###'######################################################################
###' 
###' (2) Prepare arguments for JAGS model 
###' 
###' 

### Define valid JAGS arguments
valid.ja <- c("inits", "parameters.to.save", 
              "n.chains", "n.iter", "n.burnin", "n.thin", 
              "DIC", "working.directory", "refresh", 
              "progress.bar", "digits", "RNGname", 
              "jags.module")


### JAGS parameters to save
jags.parameters.to.save <- c("mu", "sigma", 
                             "cuts", "alpha0m", "alpham", "alpha0s", 
                             "alphas", "gamma")

if (meanX && sdX) {
  jags.parameters.to.save <- c(jags.parameters.to.save, 
                               "beta_m", "beta_s")
}


### Number of chains, iterations, burn-ins, thin, etc.
jags.n.chains <- 2
jags.n.iter <- 5000
jags.n.burnin <- jags.n.iter/2
jags.n.thin <- 1
jags.DIC <- FALSE
jags.working.directory <- NULL
jags.refresh <- jags.n.iter/20
jags.progress.bar <- "text"
jags.digits <- 5
jags.RNGname <- c("Wichmann-Hill", "Marsaglia-Multicarry", 
                  "Super-Duper", "Mersenne-Twister")
jags.jags.module = c("glm", "dic")



###'######################################################################
###' 
###'  (3) Prepare values for JACS data
###' 
###' 

### Get Q matrix for mean
gridm <- seq(from = gridL[1], to = gridU[1], length = m[1])

Qm <- bs(gridm,      # the predictor variable 
         df = p[1],  # degrees of freedom 
         degree = 3, # degree of the piecewise polynomial - 3: cubic spline
         intercept = FALSE)

Qm <- matrix(c(Qm), ncol = p[1], byrow = F)


### Get Q matrix for SD
grids <- seq(from = gridL[2], to = gridU[2], length = m[2])

Qs <- bs(grids, df = p[2], degree = 3, intercept = FALSE)

Qs <- matrix(c(Qs), ncol = p[2], byrow = F)


### Get a list of JAGS data elements
jags.data <- c("G", "K", "ngk", "ng", 
               "cuts", "m", "p", "Qm", "Qs", 
               "gridm", "grids")

if (meanX && sdX) {
  jags.data <- c(jags.data, "Xm", "dimXm", 
                 "Xs", "dimXs")
}


### Get initial values
jags.inits <- vector(jags.n.chains, mode = "list")

for (i in 1:jags.n.chains) {
  
  .tmp <- list(locm = rep(as.integer(floor(0.5 * m[1])), G), 
               locs = as.integer(rep(floor(0.95 * m[2]), G)), 
               alpha0m = rnorm(1), alphamdev = rnorm(p[1]), 
               alpha0s = rnorm(1), alphasdev = rnorm(p[2]), 
               gamma = 0)
  
  if (meanX) {
    .tmp$beta_m <- rnorm(dimXm, sd = 0.1)
  }
  
  if (sdX) {
    .tmp$beta_s <- rnorm(dimXs, sd = 0.1)
  }
  
  if (K >= 4) {
    .tmp$cuts0 <- sort(seq(from = cuts12[2] + 1, to = cuts12[2] + 3, length = K - 3))
  }
  jags.inits[[i]] <- .tmp
}



###'######################################################################
###' 
###'  (4) Construct a JAGS model
###' 
###' 

cat("model\n  {\n    for(g in 1:G){\n      ngk[g,1:K] ~ dmulti(pgk[g,1:K], ng[g])\n\n      pgk[g,1] <- phi( (cuts[1] - mu[g]) / sigma[g])\n      for(k in 2:(K-1)){\n        pgk[g,k] <- phi( (cuts[k] - mu[g]) / sigma[g]) - sum(pgk[g,1:(k-1)])\n      }\n      pgk[g,K] <- 1.0 - phi( (cuts[K-1] - mu[g]) / sigma[g])\n\n      mu[g]    <-     epsilon[g,1]\n      sigma[g] <- exp(epsilon[g,2])\n\n      locm[g] ~ dcat(pFm[1:m[1]])\n      locs[g] ~ dcat(pFs[1:m[2]])\n", 
    file = modloc)
if ((!meanX && !sdX)) {
  cat("\n      epsilon[g,1]  <- (gamma * grids[locs[g]]) + gridm[locm[g]]\n      epsilon[g,2]  <- grids[locs[g]]\n", 
      file = modloc, append = TRUE)
}
if ((meanX && sdX)) {
  cat("\n      epsilon[g,1] <- inprod(Xm[g,1:dimXm], beta_m[1:dimXm]) + (gamma * grids[locs[g]]) + gridm[locm[g]]\n      epsilon[g,2] <- inprod(Xs[g,1:dimXs], beta_s[1:dimXs]) + grids[locs[g]]\n", 
      file = modloc, append = TRUE)
}
cat("    }\n", file = modloc, append = TRUE)
cat("\n    alpha0m     ~ dnorm(0.0, 0.01)\n    for(i in 1:p[1]){\n       alphamdev[i] ~ dnorm(0.0, 0.05)\n       alpham[i] <- alpha0m + alphamdev[i]\n    }\n    for(i in 1:m[1]){\n      pFm[i] <- exp(Qm[i,1:p[1]] %*% alpham[1:p[1]])\n    }\n\n    alpha0s     ~ dnorm(0.0, 0.01)\n    for(i in 1:p[2]){\n       alphasdev[i] ~ dnorm(0.0, 0.05)\n       alphas[i] <- alpha0s + alphasdev[i]\n    }\n    for(i in 1:m[2]){\n      pFs[i] <- exp(Qs[i,1:p[2]] %*% alphas[1:p[2]])\n    }\n\n    gamma ~ dnorm(0.0, 0.1)\n\n    ", 
    file = modloc, append = TRUE)
if (K >= 4) {
  cat("\n    for(k in 1:(K-3)){\n      cuts0[k] ~ dnorm(0.0, 0.01) I(cuts[2], )\n    }\n    tmp[1:(K-3)] <- sort(cuts0)\n    for(k in 1:(K-3)){\n      cuts[k+2] <- tmp[k]\n    }\n    ", 
      file = modloc, append = TRUE)
}
if (meanX) {
  cat("\n    for(i in 1:dimXm){\n      beta_m[i] ~ dnorm(0.0, 0.1)\n    }\n    ", 
      file = modloc, append = TRUE)
}
if (sdX) {
  cat("\n    for(i in 1:dimXs){\n      beta_s[i] ~ dnorm(0.0, 0.1)\n    }\n    ", 
      file = modloc, append = TRUE)
}
cat("}\n", file = modloc, append = TRUE)
if (modelfileonly) {
  return(modloc)
}



###'######################################################################
###' 
###'  Fit the constructed JAGS model & save the resulting objects
###' 
###' 

r <- jags(model.file = modloc, data = jags.data, inits = jags.inits, 
          parameters.to.save = jags.parameters.to.save, n.chains = jags.n.chains, 
          n.iter = jags.n.iter, n.burnin = jags.n.burnin, n.thin = jags.n.thin, 
          DIC = jags.DIC, working.directory = jags.working.directory, 
          refresh = jags.refresh, progress.bar = jags.progress.bar, 
          digits = jags.digits, RNGname = jags.RNGname, jags.module = jags.jags.module)

fh_hetop_extras <- list()
Finfo <- list()
Finfo$gridL <- gridL
Finfo$gridU <- gridU
Finfo$efron_p <- p
Finfo$efron_m <- m
Finfo$gridm <- gridm
Finfo$Qm <- Qm
Finfo$grids <- grids
Finfo$Qs <- Qs
fh_hetop_extras$Finfo <- Finfo
rm(Finfo)
Dinfo <- list()
Dinfo$G <- G
Dinfo$K <- K
Dinfo$ngk <- ngk
Dinfo$ng <- ng
Dinfo$fixedcuts <- fixedcuts
Dinfo$Xm <- Xm
Dinfo$Xs <- Xs
fh_hetop_extras$Dinfo <- Dinfo
rm(Dinfo)
fh_hetop_extras$waicinfo <- waic_hetop(ngk, r$BUGSoutput$sims.matrix)
ind_m <- grep("mu", colnames(r$BUGSoutput$sims.matrix))
ind_s <- grep("sigma", colnames(r$BUGSoutput$sims.matrix))
ind_c <- grep("cuts", colnames(r$BUGSoutput$sims.matrix))
ind_betam <- grep("beta_m", colnames(r$BUGSoutput$sims.matrix))
ind_betas <- grep("beta_s", colnames(r$BUGSoutput$sims.matrix))
stopifnot((length(ind_m) == G) && (length(ind_s) == G) && 
            length(ind_c == (K - 1)))
tmp <- lapply(1:nrow(r$BUGSoutput$sims.matrix), function(i) {
  mug <- r$BUGSoutput$sims.matrix[i, ind_m]
  sigmag <- r$BUGSoutput$sims.matrix[i, ind_s]
  cuts <- r$BUGSoutput$sims.matrix[i, ind_c]
  a <- sum(pg * mug)
  b <- sqrt(sum(pg * ((mug - a)^2 + sigmag^2)))
  .retval <- list(mug = (mug - a)/b, sigmag = sigmag/b, 
                  cutpoints = (cuts - a)/b)
  if (meanX && sdX) {
    .retval$beta_m <- r$BUGSoutput$sims.matrix[i, ind_betam]/b
    .retval$beta_s <- r$BUGSoutput$sims.matrix[i, ind_betas]
  }
  return(.retval)
})
fh_hetop_extras$est_star_samps <- list(mug = do.call("rbind", 
                                                     lapply(tmp, function(x) {
                                                       x$mug
                                                     })), sigmag = do.call("rbind", lapply(tmp, function(x) {
                                                       x$sigmag
                                                     })), cutpoints = do.call("rbind", lapply(tmp, function(x) {
                                                       x$cutpoints
                                                     })))
if (meanX && sdX) {
  fh_hetop_extras$est_star_samps$beta_m <- do.call("rbind", 
                                                   lapply(tmp, function(x) {
                                                     x$beta_m
                                                   }))
  fh_hetop_extras$est_star_samps$beta_s <- do.call("rbind", 
                                                   lapply(tmp, function(x) {
                                                     x$beta_s
                                                   }))
}
rm(tmp)
gc()
stopifnot(max(abs(sapply(1:nrow(r$BUGSoutput$sims.matrix), 
                         function(i) {
                           sum(pg * (fh_hetop_extras$est_star_samps$mug[i, ]^2 + 
                                       fh_hetop_extras$est_star_samps$sigmag[i, ]^2))
                         }) - 1)) < 1e-12)
fh_hetop_extras$est_star_mug <- triple_goal(fh_hetop_extras$est_star_samps$mug)
fh_hetop_extras$est_star_sigmag <- triple_goal(fh_hetop_extras$est_star_samps$sigmag)
r$fh_hetop_extras <- fh_hetop_extras
rm(fh_hetop_extras)
gc()
return(r)
}



