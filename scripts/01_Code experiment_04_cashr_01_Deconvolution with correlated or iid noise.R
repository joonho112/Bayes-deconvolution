
###'######################################################################
###'
###' Category: Code Experiment
###' 
###' Task: cashr
###'       
###'       (1) Deconvolution with correlated or iid noise
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
library(ashr)
library(cashr)
library(deconvolveR)
library(EbayesThresh)


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
###' Read in data
###'
###'

setwd(work_dir)


###' RNA-seq gene expression data on the 10^4 
###' most highly expressed genes in 119 human liver tissues
z.sel <- readRDS("datasets/cashr/z.sel.rds")
dim(z.sel)


###' Computed effective standard deviations
s <- readRDS('datasets/cashr/s.exp.rds')



###'######################################################################
###'
###' Define TRUE quantities  
###'  
###'  

###' Define true effect size prior g
###' This is equivalent to the G, in my setting
###' Normal means \theta are iid samples form the mixture
###' g(.) = 0.6 delta_0 + 0.3 N(.; 0, 1) + 0.1 N(.; mu, sigma^2)
G <- function (t) {
  0.6 * pnorm(t, 0, 0) + 0.3 * pnorm(t, 0, 1) + 0.1 * pnorm(t, 0, 3)
}


### Generate true normal means (total N = 10^4)
set.seed(777)

theta <- sample(c(
  rnorm(6e3, 0, 0),
  rnorm(3e3, 0, 1),
  rnorm(1e3, 0, 3)
))

hist(theta)
plot(density(theta))


### X grids
x.plot <- seq(-6, 6, by = 0.1)
plot(density(x.plot))
hist(x.plot)

G.plot <- G(x.plot)
plot(density(G.plot))
hist(G.plot)



###'######################################################################
###' 
###' Deconvolution
###' 
###' 

### Prepare loops
noise.label <- paste0("(", letters[1 : 5], ")")

deconv.list <- list()


### Start loops
for (i in 1:5) {
  
  
  ###'######################################################################
  ###' 
  ###' Generate a true dataset
  ###' 
  ###' 
  
  ### Get Z
  if (i <= 4) {
    Z <- z.sel[i, ]
  } else {
    set.seed(777)
    Z <- rnorm(1e4)
  }
  
  ###' Generate sample X and Z 
  ###' What's the difference between X and z
  X <- theta + s * Z   
  z <- theta + Z
  
  ### Truth
  true.data <- cbind.data.frame(
    method = "True g",  # label
    x = x.plot,         # X-grid
    cdfhat = G.plot     # CDF
  )
  
  
  ###'######################################################################
  ###' 
  ###' ashr: Methods for Adaptive shrinkage   
  ###' 
  ###' The main assumption is that the underlying distribution of effects 
  ###' is unimodal
  ###' 
  ###' 
  
  ### Fit the model for adaptive shrinkage
  fit.ashr <- ashr::ash(X,   # a p vector of estimates 
                        s,   # a p vector of corresponding standard errors 
                        
                        # how ash is run - false discovery rate 
                        method = "fdr", 
                        
                        # distribution of components in mixture
                        mixcompdist = "normal"  
                        )
  
  ### Get fitted G
  fitted_g <- get_fitted_g(fit.ashr)
  hist(fitted_g$pi)
  
  
  ### Get CDF
  ashr.plot <- as.numeric(ashr::mixcdf(ashr::get_fitted_g(fit.ashr), x.plot))
  plot(density(ashr.plot))
  
  ashr.data <- cbind.data.frame(
    method = "ashr",
    x = x.plot,
    cdfhat = ashr.plot
  )
  
  
  ###'######################################################################
  ###' 
  ###' cashr: solving Empirical Bayes Normal Means problem with correlated noise
  ###' 
  ###' => doesn't work. returns error. 
  ###' 
  ###' 
  
  ### Fit the model
  fit.cashr <- cashr::cash(X,  # a p vector of observations 
                           s)   # A p vector of standard deviations (SE)
                           
  
  ### Get CDF
  cashr.plot <- as.numeric(ashr::mixcdf(ashr::get_fitted_g(fit.cashr), x.plot))
  
  cashr.data <- cbind.data.frame(
    method = "cashr",
    x = x.plot,
    cdfhat = cashr.plot
  )
  
  
  
  ###'######################################################################
  ###' 
  ###' deconvolveR: Efron prior
  ###' 
  ###' Compute Empirical Bayes estimates using deconvolution
  ###' 
  ###' 

  ### Fit the model deconvolveR
  fit.deconvolveR <- deconvolveR::deconv(
    tau = x.plot,     # a vector of discrete support points for theta 
    X = z,            # a vector of sample values
    family = "Normal", 
    deltaAt = 0)      # the theta value where a delta function is desired
  
  
  ### Get CDF 
  deconvolveR.data <- cbind.data.frame(
    method = "deconvolveR",
    x = fit.deconvolveR$stats[, 1],
    cdfhat = fit.deconvolveR$stats[, 4]
  )
  
  
  ###'######################################################################
  ###' 
  ###' Kiefer-Wolfowitz's NPMLE (1956)
  ###' implemented by Koenker-Mizera-Gu's REBayes (2016)
  ###' 
  ###' => Doesn't work
  ###' 
  ###' 
  
  ### Define x grid
  v = seq(-6.025, 6.025, by = 0.05)
  
  
  ### Fit the model
  fit.REBayes <- REBayes::GLmix(x = X,     # Data: sample observations 
                                v = v,     # Undata: Grid values 
                                sigma = s) # vector of the Gaussian noise 
  
  
  ### Get CDF
  CDF.KW <- function(h, interp = FALSE, eps = 0.001, bw = 0.7){
    # Wasserstein distance:  ||G-H||_W
    if(interp == "biweight"){
      yk = h$x
      for (j in 1:length(yk))
        yk[j] = sum(biweight(h$x[j], h$x, bw = bw)*h$y/sum(h$y))
      H <- cumsum(yk)
      H <- H/H[length(H)]
    }
    else {
      H <- cumsum(h$y)
      H <- H/H[length(H)]
    }
    return(H)
  }
  
  REBayes.plot <- CDF.KW(fit.REBayes)
  
  REBayes.data <- cbind.data.frame(
    method = "REBayes",
    x = fit.REBayes$x,
    cdfhat = REBayes.plot
  )
  
  
  ###'######################################################################
  ###' 
  ###' EbayesThresh: Empirical Bayes Thresholding  
  ###'    
  ###'    
  
  ### Fit the model
  fit.EbayesThresh <- EbayesThresh::ebayesthresh(X,  # data vector 
                                                 sdev = s,  # standard deviation
                                                 verbose = TRUE, 
                                                 prior = "laplace", 
                                                 a = NA)
  
  ### Get CDF
  EbayesThresh.plot <- (1 - fit.EbayesThresh$w) * (x.plot >= 0) + 
    fit.EbayesThresh$w * rmutil::plaplace(x.plot, m = 0, s = 1 / fit.EbayesThresh$a)
  
  EbayesThresh.data <- cbind.data.frame(
    method = "EbayesThresh",
    x = x.plot,
    cdfhat = EbayesThresh.plot
  )
  
  
  
  ###'######################################################################
  ###' 
  ###' Collect CDFs
  ###' 
  ###' 
  
  deconv.list[[i]] <- cbind.data.frame(
    noise = noise.label[i],
    rbind.data.frame(
      true.data,
      EbayesThresh.data,
      REBayes.data,
      ashr.data,
      deconvolveR.data,
      cashr.data
    )
  )
}


###'######################################################################
###' 
###' Plotting
###' 
###' 

deconv.ggdata <- do.call(rbind.data.frame, deconv.list)
deconv.ggdata$noise <- factor(deconv.ggdata$noise,
                              levels = levels(deconv.ggdata$noise)[c(1, 2, 5, 3, 4)]
)

method.col <- c("black", scales::hue_pal()(5)[c(1, 3, 4, 2, 5)])
method.linetype <- rep(1, 6)


## plotting
deconv.plot <- ggplot(data = deconv.ggdata, aes(x = x, y = cdfhat, col = method, linetype = method)) +
  geom_line(size = 1) +
  facet_wrap(~noise, nrow = 2) +
  xlim(-5, 5) +
  scale_linetype_manual(values = method.linetype) +
  scale_color_manual(values = method.col) +
  labs(y = expression(paste("CDF of (estimated) g"))) +
  theme(plot.title = element_text(size = 15, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 15),
        legend.position = c(0.85, 0.25),
        legend.title = element_blank(),
        legend.text = element_text(size = 15))
ggsave("../figures/deconv.pdf", height = 6, width = 10)
