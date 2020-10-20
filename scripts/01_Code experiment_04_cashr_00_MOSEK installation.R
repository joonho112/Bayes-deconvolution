
###'######################################################################
###'
###' Category: Code Experiment
###' 
###' Task: cashr
###'       
###'       (1) Install MOSEK, Rmosek, REBayes, and cashr
###'       
###' Data: NULL
###' 
###' Data: 2020-03-07
###' 
###' Author: JoonHo Lee (joonho@berkeley.edu)
###' 
###' 

###' Intall the following programs, outside R
###' 
###'  (1) Rtools (the tools needed for R package development)
###'  (2) MOSEK (the optimization library we interface to)
###'      https://www.mosek.com/downloads/
###'      


### Install "Rmosek" package 
install.packages("Rmosek", type="source")


### Install mosek from Rmosek builder
library("Rmosek")
mosek_attachbuilder("C:/Program Files/Mosek/9.1/tools/platform/win64x86/bin")
install.rmosek()


### Install "REBayes" and "cashr"
library(devtools)
install.packages("REBayes")
devtools::install_github("LSun/cashr")
