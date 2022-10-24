rm(list=ls())

library(tidyverse)
library(data.table)
library(dplyr)
library(mvtnorm)
library(MatchIt)
library(futile.logger)
#library(rjags)

setwd("/Users/Stacey1/Documents/2022 Spring/Intern project manuscript/codes/DLBCL_excWkspace/source")

source("simu_cov.R")
source("add_cov.R")
source("add_time.R")
source("simu_time.R")
source("utils.R")
source("rpwexp.R")
source("match_fun.R")

setup.list <- c("L_diff","No_diff","L_rev")
ssExt.list <- c(100,200,400)
HR.list <- c(0.67,1)
for(HR in HR.list){
  for(ssExt in ssExt.list){
    ss = set_n(ssC = 140, ssE = 275, ssExt = ssExt) 
    for(setup in setup.list){
      if(setup == "L_rev"){
        #### Large diff reverse 
        covset1 = set_cov(n_cat = 2, n_cont = 0, mu_int = rep(0, 2), mu_ext = rep(0, 2), var = rep(1, 2), cov = 0.5, 
                          prob_int = c(0.75, 0.56), prob_ext = c(0.85, 0.64))
        covset2 = set_cov(n_cat = 2, n_cont = 0, mu_int = rep(0, 2), mu_ext = rep(0, 2), var = rep(1, 2), cov = 0.3, 
                          prob_int = c(0.31, 0.54), prob_ext = c(0.36, 0.5))
        covset3 = set_cov(n_cat = 3, n_cont = 0, mu_int = rep(0, 3), mu_ext = rep(0, 3), var = rep(1, 3), cov = 0.2, 
                          prob_int = c(0.84, 0.38, 0.18), prob_ext = c(0.86, 0.36, 0.15))
        covset4 = set_cov(n_cat = 0, n_cont = 1, mu_int = 63, mu_ext = 62.6, var = 11.9^2)
        
        # Enrollment pattern, drop-out, analysis start time
        c_int = set_clin(gamma = 10, e_itv = 415/10, etaC = -log(1-0.04/12),  CCOD = "fixed-first", CCOD_t = 60)
        c_ext = set_clin(gamma = ssExt/18, e_itv = 18, etaC = -log(1-0.01/12),  CCOD = "fixed-first", CCOD_t = 70)
        
      }else if(setup == "No_diff"){
        ## No diff
        covset1 = set_cov(n_cat = 2, n_cont = 0, mu_int = rep(0, 2), mu_ext = rep(0, 2), var = rep(1, 2), cov = 0.5, 
                          prob_int = c(0.75, 0.56), prob_ext = c(0.75, 0.56))
        covset2 = set_cov(n_cat = 2, n_cont = 0, mu_int = rep(0, 2), mu_ext = rep(0, 2), var = rep(1, 2), cov = 0.3, 
                          prob_int = c(0.31, 0.54), prob_ext = c(0.31, 0.54))
        covset3 = set_cov(n_cat = 3, n_cont = 0, mu_int = rep(0, 3), mu_ext = rep(0, 3), var = rep(1, 3), cov = 0.2, 
                          prob_int = c(0.84, 0.38, 0.18), prob_ext = c(0.84, 0.38, 0.18))
        covset4 = set_cov(n_cat = 0, n_cont = 1, mu_int = 63, mu_ext = 63, var = 11.9^2)
        
        # Enrollment pattern, drop-out, analysis start time
        c_int = set_clin(gamma = 10, e_itv = 415/10, etaC = -log(1-0.04/12),  CCOD = "fixed-first", CCOD_t = 60)
        c_ext = set_clin(gamma = ssExt/18, e_itv = 18, etaC = -log(1-0.01/12),  CCOD = "fixed-first", CCOD_t = 60)
        
      }else if(setup == "L_diff"){
        ## Large diff
        covset1 = set_cov(n_cat = 2, n_cont = 0, mu_int = rep(0, 2), mu_ext = rep(0, 2), var = rep(1, 2), cov = 0.5, 
                          prob_int = c(0.75, 0.56), prob_ext = c(0.65, 0.48))
        covset2 = set_cov(n_cat = 2, n_cont = 0, mu_int = rep(0, 2), mu_ext = rep(0, 2), var = rep(1, 2), cov = 0.3, 
                          prob_int = c(0.31, 0.54), prob_ext = c(0.26, 0.58))
        covset3 = set_cov(n_cat = 3, n_cont = 0, mu_int = rep(0, 3), mu_ext = rep(0, 3), var = rep(1, 3), cov = 0.2, 
                          prob_int = c(0.84, 0.38, 0.18), prob_ext = c(0.82, 0.40, 0.21))
        covset4 = set_cov(n_cat = 0, n_cont = 1, mu_int = 63, mu_ext = 63.4, var = 11.9^2)
        
        # Enrollment pattern, drop-out, analysis start time
        c_int = set_clin(gamma = 10, e_itv = 415/10, etaC = -log(1-0.04/12),  CCOD = "fixed-first", CCOD_t = 60)
        c_ext = set_clin(gamma = ssExt/18, e_itv = 18, etaC = -log(1-0.01/12),  CCOD = "fixed-first", CCOD_t = 52)
        
      }
      cov_list = c(covset1, covset2, covset3, covset4)
      
      sample_cov <- simu_cov(ssObj = ss, covObj = cov_list, HR = HR, driftHR = 1,nsim=2000)
      
      # Time-to-event distribution
      evt <- set_event(event = "weibull", shape = 0.9, lambdaC = 0.0135, beta = c(1, 0.5, 0.5, rep(0.001, 5)))
      sample_time <- simu_time(dt = sample_cov, eventObj = evt, clinInt = c_int, clinExt = c_ext)
      
      # PS matching
      #match_time <- match_cov(dt = sample_time, match = c("cov1", "cov2", "cov3", "cov4", "cov5", "cov6", "cov7", "cov8"))
      
      if(HR==1){
        save.image(paste0("Data/",setup,"_",ssExt,"_H0",".RData"))
      }else{
        save.image(paste0("Data/",setup,"_",ssExt,"_H1",".RData"))
      }
    }
  }
}

