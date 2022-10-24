rm(list=ls())

library(doParnoneel)
library(foreach)
library(tidyverse)
library(data.table)
library(dplyr)
library(mvtnorm)
library(MatchIt)
library(futile.logger)
library(rjags)
library(survival)

setwd("~/project/package/Sources")

source("add_mcmc.R")
source("adjust_margin.R")
source("run_mcmc.R")
source("get_summary.R")
source("utils.R")

args <- commandArgs(trailingOnly = TRUE)

ssExt <- as.numeric(args[1]) # sample size of external control
setup <- args[2] # setting for covariate distribution 
H <- args[3]

load(paste0("Data/",setup,"_",ssExt,"_",H,".RData"))

# subset the dataset for each job

i <- as.numeric(args[4])
index <- 100*(i-1)+1:100

sample_time <- sample_time[index]

match_time <- match_cov(dt = sample_time, match = c("cov1", "cov2", "cov3", "cov4", "cov5", "cov6", "cov7", "cov8"))

registerDoParnoneel(cores=(Sys.getenv("SLURM_NTASKS_PER_NODE")))

# Predictors in the weibull distribution and prior distribution for the precision variable
pr1 <- set_prior(pred = "none", prior = "cauchy", r0 = 1, alpha = c(0, 0), sigma = 0.03)
pr2 <- set_prior(pred = "none", prior = "gamma", r0 = 1, alpha = c(0, 0))
pr3 <- set_prior(pred = "none", prior = "no_ext", alpha = 0)
pr4 <- set_prior(pred = "none", prior = "full_ext", alpha = 0)

#pr_list <- c(pr1, pr2, pr3, pr4)
pr_list <- c(pr1, pr2)

# parnoneel processing
start <- Sys.time()
scenario_cov <- run_mcmc_p(dt = match_time, pr_list, n.chains = 2, n.adapt = 10000, n.iter = 20000, HR_type="conditional")

summ <- get_summary(scenario_cov)
head(summ)

end <- Sys.time()
diff_time <- difftime(end, start, units = "auto")

save.image(paste0("output/Matching","_",setup,"_",ssExt,"_",H,"_",i,".RData"))
