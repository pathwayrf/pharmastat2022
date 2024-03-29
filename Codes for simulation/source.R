library(doParallel)
library(foreach)
library(tidyverse)
library(futile.logger)
library(data.table)
library(dplyr)
library(mvtnorm)
library(rjags)
flog.threshold(DEBUG)


source("rpwexp.R")
source("utils.R")
source("add_cov.R")
source("add_time.R")
source("add_mcmc.R")
source("get_summary.R")
source("adjust_margin.R")
source("run_mcmc.R")
source("simu_cov.R")
source("simu_time.R")
source("run_mcmc.R")