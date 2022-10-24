# S4 class, constructors and methods for running simulations for multiple scenarios
# Last update on 102021

#' @import doParallel
#' @import foreach
#' @import rjags

`%notin%` <- Negate(`%in%`)

run_mcmc_np <- function(dt, priorObj, n.chains, n.adapt, n.iter, seed, path, HR_type="marginal"){
  
  if (missing(dt)) stop("Please provide a list of simulated time (dt).")
  if (missing(priorObj)) stop("Please provide .priorObj (priorObj).")
  n_mcmc <- valid_mcmc(n.chains, n.adapt, n.iter)
  
  if (missing(seed)){
    message("Set.seed(47)")
    seed = 47
  }
  seed_list <- seq(seed, seed + length(dt), by = 1)
  
  #flog.debug(cat("seed_list:", seed_list, "number of dataset", length(dt), "\n"))
  
  # nCluster <- parallel::detectCores()
  # print(paste(nCluster, "clusters are being used"))
  # cl <- parallel::makeCluster(nCluster)
  # doParallel::registerDoParallel(cl)
  
  
  # res_list <- foreach(i = 1:length(dt), .combine='c', .multicombine=TRUE,
  #                     # .export = c("ad_est","adjust_margin"),
  #                     # .export = c("add_mcmc","ad_est","r_post","rej_est","ad_est","adjust_margin"),
  #                     # .export = c("format_number", "format_date_time", "add_direction_data"),
  #                     # .packages = c("tidyverse", "data.table", "dplyr", "rjags"),
  #                     .verbose=FALSE) %dopar% {
  #                       
  #                       source("source.R")
  #                       
  #                       
  #                     } #loop foreach
  # 
  # parallel::stopCluster(cl)
  
  for(i in 1:length(dt)){
    seed_i = 1
    
    print(paste("-------------------", 1, "of ", length(dt), "simulated dataset with seed =", seed_i))
    
    
    if(i==1){
      res_list <- add_mcmc(dt = dt[[1]], priorObj = priorObj,
                           n.chains = n.chains, n.adapt = n.adapt,
                           n.iter = n.iter, seed = seed_i, HR_type = HR_type)
    }else{
      temp <- add_mcmc(dt = dt[[1]], priorObj = priorObj,
                       n.chains = n.chains, n.adapt = n.adapt,
                       n.iter = n.iter, seed = seed_i, HR_type = HR_type)
      res_list <- append(res_list,temp)
    }
  }
  
  
  sum_list <- lapply(res_list, function(i) {
    if(HR_type=="conditional"){
      cbind("HR" = i[['HR']], "driftHR" = i[['driftHR']],
            "prior" = i[['prior']], "pred" = i[['pred']],
            "reject" = i[['summary']]['reject'],
            "mean_HR" = i[['summary']]['mean_HR'],
            "sd_HR" = i[['summary']]['sd_HR'],
            "mean_driftHR" = i[['summary']]['mean_driftHR'],
            "sd_driftHR" = i[['summary']]['sd_driftHR'])
    }else{
      cbind("HR" = i[['HR']], "driftHR" = i[['driftHR']],
            "prior" = i[['prior']], "pred" = i[['pred']],
            "reject" = i[['summary_admargin']]['reject'],
            "mean_HR" = i[['summary_admargin']]['mean_HR'],
            "sd_HR" = i[['summary_admargin']]['sd_HR'],
            "mean_driftHR" = i[['summary']]['mean_driftHR'],
            "sd_driftHR" = i[['summary']]['sd_driftHR'],
            "p_value" = i[['summary_admargin']]['p_value']) 
    }
  })
  
  sum_dt <- as.data.frame(do.call(rbind, sum_list))
  rownames(sum_dt) <- NULL
  sum_dt$HR = as.numeric(sum_dt$HR)
  sum_dt$driftHR = as.numeric(sum_dt$driftHR)
  sum_dt$reject = as.numeric(sum_dt$reject)
  sum_dt$mean_HR = as.numeric(sum_dt$mean_HR)
  sum_dt$sd_HR = as.numeric(sum_dt$sd_HR)
  sum_dt$mean_driftHR = as.numeric(sum_dt$mean_driftHR)
  sum_dt$sd_driftHR = as.numeric(sum_dt$sd_driftHR)
  
  if (missing(path)) print("Samples from the posterior distribution from MCMC are not saved.") else {
    save(sum_dt, file = path)
    print(paste("Results the posterior distribution from MCMC are saved as", path))
  }
  sum_dt
}
