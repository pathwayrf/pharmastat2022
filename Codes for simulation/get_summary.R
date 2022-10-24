# Methods for calculating summary statistics for MCMC
# Last update on 082320

#' @import foreach

weird = 1

#' Generate summary statistics for the MCMC chains
#'
#' @param samples an object of class \code{mcmc.list}
#' @return a vector containing the mean, median, sd, reject rate for the MCMC chains
#'
#' @keywords method
rej_est = function(samples){
  summ = summary(samples)[[1]]
  quantiles <- summary(samples)[[2]]
  reject = quantiles["HR_trt_cc","97.5%"] < 1
  mean_HR  =  summ["HR_trt_cc","Mean"]
  sd_HR = summ["HR_trt_cc","SD"]
  # post_median_hr_trt_cc = quantiles["HR_trt_cc","50%"]
  mean_driftHR  =  summ["HR_cc_hc", "Mean"]
  sd_driftHR = summ["HR_cc_hc","SD"]
  print(paste("reject", reject,
              "mean_HR", mean_HR, "sd_HR", sd_HR,
              "mean_driftHR", mean_driftHR,  "sd_driftHR", sd_driftHR))
  c("reject" = reject,
    "mean_HR" = mean_HR, "sd_HR" = sd_HR,
    "mean_driftHR"  = mean_driftHR,  "sd_driftHR"  = sd_driftHR)
}


#' Generate summary statistics of a simulation scenario
#'
#' @param summ_file a \code{data.frame} containing summary statistics for the posterior samples from each simulation
#' @return a \code{data.frame} containing the mean, median, sd, reject rate, variance, sd follows standard definition, bias and mse of the simulation set
#'
#' @export
#' @keywords method
get_summary = function(dt){
  dt %>% group_by(HR, driftHR, prior, pred) %>%
    summarize(nsim = n(),
              reject = mean(reject),
              mean_HR_trt_cc = mean(mean_HR),
              sd_HR_trt_cc = mean(sd_HR),
              bias =  mean_HR_trt_cc - mean(HR),
              var = (1/(nsim - 1))*sum((mean_HR - mean_HR_trt_cc)^2),
              sd = sqrt(var),
              mse = var + bias^2,
              mean_HR_cc_hc = mean(mean_driftHR),
              sd_HR_cc_hc = mean(sd_driftHR))
}
