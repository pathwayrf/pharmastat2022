# S4 class, constructors and methods for generating posterior samples from MCMC
# Call function from utils.R
# Last update on 082320

#' @import futile.logger
#' @import survival

setClassUnion("charORNULL", c("character", "NULL"))
setClassUnion("numericORNULL", c("numeric", "NULL"))
setClassUnion("matrixORNULL", c("matrix", "NULL"))
`%notin%` <- Negate(`%in%`)

#' S4 Class for specifying prior distributions for MCMC methods
#'
#' @keywords class
.priorClass = setClass(".priorClass",
                       slots=list(pred = "character", prior = "charORNULL",
                                  r0 = "numericORNULL", alpha = "numericORNULL",  sigma = "numericORNULL" ))

#' Specify prior distributions for MCMC methods
#'

#' @param dt A \code{matrix} containing survival times, censoring, external trial indicator, treatment arm indicator, and covariates. Variables should have colnames "ext", "trt" and "cov+number" (e.g. cov1, cov2...)

#' @param prior Prior distribution for the precision parameter that controls the degree of borrowing. Half-cauchy distribution if \code{prior = "cauchy"}. No external data is included in the data if \code{prior = "no_ext"}. External control arm is assumed to have the same baseline hazards as internal control arm if \code{prior = "full_ext"}. Other options include "gamma" and "unif"
#' @param r0 Initial values for the shape of the weibull distribution for time-to-events
#' @param alpha Initial values for log of baseline hazard rate for external and internal control arms. Length of \code{alpha} should be 1 if \code{prior = "full_ext"} or \code{prior = "no_ext"}, and equal to 2 otherwise
#' @param sigma Initial values for precision parameter if \code{prior = "cauchy"}. If left \code{NULL}, default value 0.03 is used
#' @return a \code{.priorClass} class containing survival data and prior information
#'
#' @export
#' @keywords constructor
set_prior <- function(pred, prior, r0, alpha, sigma) {
  if (missing(pred)) {
    pred = "none"    
  }  

  if(missing(prior) || prior %notin% c("gamma", "cauchy", "no_ext", "full_ext", "unif")) {
    stop("Prior distribution for the precision parameter (prior) is not correctly specified. Options include gamma, cauchy, unif, no_ext, full_ext.")
  } else if (prior == "no_ext" & sum(grepl("ad", pred)) > 0) {
    stop("User choose to not include external information. Adjusting for propensity score is not an option.")
  }

  if(missing(r0)) {
    r0 = 1
    ## message('No initial values for the shape of the weibull distribution (r0) is detected. Default value 1 is used')
  }

  if(missing(alpha) || (prior %in% c("gamma", "cauchy", "unif") & length(alpha) != 2)) {
    alpha = c(0, 0)
    message('Values for log of baseline hazard rate for external and internal control arms (alpha) is not correctly specified. Default value 0 is used.')
  }
  if (missing(alpha) || (prior %in% c("no_ext", "full_ext") & length(alpha) != 1)){
    alpha = 0
    message('Values for log of baseline hazard rate for external and internal control arms (alpha) is not correctly specified. Default value 0 is used.')
  }
  if (missing(sigma)) sigma = NULL
  if (prior %in% c("cauchy", "unif") & length(sigma) != 1) {
    sigma = 0.03
    message("Initial value for precision parameter (sigma) is missing or not correctly specified. Default value 0.03 is used.")
  }
  new(".priorClass", pred = pred, prior = prior, r0 = r0, alpha = alpha, sigma = sigma)
}



#' Concatenate multiple \code{.priorClasss} class
#' @param dt A \code{matrix} containing survival times, censoring, external trial indicator, treatment arm indicator, and covariates. Variables should have colnames "ext", "trt" and "cov+number" (e.g. cov1, cov2...)
#' @param covObj a \code{.priorClasss} class with prior distributions for MCMC methods specified in \code{\link{set_prior}}
#' @return a vector of \code{.priorClasss} class
#'
#' @export
#' @keywords method
setMethod("c", signature(x = ".priorClass"), function(x, ...){
  elements = list(x, ...)
  priorClassList = list()
  for (i in 1:length(elements)){
    priorClassList[[i]] = new(".priorClass",
                              pred = slot(elements[[i]], "pred"), prior = slot(elements[[i]], "prior"),
                              r0 = slot(elements[[i]], "r0"), alpha = slot(elements[[i]], "alpha"),
                              sigma = slot(elements[[i]], "sigma"))
  }
  class(priorClassList) = ".priorClass"
  priorClassList
})


#' Generating posterior samples from MCMC
#'
#' @param priorObj a \code{.priorClass} class containing survival data and prior information
#' @param n.chains number of parallel chains for the model
#' @param n.adapt number of iterations for adaptation
#' @param n.iter number of iterations to monitor
#' @param seed the seed of R‘s random number generator. Default is 47
#' @return a single \code{mcmc.list} object.
#'
#' @export
#' @keywords method
setGeneric(name = "add_mcmc", def = function(dt, priorObj, n.chains, n.adapt, n.iter, seed, HR_type){standardGeneric("add_mcmc")})
setMethod(f = "add_mcmc", signature(dt = "matrix", priorObj = ".priorClass"),
          definition = function(dt, priorObj, n.chains, n.adapt, n.iter, seed, HR_type){

            if (missing(dt)) {
              stop ("Please provide a dataset (dt).")
            } else {
              cov_name = colnames(dt)[colnames(dt) %notin% c("driftHR", "HR", "trt", "ext", "extumatch", "time", "cnsr", "wgt")] #extract covariate names in the dataset
              if (sum(c("driftHR", "HR", "trt", "ext", "extumatch", "time", "cnsr") %notin% colnames(dt)) > 0 ) stop("Please make sure the trial data contains at least trt, ext, time and cnsr.")
              else if (sum(!grepl("cov[1234567890]", cov_name)) > 0) stop("Please make sure the covariates in the trial have the right name.")
            }
            flog.debug(cat("[add_mcmc] cov_name =", cov_name, "\n"))

            if (missing(seed)){
              message("Set.seed(47)")
              seed = 47
            }

            n_mcmc <- valid_mcmc(n.chains, n.adapt, n.iter)

            if (length(priorObj) == 1) priorObj = c(priorObj)

            lapply(seq(1, length(priorObj), by = 1), function(i){
              prior = priorObj[[i]]@prior
              pred = priorObj[[i]]@pred
              

              flog.debug(cat(">>> prior =", prior, ", pred = ", pred, "\n"))

              mcmc_res <- r_post(dt = dt,
                                 prior = prior, pred = pred,
                                 r0 = priorObj[[i]]@r0, alpha = priorObj[[i]]@alpha,
                                 sigma = priorObj[[i]]@sigma, seed = seed,

                                 n.chains = n_mcmc[['n.chains']], n.adapt = n_mcmc[['n.adapt']],
                                 n.iter = n_mcmc[['n.iter']])

              mcmc_sum <- rej_est(mcmc_res)
              flog.debug(cat(">>> mcmc_sum =", mcmc_sum, "\n"))
              cov_name = colnames(dt)[colnames(dt) %notin% c("driftHR", "HR", "trt", "ext", "time", "extumatch", "cnsr", "wgt")]
              
              if(HR_type=="marginal"){
                if (sum(grepl("all", pred)) > 0) {
                  cov = cov_name
                  #flog.debug(cat("All covariates", cov_name, "are used as predictors in the weibull distribution. Any other input is disregarded. \n"))
                }else{
                  cov = pred
                  #flog.debug(cat("Selected covariate(s)", pred, "are used as predictors in the weibull distribution. Any other input is disregarded. \n"))
                } 
                if( sum(grepl("none",pred))>0){
                  ad_sum <- NULL
                }else{
                  ad <- adjust_margin(mcmc_list = mcmc_res,dt=dt,pred=cov)
                  ad_sum <- ad_est(ad)
                }
              }else{
                ad_sum <- NULL
              }
              list("HR" = unique(dt[, 'HR']),
                   "driftHR" = unique(dt[, 'driftHR']),
                   "prior" = prior,
                   "pred" = paste(pred, collapse = " "), # print to one row
                   "mcmc.list" = mcmc_res, "summary" = mcmc_sum,"summary_admargin" = ad_sum)
            })
          })

valid_mcmc <- function(n.chains, n.adapt, n.iter){
  if (missing(n.chains) || !is.numeric(n.chains)) {
    n.chains = 2
    message("No number of parallel chains for the model is provided. Default value 2 is used")
  }
  if (missing(n.adapt) || !is.numeric(n.chains)) {
    n.adapt = 1000
    message("No Number of iterations for adaptation (n.adapt) is provided. Default value 1000 is used")
  }
  if (missing(n.iter) || !is.numeric(n.iter)) {
    n.iter = 10000
    message("No number of iterations to monitor is provided. Default value 10000 is used")
  }
  return(list("n.chains" = n.chains, "n.adapt" = n.adapt, "n.iter" = n.iter))
}
