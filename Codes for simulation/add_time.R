# S4 class, constructors and methods for simulating survival times
# Call function from utils.R
# Last update on 082320

#' @import futile.logger

setClassUnion("numericORNULL", c("numeric", "NULL"))
setClassUnion("charORNULL", c("character", "NULL"))
setClassUnion("listORNULL", c("list", "NULL"))
`%notin%` <- Negate(`%in%`)

#' S4 Class for specifying parameters for enrollment time, drop-out pattern and analysis start time
#'
#' @keywords class
.clinClass = setClass(".clinClass", slots=list(gamma ="numeric", e_itv = "numericORNULL", # enroll
                                               CCOD ="character", CCOD_t = "numeric",
                                               etaC = "numeric", etaE = "numeric", d_itv = "numericORNULL"))


#' Specify parameters for enrollment time, drop-out pattern and analysis start time
#'
#' This function allows user to specify the enrollment and drop-out rate, and the type of clinical cut-off Date. Both enrollment times and drop-out times follow piece-wise exponential distribution.
#'
#' @param gamma A vector of rate of enrollment per unit of time
#' @param e_itv A vector of duration of time periods for recruitment with rates specified in \code{gamma}. Note that the length of \code{e_itv} should be same length as \code{gamma} or 1 less.
#' @param CCOD Type of analysis start time. Analysis starts at \code{CCOD_t} months after the first or last patient's enrollment if \code{CCOD = "fixed-first"} or \code{CCOD = "fixed-last"} respectively. Analysis starts when \code{CCOD_t} events have been observed if \code{CCOD = "event"}
#' @param CCOD_t Time difference between analysis start and first patient's enrollment if \code{CCOD = "fixed-first"}. Time difference between analysis start and last patient's enrollment if \code{CCOD = "fixed-last"}. Number of events observed when analysis starts if \code{CCOD = "event"}. Patients enrolled after the analysis start time are excluded from the analysis
#' @param etaC A vector for dropout rate per unit time for control arm
#' @param etaE A vector for dropout rate per unit time for experimental arm. If left \code{NULL}, it uses the same dropout rate as eta.
#' @param d_itv A vector of duration of time periods for dropping out the study with rates specified in \code{etaC} and \code{etaE}. Note that the length of \code{d_itv} should be same length as \code{etaC} or 1 less.
#' @return A \code{.clinClass} class containing information on enrollment time, drop-out pattern and analysis start time
#'
#' @export
#' @keywords constructor
set_clin <- function(gamma, e_itv, CCOD, CCOD_t, etaC, etaE, d_itv) {
  if (missing(gamma)) stop("Please indicate the rate of enrollment per unit of time (gamma)")

  if (missing(e_itv)) e_itv = NULL
  if (length(e_itv) != length(gamma) & length(e_itv) != length(gamma) - 1) stop("The length of e_itv should be the same length as gamma or 1 less")

  if (missing(CCOD) || CCOD %notin% c("fixed-first", "fixed-last", "event")){
    stop("Please indicate type of clinical cut-off date (CCOD). Options include fixed-first, fixed-last, or event")
  }
  if (missing(CCOD_t)){
    if(CCOD == "fixed-first") stop("Please specify the time difference between analysis start and first patient's enrollment (CCOD_t).")
    else if (CCOD == "fixed-last") stop("Please specify the time difference between analysis start and first last's enrollment (CCOD_t)")
    else stop("Please specify the number of events observed when analysis starts (CCOD_t)")
  }

  if (missing(etaC)) stop("Please specify dropout rate per unit time for control arm (etaC).")
  if (missing(etaE) || length(etaC) != length(etaE)) {
    message("Dropout rate per unit time for treatment arm (etaE) is not specified correcly. Disregard this warning if this is for external trial. Otherwise the same dropout rate as eta is used.")
    etaE = etaC
  }
  if(missing(d_itv)) d_itv = NULL
  if (length(d_itv) < length(etaC) - 1) stop("Please adjust duration of time periods for dropping out of the study (d_itv). The length of d_itv should be same length as etaC or 1 less.")
  new(".clinClass", gamma =  gamma, e_itv = e_itv, CCOD = CCOD,  CCOD_t =  CCOD_t, etaE = etaE, etaC = etaC, d_itv = d_itv)
}

#' S4 Class for setting parameters for time-to-events
#'
#' @keywords class
.eventClass = setClass(".eventClass", slots = list(event = "character",
                                                   lambdaC = "numeric", shape = "numericORNULL", t_itv = "numericORNULL",
                                                   beta = "numericORNULL",
                                                   change = "listORNULL", keep = "charORNULL"))



#' Set up time-to-events
#'
#' @param event Distribution of time-to-events. Piece-wise exponential distribution for \code{event = "pwexp"} and weibull distribution for \code{event = "weibull"}
#' @param lambdaC Baseline hazard rate of internal control arm. Specify a vector for piece-wise hazard with duration specified in \code{t_itv} if \code{event = "pwexp"}
#' @param beta covariates' coefficients. \code{NULL} if no covariates are found in \code{ssObj}
#' @param shape shape of weibull distribution if \code{event = "weibull"}. \code{NULL} if \code{event = "pwexp"}
#' @param t_itv a vector indicating interval lengths where the exponential rates provided in \code{lambdaC} apply. Note that the length of \code{t_itv} is at least 1 less than that of \code{lambdaC} and that the final value rate in \code{lambdaC} applies after time `sum(t_itv)`. \code{NULL} if \code{event = "weibull"}
#' @param change a vector of covariates to be raised to the power specified in \code{power}
#' @param keep wer to raise to for covariates specified in \code{cov}. Length of \code{power} should be 1 or the same as that of \code{cov}
#' @return a \code{.eventClass} class containing time-to-events information
#'
#' @export
#' @keywords constructor
set_event <- function(event, lambdaC, beta, shape, t_itv, change, keep) {

  if (missing(event) || event %notin% c("weibull", "pwexp")) stop("Distribution of time-to-events (event) is not correctly specify. Options include weibull, and pwexp")
  if (missing(lambdaC)) stop ("Please provide the baseline hazard rate of internal control arm (lambdaC).")
  if (missing(t_itv)) t_itv = NULL
  if (missing(shape)) shape = NULL
  if (event == "weibull" & is.null(shape)) stop("Simulate time following weibull distribution. Please provide shape of the weibull distribution")
  if (event == "pwexp" & (length(t_itv) < length(lambdaC) - 1)) stop("Length of t_it should be at least 1 less than that of lambdaC")

  if (missing(change)) change = NULL
  if (missing(keep)) keep = NULL
  if(missing(beta)) beta = NULL

  new(".eventClass", event = event,
      lambdaC = lambdaC, shape = shape, t_itv = t_itv,
      beta = beta, change = change, keep = keep)
}


#' Simulate survival times
#'
#' @param eventObj a \code{.eventClass} class generated with \code{\link{set_event}} containing time-to-events information
#' @param clinInt a \code{.clinClass} class generated with \code{\link{set_clin}}. It contains information on enrollment time, drop-out pattern and analysis start time for the internal trial
#' @param clinExt a \code{.clinClass} class generated with \code{\link{set_clin}}. It contains information on enrollment time, drop-out pattern and analysis start time for the external trial
#' @param seed the seed of R‘s random number generator. Default is 47
#' @return a \code{matrix} containing survival times, censoring, external trial indicator, treatment arm indicator, and covariates
#'
#' @export
#' @keywords methods
setGeneric(name="add_time", def=function(dt, eventObj, clinInt, clinExt, seed){standardGeneric("add_time")})
setMethod(f="add_time", signature(dt = "matrix", eventObj = ".eventClass", clinInt = ".clinClass", clinExt = ".clinClass"),
          definition=function(dt, eventObj, clinInt, clinExt, seed){
            if (missing(dt))  stop("Please provide dt")
            if (missing(eventObj))  stop("Please provide eventObj")
            if (missing(clinInt)) stop("Please provide clinInt.")
            if (missing(clinExt)) stop("Please provide clinExt.")
            if (missing(seed)){
              message("Set.seed(47)")
              seed = 47
            }

            # checl
            if (c("driftHR", "HR", "trt", "ext","extumatch") %notin% colnames(dt) || sum(!grepl("cov[1234567890]", colnames(dt)[colnames(dt) %notin% c("driftHR", "HR", "trt", "ext","extumatch")])) > 0) stop("Please make sure the trial data contains the correct variable names.")

            new_cov_name = NULL
            cov_name = colnames(dt)[colnames(dt) %notin% c("driftHR", "HR", "trt", "ext","extumatch")]
            flog.debug(cat("[add_time] Number of original covariates", length(cov_name), "\n"))

            keep = eventObj@keep
            if (is.null(keep)) {
              keep = cov_name
              message(cat("All original covariates (if any):", keep, "are used for time-to-failure."))
            } else if (sum(grepl("none", keep)) > 0){
              keep = NULL
              message(cat("No original covariates are used for time-to-failure."))
            } else if (sum(keep %notin% cov_name) > 0) {
              stop("Please correctly specify the covariates to keep when simulating failure times.")
            }

            change = eventObj@change
            if (!is.null(change)){ # no change to covariates
              for (i in 1:length(change)){
                if (length(change[[i]]) != 3) {
                  stop(paste("Please correctly specify the operation for the", i, "th entry."))
                } else if (change[[i]][2] == "^") {
                  s = (change[[i]][1] %notin% colnames(dt)[colnames(dt) %notin% c("driftHR", "HR", "trt", "ext","extumatch")]) +
                    is.na(as.numeric(change[[i]][3]))
                  if (s > 0) stop(paste("Please correctly specify the ^ operation for the", i, "th entry."))
                } else if (change[[i]][2] %in% c("+", "*")) {
                  s = (change[[i]][1] %notin% colnames(dt)[colnames(dt) %notin% c("driftHR", "HR", "trt", "ext","extumatch")]) +
                    (change[[i]][3] %notin% colnames(dt)[colnames(dt) %notin% c("driftHR", "HR", "trt", "ext","extumatch")])
                  if (s > 0) stop(paste("Please correctly specify the +/* operation for the", i, "th entry."))
                } else {
                  stop("Please correctly specify the operation. The second entry should be ^, + or *.")
                }
                new_cov_name = c(new_cov_name, paste(change[[i]],collapse=""))
              }
            }

            flog.debug(cat("[add_time] Number of new covariates", length(new_cov_name), "\n"))
            flog.debug(cat("[add_time] Number of original covariates to keep", length(keep), "\n"))

            # flog.debug(cat("[add_time] Current dataset include origin covariates", cov_name, "and new covariates", new_cov_name, "\n"))

            # n_cov = sum(grepl("cov", colnames(dt)))
            n_cov = length(new_cov_name) + length(keep)
            flog.debug(cat("[add_time] Number of covariates is (length of beta should be)", n_cov, "\n"))

            beta = eventObj@beta
            if (length(beta) %notin% c(1, n_cov)){
              message("Coefficient for covariates (beta) is not recognized or correctly specified. Default value 1 is used for all covariates")
              beta = rep(1, n_cov)
            } else if (length(beta) == 1){
              message("User provides one coefficient for covariate. This value is used for all covariates")
              beta = rep(beta, n_cov)
            }


            hr = unique(dt[, 'HR'])
            dr = unique(dt[, 'driftHR'])

            time_int = g_one_t(ext = 0, dt = dt,
                               gamma = clinInt@gamma, e_itv = clinInt@e_itv,
                               etaC = clinInt@etaC, etaE = clinInt@etaE, d_itv = clinInt@d_itv,
                               CCOD = clinInt@CCOD, CCOD_t = clinInt@CCOD_t,
                               event = eventObj@event, lambda = eventObj@lambdaC,
                               shape = eventObj@shape, t_itv = eventObj@t_itv,
                               HR = hr, beta = beta,
                               change = change, keep = keep)

            lambdaExt = eventObj@lambdaC * dr

            time_ext = g_one_t(ext = 1, dt = dt,
                               gamma = clinExt@gamma, e_itv = clinExt@e_itv,
                               etaC = clinExt@etaC, etaE = clinExt@etaE, d_itv = clinExt@d_itv,
                               CCOD = clinExt@CCOD, CCOD_t = clinExt@CCOD_t,
                               event = eventObj@event, lambda = lambdaExt,
                               shape = eventObj@shape, t_itv = eventObj@t_itv,
                               HR = hr, beta = beta,
                               change = change, keep = keep)

            rbind(time_int, time_ext)

          })
