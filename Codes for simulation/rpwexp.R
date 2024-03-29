#' The Piecewise Exponential Distribution
#'
#' The piecewise exponential distribution allows a simple method to specify a distribtuion
#' where the hazard rate changes over time. It is likely to be useful for conditions where
#' failure rates change, but also for simulations where there may be a delayed treatment
#' effect or a treatment effect that that is otherwise changing (e.g., decreasing) over time.
#' \code{rpwexp()} is to support simulation both the Lachin and Foulkes (1986) sample size
#' method for (fixed trial duration) as well as the Kim and Tsiatis(1990) method
#' (fixed enrollment rates and either fixed enrollment duration or fixed minimum follow-up);
#' see \code{\link[gsDesign:nSurv]{gsDesign}}.
#'
#' Using the \code{cumulative=TRUE} option, enrollment times that piecewise constant over
#' time can be generated.
#'
#' @param n Number of observations.
#' If length(n) > 1, the length is taken to be the number required.
#' @param rate Vector containing exponential failure rates in intervals described by
#' \code{intervals}
#' @param intervals Vector containing positive values indicating interval lengths where
#' the exponential rates provided in \code{rate} apply. Note that the length of
#' \code{intervals} is 1 less than that of \code{rate} and that the final value rate
#' in \code{rate} applies after time `sum(intervals)`.
#' @param cumulative \code{FALSE} (the default) generates \code{n} independent,
#' identically distributed piecewise exponential failure rates according to the distribution
#' specified by \code{intervals} and \code{rate}. \code{TRUE} generates independent
#' inter-arrival times with the rates of arrival in each interval specified by
#' \code{intervals} determined by \code{rate}.

rpwexp <- function(n, rate=1, intervals=NULL, cumulative=FALSE){
  if(is.null(intervals)){
    if (cumulative){return(cumsum(rexp(n,rate[1])))}else
      return(rexp(n,rate[1]))}
  k <- length(rate)
  if (k==1){
    if(cumulative){return(cumsum(rexp(n,rate)))}else
      return(rexp(n,rate))
  }
  if (length(intervals) < k-1) stop("length(intervals) must be at least length(rate) - 1")
  tx <- 0
  j <- 1
  times <- array(0,n)
  timex <- cumsum(intervals)
  indx <- array(TRUE,n)
  for(i in 1:k){
    nindx <- sum(indx)
    if (nindx==0) break
    increment <- rexp(nindx,rate[i])
    if (cumulative) times[indx] <- tx + cumsum(increment)
    else times[indx] <- tx + increment
    if (i<k){
      tx <- timex[i]
      indx <- (times > timex[i])
    }
  }
  return(times)
}
