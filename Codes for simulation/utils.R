# Helper functions
# Last update on 082320

#' @import dplyr
#' @import data.table
#' @import rjags
#' @import mvtnorm
#' @import futile.logger



# Specify sample size
s_trt = function(ssC, ssE, ssExt){
  ext = cbind("ext" = rep(1, ssExt), "trt" = rep(0, ssExt))
  # int = cbind("ext" = rep(0, ssC + ssE), "trt" = sample(c(rep(1, ssE), c(rep(0, ssC)))))
  int = cbind("ext" = rep(0, ssC + ssE), "trt" = c(rep(1, ssE), c(rep(0, ssC))))
  rbind(ext, int)
}

# Simulate covariates
s_cov = function(ext, dt, n_cat, n_cont, mu, var, cov, prob){
  n = sum(dt[,'ext'] == ext)
  cov_st = sum(grepl("cov", colnames(dt))) + 1

  covariate <- NULL
  n_cov = sum(n_cont, n_cat)

  if (n_cat == 0 & n_cont == 1){
    covariate <- as.matrix(rnorm(n, mean = mu, sd = sqrt(var)), ncol = 1)
  } else if (n_cont == 0 & n_cat > 0 & sum(cov) == 0){ # if user want to simulate independent binary variables only
    covariate = sapply(1:n_cat, function(i){
      rbinom(n, 1, prob[i])
    })
  } else {
    sig <- diag(var)
    sig[lower.tri(sig)] <- cov
    sig[upper.tri(sig)] = t(sig)[upper.tri(sig)] # create a symmetrical matrix

    covariate <- rmvnorm(n, mean = mu, sigma = sig)
    if (n_cat > 0) {
      covariate[,1:n_cat] = sapply(1:n_cat, function(i){
        cut_val = qnorm(prob[i], mean = mu[i], sd = sqrt(var[i]))

        (covariate[, i] > cut_val)
      })
    }
  }

  colnames(covariate) = paste0("cov",cov_st:(cov_st + n_cov - 1))
  covariate
}

c_cov = function(dt, change, keep){
  dt2 <- dt
  new_cov_name = NULL
  for (i in 1:length(change)){
    if (change[[i]][2] == "^") {
      new_cov <- dt[, change[[i]][1]] ^ as.numeric(change[[i]][3])
    } else if (change[[i]][2] == "*") {
      new_cov <- dt[, change[[i]][1]] * dt[, change[[i]][3]]
    } else if (change[[i]][2] == "+") {
      new_cov <- dt[, change[[i]][1]] + dt[, change[[i]][3]]
    }
    dt2 <- cbind(new_cov, dt2)
    colnames(dt2)[1]<- paste(change[[i]],collapse="")
    new_cov_name = c(new_cov_name, paste(change[[i]],collapse=""))
  }
  dt3 <- dt2[, c(keep, new_cov_name, "trt", "ext", "extumatch"), drop=FALSE]
  dt3
}


# Simulate time-to-events following weibull distribution
s_wb_t = function(dt, lambda, HR, beta, shape){
  n = nrow(dt)
  n_cov = sum(grepl("cov", colnames(dt))) #check if the data include any covariate
  cov = dt[, colnames(dt) %notin% c("trt", "ext", "extumatch"), drop = FALSE]
  
  if (n_cov == 0) {
    log_k = - (log(lambda) + log(HR) * dt[, "trt"])
  } else {
    log_k = - (log(lambda) + log(HR) * dt[, "trt"] + cov %*% beta)/shape #log(k), k is the scale
  }
  sapply(log_k, function(i){rweibull(1, shape, exp(i))})
}


# Simulate time-to-events following piecewise-exponential distribution
s_exp_t = function(dt, lambda, HR, beta, t_itv){

  n = nrow(dt)
  n_cov = sum(grepl("cov", colnames(dt))) #check if the data include any covariate
  if (n_cov == 0) lambda_x = lambda * exp(log(HR) * dt[, "trt"])
  else lambda_x = lambda * exp(log(HR) * dt[, "trt"] + dt[, colnames(dt) %notin% c("driftHR", "HR", "trt", "ext", "extumatch")] %*% beta)
  sapply(lambda_x, function(i) rpwexp(1, rate = i, intervals = t_itv, cumulative=FALSE))
}


# Simulate enrollment time following piecewise-exponential distribution
s_exp_e = function(n, gamma, e_itv){

  aR <- e_itv
  if (n-sum(gamma[1:length(gamma)-1]*e_itv[1:length(gamma)-1])<=0) {
    stop ("[s_exp_e] User requested to keep the enrollment rate fixed but the total enrollment is already greater than the specified sample size prior to the last interval of R.")
  } else{

    aR[length(gamma)] <- (n - sum(gamma[1:length(gamma)-1]*e_itv[1:length(gamma)-1]))/gamma[length(gamma)]
    sample(rpwexp(n, rate=gamma, intervals = e_itv[1:length(gamma)-1], cumulative=TRUE))
  }

}


# Simulate censored survival times
g_one_t <- function(ext, dt,
                    gamma, e_itv, # enroll
                    event, lambda, HR, beta, shape, t_itv, # event
                    etaC, etaE, d_itv, # drop out
                    CCOD, CCOD_t, # CCOD
                    change, keep){

  d <- dt[dt[, "ext"] == ext,] # matrix

  x <- data.table(d)
  n = nrow(x)
  ssC = sum(x$trt == 0);  ssE = sum(x$trt == 1)

  x[, enterT := s_exp_e(n, gamma, e_itv)] #duration

  # update matrix

  d_change = d[, colnames(d) %notin% c("driftHR", "HR")]
  if (!is.null(change)) d_change = c_cov(d[, colnames(d) %notin% c("driftHR", "HR")], change, keep) #update the matrix

  # survival time
  if (event  == "weibull") x[, time0 := s_wb_t(d_change, lambda, HR, beta, shape)]
  else x[, time0 := s_exp_t(d_change, lambda, HR, beta, t_itv)]

  # adding LTFU
  x[x$"trt" == 0, ltfuT := rpwexp(ssC, etaC, intervals = d_itv, cumulative=FALSE)]
  x[x$"trt" == 1, ltfuT := rpwexp(ssE, etaE, intervals = d_itv, cumulative=FALSE)]
  x[, time1:= ifelse(time0 > ltfuT, ltfuT, time0)]
  x[, cnsr1 := ifelse(time0 > ltfuT, 1, 0)]
  x[, ct := enterT + time1] # ct: calendar time a subject had evt/censoring

  # add analysis start time
  if (CCOD == "fixed-first") { #pre-fixed, calculate from enrollment of first patients
    x[, time := ifelse(ct > CCOD_t, CCOD_t - enterT, x$time1)]
    x[, cnsr := ifelse(ct > CCOD_t, 1, x$cnsr1)]
    x = x[x$time >= 0, ]
  } else if (CCOD == "fixed-last") {
    st = max(x$enterT) + CCOD_t
    x[, time := ifelse(ct > st, st - enterT, x$time1)]
    x[, cnsr := ifelse(ct > st, 1, x$cnsr1)]
    x[, cnsr := ifelse(ct > CCOD_t, 1, x$cnsr1)]
  }
  else { # analysis starts when st events have been observed
    st = sort(x$ct[x$cnsr1 == 0])[pmin(sum(x$cnsr1 == 0), CCOD_t)] # compete leaving with event
    if (identical(st, numeric(0))) stop ("Please decrease the number of events required to start the analysis (CCOD_t)")
    x[, time := ifelse(ct > st, st - enterT, x$time1)]

    x[, cnsr := ifelse(ct > st, 1, x$cnsr1)]

    print(paste("User chooses to let number of events drive the CCOD", sum(x$time < 0),
                "patients were enrolled in the trial after ", st,
                "days and were excluded from the analysis"))

    x = x[x$time >= 0, ]
  }
  x_fin <- x[, -c("enterT", "time0", "ltfuT", "time1", "cnsr1", "ct")]
  as.matrix(sapply(x_fin, as.numeric))
}


# Calculate propensity score
# Function r_post below calls this
c_ps <- function(dt, cov_name){
  dat = as.data.frame(dt)
  dat$int = 1 - dat$ext
  var_col = c(cov_name, "int")
  dat2 = dat[, var_col]
  mod <- glm(int ~ ., family = binomial(), data = dat2)
  cbind(dt, "wgt" = mod$fitted.values)
}


# Generate posterior sample
r_post <- function(dt,
                   prior, pred,
                   r0, alpha, sigma,
                   n.chains, n.adapt, n.iter, seed){
  cov_name = colnames(dt)[colnames(dt) %notin% c("driftHR", "HR", "trt", "ext", "time", "extumatch", "cnsr", "wgt")]

  if (sum(grepl("none", pred)) > 0) {
    x = dt[,c("trt"),drop=FALSE]
    n_cov = 0
    #flog.debug("No covariates are used as predictors in the weibull distribution. Any other input is disregarded.\n")
  } else if (sum(grepl("all", pred)) > 0) {
    x = dt[,c("trt", cov_name)]
    n_cov = length(cov_name)
    #flog.debug(cat("All covariates", cov_name, "are used as predictors in the weibull distribution. Any other input is disregarded. \n"))
  } else if (sum(grepl("ps", pred)) == 0) {
    x = dt[,c("trt", pred)]
    n_cov = length(pred)
    #flog.debug(cat("Selected covariate(s)", pred, "are used as predictors in the weibull distribution. Any other input is disregarded. \n"))
  } else if (sum(grepl("ps", pred)) > 0 & ("wgt" %in% colnames(dt))) { # column wgt already existed
    x = dt[,c("trt", "wgt")]
    n_cov = 1
    #flog.debug("Weight already provided in the data. It is used as predictor directly. \n")
  } else if (sum(grepl("ps", pred)) > 0 & length(pred) == 1){ # column wgt does not exist
    dt_ps <- c_ps(dt, cov_name)
    x = dt_ps[,c("trt", "wgt")]
    n_cov = 1
    #flog.debug(cat("Propensity score is calculated based on all coariates", cov_name, "are used as predictors in the weibull distribution.\n"))
  } else if (sum(grepl("ps", pred)) > 0 & length(pred) > 1){
    dt_ps <- c_ps(dt, pred[pred != "ps"])
    x = dt_ps[,c("trt", "wgt")]
    n_cov = 1
    #flog.debug(cat("Propensity score calculated based selected covariate(s)", pred[pred != "ps"],
    #               "ise used as predictors in the weibull distribution.\n"))
  }

  data_list = list(N = nrow(dt), timev = dt[,'time'], event = 1-dt[,'cnsr'], ext = dt[,'ext'] +1, n_cov = n_cov, x = x)
  inits_list = list(.RNG.name="base::Super-Duper", .RNG.seed = seed, r0 = r0, alpha = alpha, beta = rep(0, n_cov + 1)) #.RNG.seed = 1,

  out_list = c("alpha", "beta", "r0", "HR_trt_cc", "HR_cc_hc")
  
  if (prior == "full_ext") {
    prior = "no_ext"

  } else if (prior == "no_ext") { # no borrow
    dt2 = dt[dt[, "ext"] == 0,] #remove external trial data

    if (sum(grepl("none", pred)) > 0) {
      x = dt2[,c("trt"),drop=FALSE]
    } else if (sum(grepl("all", pred)) > 0) {
      x = dt2[,c("trt", cov_name)]
    } else {
      x = dt2[,c("trt", pred)]
    }

    data_list = list(N = nrow(dt2), timev = dt2[,'time'], event = 1-dt2[,'cnsr'], ext = dt2[,'ext'] +1, n_cov = n_cov, x = x)

  } else if (prior == "cauchy" | prior == "unif") {
    
      dt2 = dt[dt[, "extumatch"] == 0,] #remove unmatched external trial data
    
      if (sum(grepl("none", pred)) > 0) {
        x = dt2[,c("trt"),drop=FALSE]
      } else if (sum(grepl("all", pred)) > 0) {
        x = dt2[,c("trt", cov_name)]
      } else {
        x = dt2[,c("trt", pred)]
      }
      
    data_list = list(N = nrow(dt2), timev = dt2[,'time'], event = 1-dt2[,'cnsr'], ext = dt2[,'ext'] +1, n_cov = n_cov, x = x)
      
      inits_list = list(.RNG.name = "base::Super-Duper", .RNG.seed = seed, r0 = r0, alpha = alpha, beta = rep(0, n_cov + 1), sigma = sigma) #.RNG.seed = 1,
    out_list = c("alpha", "beta", "HR_trt_cc", "r0", "HR_cc_hc", "prec")

  } else if (prior == "gamma") {
     
      dt2 = dt[dt[, "extumatch"] == 0,] #remove unmatched external trial data
    
      if (sum(grepl("none", pred)) > 0) {
        x = dt2[,c("trt"),drop=FALSE]
      } else if (sum(grepl("all", pred)) > 0) {
        x = dt2[,c("trt", cov_name)]
      } else {
        x = dt2[,c("trt", pred)]
      }
      
    data_list = list(N = nrow(dt2), timev = dt2[,'time'], event = 1-dt2[,'cnsr'], ext = dt2[,'ext'] +1, n_cov = n_cov, x = x)

      out_list = c("alpha", "beta", "HR_trt_cc", "r0", "HR_cc_hc", "tau")
  }
  
  prior_txt = paste0("model/",prior,".txt")
  #prior_txt = system.file("model", paste0(prior, ".txt"), package = "psborrow")
  
  jags.out <- jags.model(prior_txt,
                         data = data_list,
                         inits = inits_list,
                         n.chains = n.chains, n.adapt = n.adapt, quiet=FALSE)

  # set.seed(314)
  jags.sample <- coda.samples(jags.out, out_list, n.iter = n.iter)
  print(summary(jags.sample))
  jags.sample
}