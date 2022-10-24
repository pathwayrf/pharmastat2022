## adjusted marginal with bootstrap

library(survival)

ad_est <- function(out){
  mean_HR <- out[1]
  reject <- unname(ifelse(out[4]<1,1,0))
  c("reject"=reject,"mean_HR"=mean_HR,"sd_HR"=out[2],"p_value"=out[5])
}


HR_new <- function(data, indices, N_new_half, mcmc_list, pred){#N_new_half, survival times be simulated, i.e. m in the paper
  
  dt <- data
  
  df <- summary(mcmc_list)
  coef <- df[[1]][grep("beta",rownames(df[[1]])),1]
  r0 <- df[[1]][grep("r0",rownames(df[[1]])),1]
  
  ## baseline hazard
  dt.int <- data.frame(dt) #RCT
  dt.int <- dt.int[dt.int$ext==0,]
  surv.formula <- paste("Surv(time, 1-cnsr) ~ trt +", paste0(pred, collapse="+"))
  #print(surv.formula)
  cox.fit <- coxph(as.formula(surv.formula), data = dt.int) #Conditional HR
  bh <- basehaz(cox.fit, centered = FALSE)
  
  #bh$hazard <- sapply(bh$time,f.h,r0=r0)
  
  ####################Adjust Marginal survival curve
  d <- as.matrix(dt.int[,pred])
  N <- dim(d)[1]
  
  s1 <- matrix(0, nrow = nrow(bh), ncol = N)#all exposure, i,e, when all x==1 
  s0 <- matrix(0, nrow = nrow(bh), ncol = N)#all non exposure, i,e, when all x==0
  
  for (i in 1:N) {
    s1[, i] <- exp(-bh$hazard) ^ (exp(coef[1] + d[i,] %*% coef[-1]))
    s0[, i] <- exp(-bh$hazard) ^ (exp(d[i,] %*% coef[-1]))
  }
  
  s1_cond <- vector(mode="numeric", length = nrow(bh))
  s0_cond <- vector(mode="numeric", length = nrow(bh))
  
  s1_Mean <- rowMeans(s1)
  s0_Mean <- rowMeans(s0)
  
  s1_cond[1] <- s1_Mean[1]
  s0_cond[1] <- s0_Mean[1]
  
  for (j in 2:nrow(bh)) {
    s1_cond[j] <- s1_Mean[j] / s1_Mean[j - 1]
    s0_cond[j] <- s0_Mean[j] / s0_Mean[j - 1]
  }
  
  #############################Simulated Data 
  N_1 <- N_new_half#100#round(N/2),100,500,1000,2500,5000,10000,25000
  N_0 <- N_new_half#100#N-N_1,100,500,1000,2500,5000,10000
  ######Exposure
  #eventtime of dead
  newd1 <- matrix(0, nrow = nrow(bh), ncol = N_1)
  ######Exposure
  newd1=matrix(0, nrow = nrow(bh), ncol = N_1)
  for (i in 1:N_1) {
    newd1[,i] <- rbinom(nrow(bh), 1, s1_cond)
  }
  
  index_event1 <- apply(newd1, 2, function(x) which(x == 0)[1])
  
  eventtime1 <- bh$time[index_event1]
  eventtime1[which(is.na(index_event1))] <- max(bh$time)
  status1 <- rep(1, N_1)
  status1[which(is.na(index_event1))] <- 0
  ########UnExposure
  newd0=matrix(0, nrow = nrow(bh), ncol = N_0)
  for (i in 1:N_0) {
    newd0[,i] <- rbinom(nrow(bh), 1, s0_cond)
  }
  
  index_event0 <- apply(newd0, 2, function(x) which(x == 0)[1])
  
  eventtime0 <- bh$time[index_event0]
  eventtime0[which(is.na(index_event0))] <- max(bh$time)
  status0 <- rep(1, N_0)
  status0[which(is.na(index_event0))] <- 0
  #####New data
  newdata <- data.frame(rbind(cbind(eventtime1, status1, rep(1, N_1)), 
                              cbind(eventtime0, status0, rep(0, N_0))))
  names(newdata) <- c("eventtime", "status", "x")
  cox.fit.new <- coxph(Surv(eventtime, status) ~ x, newdata) 
  ##########Results of "AdjustMarginal" HR
  return(summary(cox.fit.new)$coefficient[2])
}

### Calculate adjusted marginal HR
adjust_margin <- function(mcmc_list,dt,pred,m=500){
  
  #############################Simulate Data based on Adjust Marginal survival curve
  HR_new.boot <- boot(dt, HR_new, R = 1000, N_new_half = m, mcmc_list=mcmc_list, pred=pred)
  HR_new.boot.ci <- boot.ci(HR_new.boot, conf = 0.95, type = c("perc"))
  
  #### Convert confidence interval to a p-value
  lower <- HR_new.boot.ci$percent[,4]
  upper <- HR_new.boot.ci$percent[,5]
  z <- mean(HR_new.boot$t-1)*(2*1.96)/(upper-lower)
  p.value <- pnorm(z)
  
  c(mean(HR_new.boot$t),sd(HR_new.boot$t),lower, upper, p.value)
}
