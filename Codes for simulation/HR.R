
#######################################################################################
#
#                                                                                     
#   Filename    :	HR.R  
#   Description : Simulation for right-censored time-to-event outcome 
#                                                                                     
#   Project     :       Biometrical Journal article "Making apples from oranges: comparing non-collapsible effect estimators
#                                      and their standard errors after adjustment for different covariate sets"                                                             
#   Authors     :       Rhian Daniel, Jingjing Zhang, Daniel Farewell                                                                
#   Date        :       26.09.2020
#   Purpose     :       As described in BiomJ article
#																				  
#   R Version   :       3.5.2 (2018-12-20)                                                                
#
#
#   Input data files  :    --                                                           
#   Output data files :    "HR_N_1000_Coefficient_0_1.tiff", "HR_N_1000_Coefficient_0_1.csv",
#                          "HR_N_1000_Coefficient_1_0.tiff", "HR_N_1000_Coefficient_1_0.csv"......
#
#   Required R packages :  simsurv, survival, ggplot2, plyr, ipw, survey, RColorBrewer,...
#
#
########################################################################################
f.HR <- function(iConfounding, m, N_repeated, N, file_path, Coefficient){

  ######
  results_coef <- matrix(0, N_repeated + 5, 4)
  colnames(results_coef) <- c("Unadjusted", "IPTW", "Adjusted marginal", "Conditional")
  status_new_all <- list()
  
  set.seed(1500)
  t_0_repeated <- matrix(runif(N * N_repeated, min = 0, max = 2), N_repeated, N)
  c1 <- rnorm(N)#Covariate c1
  
  for (sim in 1:N_repeated)
  {
    print(sim)
    
    ######################Simulate initial survival data
    set.seed(sim)
    
    if (iConfounding){
      x <- rbinom(N, 1, 1 / (1 + exp(-1 * c1)))#confounding
    } else{
      x <- rbinom(N, 1, 0.5)
    }
    
    xc <- data.frame( x = x, c1 = c1)
    d <- as.data.frame(cbind(simsurv(lambdas = 0.1, gammas = 1.5, 
                                     betas = c(x = Coefficient[1], c1 = Coefficient[2]), 
                                     x = xc), xc))
    ##add censoring
    t_0 <- t_0_repeated[sim,]
    t_1 <- d$eventtime
    t_e <- t_0 + t_1
    
    t <- pmin(t_e, 10) - t_0
    status_new <- as.numeric(t_e <= 10)
    d$status <- status_new
    status_new_all[[sim]] <- status_new
    #########################Conditional & Unadjusted HR 
    cox.fit <- coxph(Surv(eventtime, status) ~ x + c1 , d) #Conditional HR
    results_coef[sim,4] <- cox.fit$coefficients[1]
    bh <- basehaz(cox.fit, centered = FALSE)
    
    cox.fit.unadj <- coxph(Surv(eventtime, status) ~ x, d) #Unadjusted HR
    results_coef[sim, 1]<-cox.fit.unadj$coef
    ###########################IPW
    d$y <- Surv(d$eventtime, d$status)
    
    temp <- ipwpoint(exposure = x, family = "binomial", 
                     link = "logit",numerator = ~ 1, 
                     denominator = ~ c1, data = d)
    
    d$sw <- temp$ipw.weights
    msm <- svycoxph(y ~ x, 
                    design = svydesign(id=~ 1, weights = ~ sw, data = d))
    results_coef[sim, 2] <- msm$coefficients
    ####################Adjust Marginal survival curve
    s1 <- matrix(0, nrow = nrow(bh), ncol = N)#all exposure, i,e, when all x==1 
    s0 <- matrix(0, nrow = nrow(bh), ncol = N)#all non exposure, i,e, when all x==0
    
    for (i in 1:N) {
      s1[, i] <- exp(-bh$hazard) ^ (exp(cox.fit$coef[1] + cox.fit$coef[2] * d$c1[i]))
      s0[, i] <- exp(-bh$hazard) ^ (exp(cox.fit$coef[2] * d$c1[i]))
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
    
    
    #############################Simulate Data based on Adjust Marginal survival curve
    N_1 <- m#round(N/2)#user defined from 500, 2500, 5000
    N_0 <- m#N-N_1
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
    results_coef[sim, 3]<-cox.fit.new$coef
  }
  
  results_coef[N_repeated + 2, ]<-colMeans(results_coef[1:N_repeated, ])
  results_coef[N_repeated + 3, ]<-apply(results_coef[1:N_repeated, ], 2, sd)
  results_coef[N_repeated + 4, ]<-results_coef[N_repeated + 3, ] / sqrt(N_repeated)
  results_coef[N_repeated + 5, ]<-results_coef[N_repeated + 3, ] / sqrt(2 * (N_repeated - 1))
  ######################################Visualization
  dfresults <- data.frame(Methods = 
                            factor(rep(c("Unadjusted", "IPTW", "Adjusted marginal", "Conditional"), each = N_repeated),
                                   levels=c("Unadjusted", "IPTW", "Adjusted marginal", "Conditional")),
                          Value = c(results_coef[c(1:N_repeated), 1],results_coef[c(1:N_repeated), 2],
                                    results_coef[c(1:N_repeated), 3],results_coef[c(1:N_repeated), 4]))
  
  mu <- ddply(dfresults, "Methods", summarise, grp.mean = mean(Value))
  
  plot_title <- c()
  file_name <- c()
  if (iConfounding){
    plot_title <- "Treatment effect, Confounder Effect (Coefficients:(1,1)*)"
    file_name <- paste("HR_", "Nsim_", N_repeated, "_2m_", 2 * m, "_Coefficient_",
                       Coefficient[1], "_", Coefficient[2], "_confounding", sep = "")
  }else if (identical(Coefficient, c(0, 1))){
    plot_title <- "Null treatment effect, Covariate Effect (Coefficients:(0,1))"
    file_name <- paste("HR_", "Nsim_", N_repeated, "_2m_", 2 * m,"_Coefficient_",
                       Coefficient[1], "_", Coefficient[2], sep = "")
  }else if (identical(Coefficient, c(1, 0))){
    plot_title <- "Treatment effect, Null covariate Effect (Coefficients:(1,0))"
    file_name <- paste("HR_", "Nsim_", N_repeated, "_2m_", 2 * m, "_Coefficient_",
                       Coefficient[1], "_", Coefficient[2], sep = "")
  }else{
    plot_title <- "Treatment effect, Covariate Effect (Coefficients:(1,1))"
    file_name <- paste("HR_", "Nsim_", N_repeated, "_2m_", 2 * m, "_Coefficient_",
                       Coefficient[1], "_", Coefficient[2], sep = "")
  }
  
  
  p<-ggplot(dfresults, aes(x = Value, color = Methods)) + 
    geom_density() +
    geom_vline(data = mu, aes(xintercept = grp.mean, color = Methods), linetype="dashed") + 
    scale_color_manual(values = brewer.pal(n = 8, name = "Dark2")[c(1,3,4,6)]) +
    labs(title = plot_title, x="log(HR)", y = "Density")+
    theme(legend.position = "top")
  

  ggsave(paste(file_path, "/results/", file_name, ".tiff", sep = ""), 
         units="in", width=8, height=6, dpi=300, compression = 'none')
  
  ################Save results
  write.table(results_coef, 
              file=paste(file_path, "/results/", file_name, ".csv", sep=""), 
              row.names=FALSE, col.names=TRUE, sep = ",")
  

  
  
  
}