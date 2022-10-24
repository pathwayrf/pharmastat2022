library(MatchIt)
library(dplyr)
m_cov <- function(dt, match.formula){
  dt2 <- data.frame(dt)
  dt2$int = 1 - dt2$ext
  
  m.out <- matchit(as.formula(match.formula), method = "nearest", ratio = 1, caliper = 0.2, data = dt2)
  #matchit(ext ~ cov1 + cov2, method = "nearest", ratio = 1, caliper = 0.2, data = sample_cov.df) # match on all covariates
  
  dt.match <- match.data(m.out,drop.unmatched = F)
  
  dt.match$extumatch <- ifelse(dt.match$weights>0 | dt.match$int==1,0,1)  #1=unmatch
  print(paste0(nrow(dt.match[dt.match$extumatch==1,]), " patients from the external arm are excluded."))
  
  match.all <- as.matrix(dt.match[,colnames(dt.match) %notin% c("distance", "weights", "subclass", "int")])
  #print(paste0("After matching, we have ", nrow(match.all), " patients in total."))
  return(match.all)
}


match_cov <- function(dt, match) {
  if (missing(dt)) stop("Please provide a list of simulated data (dt).")
  if (missing(match)) {
    stop("Please provide covariates name to match on.")
  } else if (sum(!match %in% colnames(dt[[1]])[!colnames(dt[[1]]) %in% c("driftHR", "HR", "trt", "ext")]) > 0) {
    stop("Please make sure the trial data contains the correct variable names.")
  }
  match.formula <- paste("int ~", paste0(match, collapse=' + '))
  print(match.formula)
  res_list <- lapply(seq(1, length(dt), by = 1), function(i){
    m_cov(dt[[i]], match.formula)
  }) #loop foreach
  
  res_list
}

