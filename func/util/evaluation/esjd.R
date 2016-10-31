ESJD <- function(rec, config){
  # evaluate effective sample size/total sample size
  target <- rec$U[, 1, 1]
  
  sjd <- (target[-1] - target[-length(target)])^2
  msjd <- cumsum(sjd)/1:length(sjd)
  
  if (!is.null(config$sampler$frog_step)){
    msjd <- msjd/sqrt(config$sampler$frog_step)
  }
    
  msjd
}
