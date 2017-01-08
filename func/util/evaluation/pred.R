library(abind)

predMeanError <- function(rec, true_theta, num = 100){
  # use the first half as burn in
  num <- min(dim(rec$U)[1], num)
  iter_max <- dim(rec$Theta)[1] - 1
  iter <- 
    seq(iter_max * 0.25, iter_max, length.out = num) %>% round
  
  # evaluate effective sample size/total sample size
  if (is.null(rec$Theta)){ # calculate Theta if not exist
    Theta <- 
      lapply(1:dim(rec$U)[1], 
             function(i) rec$U[i, ,] %*% t(rec$V[i, , ])) %>% abind(along = 0)
  } else {
    Theta <- rec$Theta
  }
  
  mean_Pred <- 
    sapply(2:length(iter), function(i)
      sqrt(mean((apply(Theta[iter[1]:iter[i], , ], 2:3, mean) - true_theta)^2)/
        var(as.vector(true_theta))
        )
    )
  
  time_Pred <- rec$time[iter[-1]]

  cbind(time_Pred - min(time_Pred), mean_Pred)
}
