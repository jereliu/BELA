essMatrix <- function(rec, num = 1e2){
  # use the first half as burn in
  iter <- 
    seq(dim(rec$U)[1]/2, dim(rec$U)[1], length.out = num) %>% round
  
  # evaluate effective sample size/total sample size
  ess_U <- sapply(2:length(iter),
                  function(i) effectiveSize(rec$U[iter[1]:iter[i], 1, 1]))
  ess_V <- sapply(2:length(iter),
                  function(i) effectiveSize(rec$V[iter[1]:iter[i], 1, 1]))
  ess_time <- rec$time[iter[-1]]
  
  list(time = ess_time - min(ess_time), U = ess_U, V = ess_V)
}

essMatrix_whole <- function(rec, num = 10){
  # use the first half as burn in
  iter_cut <- 
    seq(dim(rec$U)[1]/2, dim(rec$U)[1], dim(rec$U)[1]/(2*num))
  
  # evaluate effective sample size/total sample size
  ess_U <- numeric(num)
  ess_V <- numeric(num)
  
  pb <- txtProgressBar(
    min = 1, max = num, style = 3)
  
  for (i in 1:(length(iter_cut)-1)){
    setTxtProgressBar(pb, i)
    
    iter_pre <- iter_cut[i]
    iter <- iter_cut[i + 1]
    
    ess_U[i] <- 
      outer(1:dim(rec$U)[2], 1:dim(rec$U)[3],
            Vectorize(function(i, j) 
              effectiveSize(rec$U[iter_pre:iter, i, j])
            )) %>% 
      mean %>% divide_by(iter - iter_pre)
    
    ess_V[i] <- 
      outer(1:dim(rec$V)[2], 1:dim(rec$V)[3],
            Vectorize(function(i, j) 
              effectiveSize(rec$V[iter_pre:iter, i, j])
            )) %>% 
      mean %>% divide_by(iter - iter_pre)
  }
  
  list(U = ess_U, V = ess_V)
}