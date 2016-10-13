rv_coef <- function(X, Y){
  # center column vector
  X <- t(t(X) - colMeans(X))
  Y <- t(t(Y) - colMeans(Y))
  
  S_xy <- t(X) %*% Y
  S_xx <- t(X) %*% (X)
  S_yy <- t(Y) %*% (Y) 
  
  COVV <- Trace(S_xy %*% t(S_xy))
  VAV_X <- Trace(S_xx %*% S_xx) 
  VAV_Y <- Trace(S_yy %*% S_yy) 
  
  COVV/sqrt(VAV_X * VAV_Y)
}

#### extra functions ####
nuclear <- function(X) sum(svd(X)$d)