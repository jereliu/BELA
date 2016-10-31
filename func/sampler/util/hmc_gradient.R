grad_hmc_U <- function(T, U, V, lambda, d1){
  (-T + d1(U %*% t(V))) %*% V + 2*lambda * U
}

grad_hmc_V <- function(T, U, V, lambda, d1){
  t(-T + d1(U %*% t(V))) %*% U + 2*lambda * V
}