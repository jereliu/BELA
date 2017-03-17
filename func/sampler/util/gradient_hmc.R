# conditional gradient for HMC

grad_hmc_U <- function(T_suff, U, V, lambda, d1){
  (-T_suff + d1(U %*% t(V))) %*% V + lambda * U
}

grad_hmc_V <- function(T_suff, U, V, lambda, d1){
  t(-T_suff + d1(U %*% t(V))) %*% U + lambda * V
}

lik_target <- function(U, V, T_suff, lambda, dist_family){
  dist_family$negloglik(T_suff, U %*% t(V)) + 
    (lambda/2) * (sum(U^2) + sum(V^2))
}

lik_target_S <- function(S, T_suff, lambda, dist_family){
  U <- S[1:nrow(T_suff), ]
  V <- S[(nrow(T_suff)+1):nrow(S), ]
  
  dist_family$negloglik(T_suff, U %*% t(V)) + 
    (lambda/2) * (sum(U^2) + sum(V^2))
}

