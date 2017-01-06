gradient_G <- function(lambda_G1, lambda_G2, prior, info){
  
  # sufficient stats
  n_ij <- 2*info$stat$n_ij
  s_j <- prior$Q$Sg_cond %>% 
    rep(info$stat$I) %>% matrix(nrow = info$stat$J)
  
  cond_mu <- info$stat$cond_mu
  TxSigma <- info$stat$TxSigma
  
  # 
  grad_G1  <- 
    (2 * n_ij + 1 - lambda_G1) * trigamma(lambda_G1) - 
    (TxSigma - 0.5/s_j) * (1 + 2*lambda_G1) / (lambda_G2^2) +
    cond_mu/(s_j * lambda_G2) + 1
  
  grad_G2  <- 
    -(2 * n_ij + 1)/lambda_G2 - 
    (2 * TxSigma - 1/s_j) * (lambda_G1 + lambda_G1^2) / (lambda_G2^3) -
    (cond_mu/s_j) * (lambda_G1/(lambda_G2^2))
  
  # return
  list(G1 = as.matrix(grad_G1), 
       G2 = as.matrix(grad_G2))
}

ELBO_G <- function(lambda_G1, lambda_G2, prior, info){
  # sufficient stats
  n_ij <- 2*info$stat$n_ij
  s_j <- prior$Q$Sg_cond %>% 
    rep(info$stat$I) %>% matrix(nrow = info$stat$J)
  
  cond_mu <- info$stat$cond_mu
  TxSigma <- info$stat$TxSigma
  
  # 
  2 * n_ij * (digamma(lambda_G1) - log(lambda_G2)) - 
    (TxSigma - 0.5/s_j) * (lambda_G1 + lambda_G1^2)/(lambda_G2^2) + 
    cond_mu*lambda_G1/(s_j * lambda_G2) - 
    ((lambda_G1 - 1)*digamma(lambda_G1) + log(lambda_G2) - 
       lgamma(lambda_G1) - lambda_G1)
}

ELBO_G_vec <- function(lambda, i, j, prior, info){
  lambda_G1 <- lambda[1]
  lambda_G2 <- lambda[2]
  
  # sufficient stats
  n_ij <- 2*info$stat$n_ij
  s_j <- prior$Q$Sg_cond %>% 
    rep(info$stat$I) %>% matrix(nrow = info$stat$J)
  
  cond_mu <- info$stat$cond_mu
  TxSigma <- info$stat$TxSigma
  
  # 
  2 * n_ij * (digamma(lambda_G1) - log(lambda_G2)) - 
    (TxSigma - 0.5/s_j) * (lambda_G1 + lambda_G1^2)/(lambda_G2^2) + 
    cond_mu*lambda_G1/(s_j * lambda_G2) - 
    ((lambda_G1 - 1)*digamma(lambda_G1) + log(lambda_G2) - 
       lgamma(lambda_G1) - lambda_G1)
}