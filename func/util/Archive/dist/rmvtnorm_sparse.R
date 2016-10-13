library(mvtnorm)

rmvtnorm_sparse <- 
  function(n_OTU, n_sample, 
           Sigma = NULL,
           n_block = 5, theta = 0.01){
    # either provide Sigma, or generate a sparse Sigma
    if (is.null(Sigma)){
      Sigma <- 
        cov_sparse(
          n_block, rep(n_sample/n_block, n_block), 
          theta) %>% abs
    }
    
    # generate mtv normal
    Q <- rmvnorm(n_OTU, rep(0, n_sample), Sigma)
    list(Q = Q, Sigma = Sigma)
  }