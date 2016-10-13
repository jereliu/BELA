library(VineCopula)

gradient_Q <- 
  function(par, prior, info, mode, verbose)
  {
    # DOCUMENTATION:
    # Stochastic Gradient for variational parameter:
    #   sigma_i ~ Beta(a + n^i + lambda, b)
    # with fomula E{ d_T[logp - logq] * d_l T(u) }
    
    # INPUT:
    # par$marginal$Q: lambda_Q paramteres to update
    # par$n_sample: number of gradient samples
    # mode: type of parameter to update: marginal or copula
    
    # RETURN:
    # par$marginal$sigma: updated lambda_sigma (J x 1 vector)
    
    if (verbose) cat("Q..")
    
    #### 0. Define global variable ####
    # extract relevant parameter
    I <- info$stat$I
    J <- info$stat$J
    n_i <- info$stat$n_i
    
    n_sample <- info$oper$n_sample
    a <- prior$sigma$a
    b <- prior$sigma$b
    
    var_par <- par[[mode]]$Q
    
    # initialize container
    gradient_list <-
      array(NaN, dim = c(n_sample, I, J))
    
    #### 1. Compute Gradient Samples
    for (i in 1:n_sample) {
      #### 1.1 Compute d_l Q(u) ####
      info$sample$u <- 
        matrix(runif(I*J), nrow = I)
      info$sample$Q <- 
        sapply(1:J,
               function(j) 
                 genfunc_Q(var_par[, j], 
                           info$sample$u[, j]))
      
      # compute noisy implicit gradient
      deriv_gen <- 
        genfunc_deriv_Q(var_par, info)
      
      #### 1.2. Compute d_T [logp - logq] ####
      deriv_diff <- 
        2 * (
          (var_par - 
             info$sample$sigma %*% 
             t(info$sample$T)) * 
            info$sample$Q 
          - 
            info$sample$Q %*% 
            prior$Q$Sigma_inv
        )
      
      #### 1.3. Calc product, store ####
      gradient_list[i, , ] <- 
        deriv_diff * deriv_gen
    }
    
    #### 3. Average n Output ####
    apply(gradient_list, 1, mean)
  }