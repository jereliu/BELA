gradient_sigma <- 
  function(par, prior, info, mode, verbose)
  {
    # DOCUMENTATION:
    # Stochastic Gradient for variational parameter:
    #   sigma_i ~ Beta(a + n^i + lambda, b)
    # with fomula E{ d_T[logp - logq] * d_l T(u) }
    
    # INPUT:
    # par$marginal$sigma: lambda_sigma paramteres to update
    # par$n_sample: number of gradient samples
    # mode: type of parameter to update: marginal or copula
    
    # RETURN:
    # par$marginal$sigma: updated lambda_sigma (J x 1 vector)
    
    if (verbose) cat("sigma..")
    
    #### 0. Define global variable ####
    # extract relevant parameter
    I <- info$stat$I
    n_i <- info$stat$n_i
    n_sample <- info$oper$n_sample
    a <- prior$sigma$a
    b <- prior$sigma$b
    
    var_par <- par[[mode]]$sigma
    # initialize container
    gradient_list <-
      array(NaN, dim = c(n_sample, I))
    
    #### 1. Compute Gradient Samples
    for (i in 1:n_sample) {
      #### 1.1 Compute d_l T(u) ####
      info$sample$u <- runif(I)
      info$sample$sigma <-
        genfunc_sigma(
          var_par, info$sample$u, n_i, a, b)
      # compute noisy implicit gradient
      deriv_gen <- 
        genfunc_deriv_sigma(var_par, info, prior)
      
      #### 1.2. Compute d_T [logp - logq] ####
      deriv_diff <-
        par$marginal$sigma / info$sample$sigma - 
        rowSums((info$sample$Q)^2 * info$sample$T)
      
      #### 1.3. Calc product, store ####
      gradient_list[i, ] <- 
        deriv_diff * deriv_gen
    }
    
    #### 3. Average n Output ####
    colMeans(gradient_list)
  }