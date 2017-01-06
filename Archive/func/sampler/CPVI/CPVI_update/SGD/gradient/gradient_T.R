gradient_T <- 
  function(par, prior, info, mode, verbose)
  {
    # DOCUMENTATION:
    # Stochastic Gradient for variational parameter:
    #   T_j ~ Gamma(n^j, lambda_{T, j})
    # with fomula E{ d_T[logp - logq] * d_l T(u) }
    
    # INPUT:
    # par$marginal$T: lambda_T paramteres to update
    # par$n_sample: number of gradient samples
    # mode: type of parameter to update: marginal or copula
    
    # RETURN:
    # par$marginal$T: updated lambda_T (J x 1 vector)
    
    if (verbose) cat("T..")
    
    #### 0. Define global variable ####
    # extract relevant parameter
    J <- info$stat$J
    n_j <- info$stat$n_j
    n_sample <- info$oper$n_sample

    var_par <- par[[mode]]$T
    # initialize container
    gradient_list <- 
      array(NaN, dim = c(n_sample, J))
    
    #### 1. Compute Gradient Samples
    for (i in 1:n_sample) {
      info$sample$u <- runif(J)
      info$sample$T <- 
        genfunc_T(var_par, info$sample$u, n_j)
      
      #### 1.1 Compute d_l T(u) ####
      # compute noisy implicit gradient
      deriv_gen <- 
        genfunc_deriv_T(var_par, info, prior)
      
      #### 1.2. Compute d_T [logp - logq] ####
      deriv_diff <- 
        par$marginal$T - 
        colSums(info$sample$sigma * (info$sample$Q)^2)
      
      #### 1.3. Calc product, store ####
      gradient_list[i, ] <- 
        deriv_diff * deriv_gen
    }
    
    #### 3. Average n Output ####
    colMeans(gradient_list)
  }