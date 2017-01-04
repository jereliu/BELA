# TODO: Dynamically adjust step size
glrm_optimizer <- 
  function(Y, lambda, family_name, 
           init, config, rec, info)
  {
    #warning("V not updated")
    # unpack family properties
    n <- info$n
    p <- info$p
    k <- info$k
    true_theta <- info$true_theta
    
    family <- glrm_family(family_name)
    T_suff <- family$sufficient(Y)
    d1 <- family$partition$d
    d2 <- family$partition$d2
    negloglik <- family$negloglik
    
    # unpack mcmc parameters
    record_freq <- config$record_freq
    time_max <- config$time_max
    
    iter_max <- config$optimiz$iter_max
    optim_tol <- config$optimiz$optim_tol  
    step_size <- config$optimiz$step_size
    
    U_cur <- init$U
    V_cur <- init$V
    Ru_cur <- matrix(rnorm(n*k), ncol = k)
    Rv_cur <- matrix(rnorm(p*k), ncol = k) 
    
    time0 <- proc.time()[3]
    
    # initiate sampler
    pb <- txtProgressBar(min = 1, max = iter_max, style = 3)
    
    for (iter in 1:iter_max) {
      setTxtProgressBar(pb, iter)
      
      U_old <- U_cur; V_old <- V_cur
      Ru_old <- Ru_cur; Rv_old <- Rv_cur
      Theta_old <- U_old %*% t(V_old)
      
      # generate A', A'', B
      A_d1 <- d1(Theta_old)
      A_d2 <- d2(Theta_old)
      
      # sample for U, then acc/rej
      step_size_true <- 
        step_size * 
        (1/(2 * (apply(A_d2, 1, max) + 2 * lambda))) *
        matrix(1, n, k)
      
      grad <- grad_hmc_U(T_suff, U_cur, V_cur, lambda/2, d1)
      U_cur <- 
        U_cur - 
        0.5 * step_size_true^2 * grad + 
        step_size_true * Ru_cur
      Ru_cur <- 
        Ru_cur - 0.5 * step_size_true * grad
      
      # sample for V, then acc/rej
      step_size_true <-
        step_size *
        t((1/(2 * (apply(A_d2, 2, max) + 2 * lambda))) *
            matrix(1, k, p))

      grad <-
        grad_hmc_V(T_suff, U_cur, V_cur, lambda/2, d1)

      V_cur <-
        V_cur - step_size_true * grad +
        step_size_true * Rv_cur
      Rv_cur <- 
        Rv_cur - 0.5 * step_size_true * grad
      
      
      # record
      if (iter %% record_freq == 0){
        rec$U[iter/record_freq, , ] <- U_cur
        rec$V[iter/record_freq, , ] <- V_cur
        
        # rec$acc[iter/record_freq, ] <- 
        #   c(mean(acc_U), mean(acc_V))
        rec$obj[iter/record_freq] <- 
          negloglik(T_suff, U_cur %*% t(V_cur)) + 
          lambda * (sum(V_cur^2) + sum(U_cur^2))
        
        rec$time[iter/record_freq] <-
          (proc.time()[3] - time0)/60
        
        if (iter/record_freq > 1){
          rec$error[iter/record_freq] <-
            abs(rec$obj[iter/record_freq] - rec$obj[iter/record_freq-1])
          
          if (rec$error[iter/record_freq] < optim_tol){
            cat("\n Convergence! ヽ(　･∀･)ﾉ")
            rec$output <- list(U = U_cur, V = V_cur)
            return(rec)
          }
        }
      }
    }
    
    plot(rec$obj, type = "l")
    # return
    cat(paste0("\n max iteration (", iter_max, ") reached with no convergence T_T"))
    rec$output <- list(U = U_cur, V = V_cur)
    rec
  }