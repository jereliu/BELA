# TODO: Dynamically adjust step size

glrm_sampler_hmc <- 
  function(Y, lambda, family_name, 
           init, config, rec, info){
    # unpack family properties
    n <- info$n 
    p <- info$p
    k <- info$k
    true_theta <- info$true_theta
    
    family <- glmr_family(family_name)
    T <- family$sufficient(Y)
    d1 <- family$partition$d
    d2 <- family$partition$d2
    negloglik <- family$negloglik
    
    # unpack mcmc parameters
    record_freq <- config$record_freq
    time_max <- config$time_max
    
    iter_max <- config$sampler$iter_max
    step_size <- config$sampler$step_size
    frog_step <- config$sampler$frog_step  
    mmtm_freq <- config$sampler$mmtm_freq 
    rotn_freq <- config$sampler$rotn_freq  
    
    U_cur <- init$U
    V_cur <- init$V
    Ru_cur <- matrix(rnorm(n*k), ncol = k)
    Rv_cur <- matrix(rnorm(p*k), ncol = k) 
    
    time0 <- proc.time()[3]
    
    # initiate sampler
    pb <- txtProgressBar(min = 1, max = iter_max, style = 3)
    
    for (iter in 1:iter_max) {
      setTxtProgressBar(pb, iter)
      
      if (iter %% mmtm_freq == 0) {
        Ru_cur <- matrix(rnorm(n*k), ncol = k)
        Rv_cur <- matrix(rnorm(p*k), ncol = k) 
      }
      
      if (iter %% rotn_freq == 0) {
        R <- rortho(k)
        U_cur <- U_cur %*% R
        V_cur <- V_cur %*% R
      }
      
      U_old <- U_cur; V_old <- V_cur
      Ru_old <- Ru_cur; Rv_old <- Rv_cur
      Theta_old <- U_old %*% t(V_old)
      
      # generate A', A'', B
      A_d1 <- d1(Theta_old)
      A_d2 <- d2(Theta_old)
      
      # sample for U, then acc/rej
      for (f_step in 1:frog_step){
        step_size_true <- 
          step_size * 
          (1/(2 * (apply(A_d2, 1, max) + 2 * lambda))) *
          matrix(1, n, k)
        
        grad <- grad_hmc_U(T, U_cur, V_cur, lambda, d1)
        U_cur <- 
          U_cur - 
          0.5 * step_size_true^2 * grad + 
          step_size_true * Ru_cur
        Ru_cur <- 
          Ru_cur - 0.5 * step_size_true * grad
      }
      
      # metroplis step for U
      # acc_U <- 1
      acc_prob <-
        (negloglik(T, U_old %*% t(V_old)) - 
           negloglik(T, U_cur %*% t(V_old)) +
           lambda *
           (sum(U_old^2) - sum(U_cur^2)) +
           0.5 *
           (sum(Ru_old^2) - sum(Ru_cur^2))
        ) %>% min(0, .) %>% exp
      
      acc_U <- (runif(1) < acc_prob)
      U_cur <- U_cur * acc_U + U_old * (1 - acc_U)
      Ru_cur <- Ru_cur * acc_U + Ru_old * (1 - acc_U)
      
      # sample for V, then acc/rej
      V_update <- FALSE
      if (V_update){
        for (f_step in 1:frog_step){
          step_size_true <- 
            step_size * 
            t((1/(2 * (apply(A_d2, 2, max) + 2 * lambda))) *
                matrix(1, k, p))
          
          grad <- 
            grad_hmc_V(T, U_cur, V_cur, lambda, d1)
          
          V_cur <- 
            V_cur - 
            0.5 * step_size_true^2 * grad + 
            step_size_true * Rv_cur
          Rv_cur <- 
            Rv_cur - 0.5 * step_size_true * grad
        }
      }
      
      # metroplis step for V
      # acc_V <- 1
      acc_prob <-
        (negloglik(T, U_cur %*% t(V_old)) - 
           negloglik(T, U_cur %*% t(V_cur)) +
           lambda *
           (sum(V_old^2) - sum(V_cur^2)) +
           0.5 *
           (sum(Rv_old^2) - sum(Rv_cur^2))
        ) %>% min(0, .) %>% exp
      
      acc_V <- (runif(1) < acc_prob)
      V_cur <- V_cur * acc_V + V_old * (1 - acc_V)
      Rv_cur <- Rv_cur * acc_V + Rv_old * (1 - acc_V)
      
      # record
      if (iter %% record_freq == 0){
        rec$U[iter/record_freq, , ] <- U_cur
        rec$V[iter/record_freq, , ] <- V_cur
        rec$Theta[iter/record_freq, , ] <- U_cur %*% t(V_cur)
        
        rec$acc[iter/record_freq, ] <-
          c(mean(acc_U), mean(acc_V))
        rec$error[iter/record_freq] <-
          mean((U_cur %*% t(V_cur) - true_theta)^2) %>% sqrt
        rec$obj[iter/record_freq] <- 
          negloglik(T, U_cur %*% t(V_cur)) + 
          lambda * (sum(V_cur^2) + sum(U_cur^2))
        
        rec$time[iter/record_freq] <-
          (proc.time()[3] - time0)/60
        
        # stop if time used up
        if ((proc.time()[3] - time0)/60 > time_max){
          break
        }
      }
      
    }
    
    # return
    rec
  }