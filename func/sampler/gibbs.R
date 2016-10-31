glrm_sampler_gibbs <- 
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
    iter_max <- config$sampler$iter_max
    record_freq <- config$record_freq
    time_max <- config$time_max
    
    U_cur <- init$U
    V_cur <- init$V
    Theta_cur <- U_cur %*% t(V_cur)
    time0 <- proc.time()[3]
    
    # initiate sampler
    pb <- txtProgressBar(min = 1, max = iter_max, style = 3)
    
    for (iter in 1:iter_max) {
      setTxtProgressBar(pb, iter)
      
      U_old <- U_cur; V_old <- V_cur
      
      # generate A', A'', B
      Theta_old <- U_old %*% t(V_old)
      A_d1 <- d1(Theta_old)
      A_d2 <- d2(Theta_old)
      B <- T + A_d1 - A_d2 * Theta_old
      
      # sample for U, then acc/rej
      lnr_coef_u <- B %*% V_old    # n x k, each row lnr coef for u_i
      U_cur <- 
        lapply(1:nrow(A_d2), 
               function(i){
                 sigma <- 
                   solve(
                     0.5 * t(V_old) %*% diag(A_d2[i, ]) %*% V_old + 
                       lambda * diag(k)
                   )
                 mu <- 0.5 * sigma %*% lnr_coef_u[i, ]
                 rmvnorm(1, mean = mu, sigma = sigma)
               }
        ) %>% do.call("rbind", .)
      
      # metroplis step for U
      acc_U <- NULL
      acc_prob <- 
        acc_prob_U(U_cur, U_old, V_cur, V_old, 
                   lambda, family)
      acc_U <- matrix(runif(n*k), n, k) < acc_prob
      U_cur <- U_cur * acc_U + U_old * (1 - acc_U)
      
      # sample for V, then acc/rej
      Theta_cur <- U_old %*% t(V_cur)
      A_d1 <- d1(Theta_cur)
      A_d2 <- d2(Theta_cur)
      B <- T + A_d1 - A_d2 * Theta_cur
      lnr_coef_v <- t(B) %*% U_cur # p x k, each row lnr coef for v_j
      
      V_update <- FALSE
      if (V_update){
        
        V_cur <-
          lapply(1:ncol(A_d2),
                 function(j){
                   sigma <- solve(
                     0.5 * t(U_cur) %*% diag(A_d2[, j]) %*% U_cur +
                       lambda * diag(k))
                   mu <- 0.5 * sigma %*% lnr_coef_v[j, ]
                   rmvnorm(1, mean = mu, sigma = sigma)
                 }
          ) %>% do.call("rbind", .)
      }
      
      # metroplis step for V
      acc_V <- NULL
      acc_prob <- 
        acc_prob_V(U_cur, U_old, V_cur, V_old, 
                   lambda, family)
      acc_V <- matrix(runif(p*k), p, k) < acc_prob
      V_cur <- V_cur * acc_V + V_old * (1 - acc_V)
      
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
      }
      
      # stop if time used up
      # if ((proc.time()[3] - time0)/60 > time_max){
      #   break
      # }
    }
    
    # return
    rec
  }