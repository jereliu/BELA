# Gibbs Sampler, rejection by row

glrm_sampler_gibbs <- 
  function(Y, lambda, family_name, 
           init, config, rec, info){
    par_update <- c("U", "V")[1]
    set.seed(config$sampler$samp_seed)
    # unpack family properties
    n <- info$n 
    p <- info$p
    k <- info$k
    true_theta <- info$true_par$theta
    
    family <- glrm_family(family_name)
    T_suff <- family$sufficient(Y)
    d1 <- family$partition$d
    d2 <- family$partition$d2
    negloglik <- family$negloglik
    
    # unpack mcmc parameters
    iter_max <- config$sampler$iter_max
    record_freq <- config$record_freq
    time_max <- config$time_max
    rotn_freq <- config$sampler$rotn_freq
    
    rec$U[1, , ] <- U_cur <- init$U
    rec$V[1, , ] <- V_cur <- init$V
    rec$Theta[1, , ] <- Theta_cur <- U_cur %*% t(V_cur)
    time0 <- proc.time()[3]
    
    # initiate sampler
    pb <- txtProgressBar(min = 1, max = iter_max, style = 3)
    
    for (iter in 1:iter_max) {
      setTxtProgressBar(pb, iter)
      # if (iter %% rotn_freq == 0) {
      #   R <- rortho(k)
      #   U_cur <- U_cur %*% R
      #   V_cur <- V_cur %*% R
      # }
      
      U_old <- U_cur; V_old <- V_cur
      
      ####  U  ################
      acc_U <- rep(NaN, n)
      if ("U" %in% par_update){
        # Loops to Sample U
        for (i in 1:n){
          U_old <- U_cur
          # generate A', A'', B
          Theta_old <- U_old %*% t(V_old)
          A_d1 <- d1(Theta_old)
          A_d2 <- d2(Theta_old)
          B <- T_suff - A_d1 + A_d2 * Theta_old
          lnr_coef_u <- B %*% V_old    # n x k, each row lnr coef for u_i
          
          # sample for i^th row of U
          U_prop <- U_cur
          sigma <-
            solve(
              t(V_cur) %*% diag(A_d2[i, ]) %*% V_cur +
                lambda * diag(k)
            )
          mu <-  sigma %*% lnr_coef_u[i, ]
          
          # U_prop[i, ] <-
          #   rmvnorm(1,
          #           mean = mu,
          #           sigma = sigma * diag(k))
          U_prop[i, ] <-
            rmvnorm(1,
                    mean = mu,
                    sigma = sigma)
          
          # metroplis step for U
          acc_prob <-
            acc_prob_U(U_prop, U_cur, V_cur, V_old, i,
                       lambda, family, T_suff)
          # warning("no rejection for U)
          acc_U[i] <- (runif(1) < acc_prob)
          
          if (acc_U[i])
            U_cur <- U_prop
        }
      }
      
      ####  V  #################
      acc_V <- rep(NaN, p)
      if ("V" %in% par_update){
        # warning("only V updated")
        # Loops to Sample V
        for (j in 1:p){
          #cat(paste0("iter ", iter, " j=", j))
          #if (j == 39) debugonce(acc_prob_V)
          
          V_old <- V_cur
          
          # generate A', A'', B
          Theta_cur <- U_cur %*% t(V_cur)
          A_d1 <- d1(Theta_cur)
          A_d2 <- d2(Theta_cur)
          B <- T_suff - A_d1 + A_d2 * Theta_cur
          lnr_coef_v <- t(B) %*% U_cur # p x k, each row lnr coef for v_j
          
          # sample for i^th row of U
          V_prop <- V_cur
          sigma <- solve(
            t(U_cur) %*% diag(A_d2[, j]) %*% U_cur +
              lambda * diag(k))
          mu <- sigma %*% lnr_coef_v[j, ]
          V_prop[j, ] <-
            rmvnorm(1,
                    mean = mu,
                    sigma = sigma)
          
          # V_prop[j, ] <-
          #   rmvnorm(1,
          #           mean = V_cur[j, ],
          #           sigma = sigma * diag(k))
          
          # metroplis step for V
          acc_prob <-
            acc_prob_V(U_cur, U_cur, V_prop, V_old, j,
                       lambda, family, T_suff)
          
          # warning("no rejection for V")
          acc_V[j] <- (runif(1) < acc_prob)
          
          if (acc_V[j])
            V_cur <- V_prop
          #cat(paste0(" acc=", acc_V[j], "\n"))
        }
      }
      
      # record
      if (iter %% record_freq == 0){
        rec$U[iter/record_freq + 1, , ] <- U_cur
        rec$V[iter/record_freq + 1, , ] <- V_cur
        rec$Theta[iter/record_freq + 1, , ] <- U_cur %*% t(V_cur)
        
        rec$acc[iter/record_freq, ] <- 
          c(mean(acc_U), mean(acc_V))
        # rec$error[iter/record_freq] <- 
        #   mean((U_cur %*% t(V_cur) - true_theta)^2) %>% sqrt
        rec$obj[iter/record_freq] <- 
          negloglik(T_suff, U_cur %*% t(V_cur)) + 
          (lambda/2) * (sum(V_cur^2) + sum(U_cur^2))
        
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