# Gibbs Sampler, rejection by row
# Note: Gibbs for Poisson case is sensitive to 
#         initialization, if not properly initialized
#         then proposals for some U[i]/V[j] that corresponds to extreme 
#         Y observation will always be rejected 

glrm_sampler_gibbs <- 
  function(Y, lambda, 
           family_name = c("gaussian", "poisson"), 
           prior_name = c("gaussian"),
           init, config, rec, info)
  {
    set.seed(config$sampler$samp_seed)
    family_name <- match.arg(family_name)
    prior_name <- match.arg(prior_name)
    
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
    parm_updt <- config$sampler$parm_updt
    
    rec$U[1, , ] <- U_old <- U_cur <- init$U
    rec$V[1, , ] <- V_old <- V_cur <- init$V
    rec$Theta[1, , ] <- Theta_cur <- U_cur %*% t(V_cur)
    time0 <- proc.time()[3]
    
    acc_U <- matrix(NaN, iter_max, n)
    acc_V <- matrix(NaN, iter_max, p)
    
    # initiate sampler
    pb <- txtProgressBar(min = 1, max = iter_max, style = 3)
    
    for (iter in 1:iter_max) {
      setTxtProgressBar(pb, iter)
      # if (iter %% rotn_freq == 0) {
      #   R <- rortho(k)
      #   U_cur <- U_cur %*% R
      #   V_cur <- V_cur %*% R
      # }
      
      ####  U  ################
      if ("U" %in% parm_updt){
        # Loops to Sample U
        U_old <- U_cur
        # generate A', A'', B
        Theta_cur <- U_cur %*% t(V_old)
        A_d1 <- d1(Theta_cur)
        A_d2 <- d2(Theta_cur)
        B <- T_suff - A_d1 + A_d2 * Theta_cur
        lnr_coef_u <- B %*% V_cur    # n x k, each row lnr coef for u_i
        
        for (i in 1:n){
          # sample for i^th row of U
          sigma_i <-
            ginv(
              t(V_cur) %*% diag(A_d2[i, ]) %*% V_cur +
                lambda * diag(k)
            )
          mu_i <- sigma_i %*% lnr_coef_u[i, ]
          
          U_prop <- U_old
          U_prop[i, ] <-
            rmvnorm(1, mean = mu_i, sigma = sigma_i)
          
          # metroplis step for U
          acc_prob <-
            acc_prob_U(U_prop, U_old, V_cur, V_old, i,
                       lambda, family, T_suff)
          # warning("no rejection for U)
          acc_U[iter, i] <- acc_prob
          
          if (runif(1) < acc_prob){
            U_cur[i, ] <- U_prop[i, ]
          }
        }
      }
      
      ####  V  #################
      if ("V" %in% parm_updt){
        V_old <- V_cur
        # generate A', A'', B
        Theta_cur <- U_cur %*% t(V_old)
        A_d1 <- d1(Theta_cur)
        A_d2 <- d2(Theta_cur)
        B <- T_suff - A_d1 + A_d2 * Theta_cur
        lnr_coef_v <- t(B) %*% U_cur # p x k, each row lnr coef for v_j
        
        # Loops to Sample V
        for (j in 1:p){
          # sample for i^th row of U
          sigma_j <- ginv(
            t(U_cur) %*% diag(A_d2[, j]) %*% U_cur +
              lambda * diag(k))
          mu_j <- sigma_j %*% lnr_coef_v[j, ]
          
          V_prop <- V_old
          V_prop[j, ] <-
            rmvnorm(1,
                    mean = mu_j,
                    sigma = sigma_j)
          
          # V_prop[j, ] <-
          #   rmvnorm(1,
          #           mean = V_cur[j, ],
          #           sigma = sigma * diag(k))
          
          # metroplis step for V
          acc_prob <-
            acc_prob_V(U_cur, U_cur, V_prop, V_old, j,
                       lambda, family, T_suff)
          
          # warning("no rejection for V")
          acc_V[iter, j] <- acc_prob
          
          if (runif(1) < acc_prob){
            V_cur[j, ] <- V_prop[j, ]
          }
        }
      }
      
      # record
      if (iter %% record_freq == 0){
        rec$U[iter/record_freq + 1, , ] <- U_cur
        rec$V[iter/record_freq + 1, , ] <- V_cur
        rec$Theta[iter/record_freq + 1, , ] <- U_cur %*% t(V_cur)
        
        rec$acc[iter/record_freq, ] <- 
          c(acc_U[iter, ], acc_V[iter, ])
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