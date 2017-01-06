# Slice Sampler, sample by row, use rapid-exploring random tree

glrm_sampler_stein <-
  function(Y, lambda, family_name, 
           init, config, rec, info, 
           # whether visualize tree when dim = 2
           visual_2d = FALSE){
    # unpack family properties
    n_particle <- config$sampler$n_particle
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
    step_size <- config$sampler$step_size
    auto_corr <- 0.9
    
    iter_idx <- 
      c(seq(1, round(iter_max/4), length.out = 50),
        seq(round(iter_max/4)+100, iter_max, length.out = 50)) %>%
      round %>% unique
    
    
    # initiate parameters
    S_cur <- 
      rnorm(n_particle * (n + p) * k, sd = 1) %>% 
      matrix(nrow = n_particle)
    W_cur <- matrix(step_size, nrow = n_particle, ncol = (n + p) * k)
    # adaptive parameters
    fudge_factor <- 1e-6
    G_cur <- matrix(0, nrow = n_particle, ncol = (n + p) * k)
    T_cur <- matrix(step_size, nrow = n_particle, ncol = (n + p) * k) 
    adapt_flag <- TRUE
    first_flag <- TRUE
    time0 <- proc.time()[3]
    
    # initiate sampler
    pb <- txtProgressBar(min = 1, max = iter_max, style = 3)
    rec$Theta <- array(NaN, dim = c(length(iter_idx), n, p))
    grad_hist <- -rep(1, 500)
    grad_norm <- 0
    opt_method <- "AdaGrad"
    
    for (iter in 1:iter_max) {
      # setTxtProgressBar(pb, iter)
      grad_cur <- 
        grad_ksd(S_cur, alpha = -1, 
                 T_suff, lambda, dist_family = family)$grad
      grad_adj <- grad_cur
      # (learning rate) update 
      
      if ((mean(abs(grad_cur)) < 1) & adapt_flag){
        # funky updates to accelerate convergence, not always work well
        if (opt_method == "vanilla"){
          W_cur <- matrix(step_size, nrow = n_particle, ncol = (n + p) * k)
        } else if (opt_method == "RMSprop"){
          # RMSprop
          G_cur <- auto_corr * G_cur + (1 - auto_corr) * grad_cur^2
          W_cur <- step_rate / sqrt(G_cur + 1e-6)
        } else if (opt_method == "AdaGrad"){
          # Adagrad
          G_cur <- auto_corr * G_cur + (1 - auto_corr) * grad_cur^2
          grad_adj <- grad_cur/ (sqrt(G_cur) + 1e-6)
        } else if (opt_method == "AdaDelta") {
          # # AdaDelta
          G_cur <- auto_corr * G_cur + (1 - auto_corr) * grad_cur^2
          W_cur <- sqrt(T_cur + fudge_factor)/sqrt(G_cur + fudge_factor)
        } else if (opt_method == "momentum"){
          # momentum
          G_cur <- (1 - auto_corr) * G_cur +  auto_corr * grad_cur
          grad_adj <- G_cur
        }
        # adapt_flag <- FALSE
      }
      
      # parameter update
      S_cur <- S_cur + W_cur * grad_adj 
      
      grad_diff <- mean(abs(grad_cur)) - grad_norm
      grad_hist <- c(grad_hist[-1], grad_diff)
      
      grad_norm <- grad_diff + grad_norm
      print(paste0(iter, ":", round(grad_norm, 5)))
      
      # aftermath update
      if ((mean(abs(grad_cur)) < 1) & adapt_flag){
        if (opt_method == "AdaDelta") {
          if (first_flag) {
            T_cur <- T_cur * 0
            first_flag <- FALSE
          }
          T_cur <- auto_corr * T_cur + (1 - auto_corr) * (W_cur * grad_adj)^2
        }
        if ((sum(grad_hist>0) > 20) & (grad_norm > 1e-3) & (opt_method != "vanilla")) {
          # if stuck for a while, then switch back to vanilla
          print(paste(opt_method, "appear to be stuck, switching to vanilla.."))
          orig_method <- opt_method
          opt_method <- "vanilla"
          G_cur <- matrix(0, nrow = n_particle, ncol = (n + p) * k)
          T_cur <- matrix(0, nrow = n_particle, ncol = (n + p) * k)
        }
        # if (sum(grad_hist>0)==0 & opt_method == "vanilla") {
        #   print(paste("gradient stabled, switch back to", orig_method))
        #   opt_method <- orig_method
        #   W_cur <- W_cur * 0.9 # shrink step size a bit
        # }
      }
      
      
      # record
      if (iter %in% iter_idx){
        rec_id <- which(iter_idx == iter)
        Theta_est <- 
          apply(S_cur, 1, 
                function(S){
                  S <- matrix(S, ncol = k)
                  U <- S[1:n, ]
                  V <- S[(n+1):(n+p), ]
                  U %*% t(V)
                }
          ) %>% apply(1, mean)
        
        rec$Theta[rec_id, , ] <- matrix(Theta_est, ncol = p)
        rec$obj[rec_id] <- 
          mixing_stein(grad_neglik_S,
                       S_cur, theta_row = sum(dim(T_suff)),
                       nboot = 0,
                       T_suff = T_suff,
                       lambda = lambda,
                       dist_family = family)$ksd
        rec$grad[rec_id] <- grad_norm
        rec$time[rec_id] <- (proc.time()[3] - time0)/60
        
        if (grad_norm <= 1e-3)
          break
      }
      
      
      # stop if time used up
      # if ((proc.time()[3] - time0)/60 > time_max){
      #   break
      # }
    }
    
    # return
    rec
  }
