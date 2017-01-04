default <- FALSE

if(default){
  visual = FALSE
  freq_slicing = FALSE
  dir_hess = FALSE
  hit_n_run = TRUE
  rej_slices = TRUE
}

# Slice Sampler, sample by row using doubling procedure on line sampler
glrm_sampler_slice_entire_propslice <- 
  function(Y, lambda, family_name, 
           init, config, rec, info, 
           # whether visualize tree when dim = 2
           visual = FALSE, 
           freq_slicing = FALSE,
           dir_hess = FALSE,
           hit_n_run = TRUE,
           rej_slices = TRUE){
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
    
    target_lik <- function(U, V, T_suff, lambda){
      negloglik(T_suff, U %*% t(V)) + 
        (lambda/2) * (sum(V^2) + sum(U^2))
    }
    
    # unpack mcmc parameters
    iter_max <- config$sampler$iter_max
    record_freq <- config$record_freq
    time_max <- config$time_max
    rotn_freq <- config$sampler$rotn_freq
    edge_max <- config$sampler$edge_max 
    if (rej_slices) {
      line_step <- config$sampler$line_step
    } else {
      line_step <- 1
    }
    
    U_cur <- init$U
    V_cur <- init$V
    Theta_cur <- U_cur %*% t(V_cur)
    time0 <- proc.time()[3]
    
    # initiate sampler
    pb <- txtProgressBar(min = 1, max = iter_max, style = 3)
    slice_lik_U_cur <- NULL
    
    vol_ratio <- c()
    cum_vol_ratio <- 1
    
    for (iter in 1:iter_max) {
      setTxtProgressBar(pb, iter)
      # if (iter %% rotn_freq == 0) {
      #   R <- rortho(k)
      #   U_cur <- U_cur %*% R
      #   V_cur <- V_cur %*% R
      # }
      
      ####  propose slice  ################
      slice_lik_U_old <- slice_lik_U_cur
      slice_lik_U_cur <- 
        target_lik(
          U_cur, V_cur, T_suff, lambda) + rexp(1)
      
      ####  propose Theta  ################
      rec_step <- NULL
      rec_step$U <- array(NaN, dim = c(line_step, n, k))
      rec_step$V <- array(NaN, dim = c(line_step, p, k))
      
      for (step_id in 1:line_step){
        ####  Theta:U  ################
        acc_U <- rep(NaN, n)
        Theta_cur <- U_cur %*% t(V_cur)
        
        # 1. select direction
        #(TODO: optionally, 
        #       sample principal direction using 
        #       eigenvectors of hessian matrix)
        U_dir <- matrix(rmvnorm(1, rep(0, n*k)), nrow = n)
        U_dir <- U_dir/sqrt(sum(U_dir^2))
        
        # 2. sample
        U_cur <- 
          sampleLevelSet_entire(
            var = "U",
            U_cur, V_cur, T_suff, lambda,
            U_dir, V_dir = NULL, 
            slice_lik_U_cur, target_lik, 
            sample_scale = 5
          )
        
        rec_step$U[step_id, , ] <- U_cur
          
        ####  Theta:V  #################
        acc_V <- rep(NaN, p)
        
        # code to update V_cur
        V_old <- V_cur
        
        rec_step$V[step_id, , ] <- V_cur
      }
      
      ####  reject slice and corresponding Theta ################
      # rejection step after sampler
      if (rej_slices & (iter > 1)) {
        f <- function(step_id){
          target_lik(rec_step$U[step_id, , ],
                     rec_step$V[step_id, , ], 
                     T_suff, lambda) - 
            slice_lik_U_old
        }
        vol_ratio_cur <- mean(sapply(1:line_step, f) < 0)
        vol_ratio_cum <- vol_ratio_cum * vol_ratio_cur
        vol_ratio_list <- vol_ratio_cur
        
        break
        if ((vol_ratio_cum < 0.55) & (vol_ratio_cum > 0.8)){
          # WHAT to do if accept?
          
        } else {
          # WHAT to do if reject?
        }
      }
      ##########################
      # record
      if (iter %% record_freq == 0){
        rec$U[iter/record_freq, , ] <- U_cur
        rec$V[iter/record_freq, , ] <- V_cur
        rec$Theta[iter/record_freq, , ] <- U_cur %*% t(V_cur)
        
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