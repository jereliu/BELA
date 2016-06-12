library(R.utils) # time out feature

update_MFVI <- 
  function(lambda, prior, info, 
           par_name = "T", method, 
           verbose = FALSE)
  {
    # DOCUMENTATION
    # returns updated parameter par_name
    
    # INPUT
    # par_name: name of parameter to be updated
    
    ##### 1. assign update function based on parameter name ####
    update_func <- 
      # get name
      paste0("update_", method, "_", par_name) %>% 
      # find function
      parse(text = .) %>% eval
    
    #### 2. update ####
    update_list <- 
      update_func(lambda, prior, info, verbose = verbose)
    
    #### 3. return ####
    update_list
  }


update_MFVI_T <- 
  function(lambda, prior, info, verbose = verbose){
    # update mean estimates
    if (is.null(info$mean[["sigma"]])) {
      info$mean$sigma <- 
        mean_sigma(lambda, info, prior, 
                   verbose = verbose)
    }
    if (is.null(info$mean[["Q2"]])){
      info$mean$Q2 <- 
        mean_Q2(lambda, info, prior, 
                verbose = verbose)
    }
    
    # compute
    lambda$T <- 
      # compute E(sigma) * E(|Q_ij|^2)
      rep(info$mean$sigma, info$stat$J) %>% 
      matrix(nrow = info$stat$J, byrow = TRUE) %>%    
      multiply_by(info$mean$Q2) %>% 
      # add over OTU
      rowSums
    
    # return
    list(lambda = lambda, info = info)
  }

update_MFVI_sigma <- 
  function(lambda, prior, info, verbose = verbose){
    # update mean estimates
    if (is.null(info$mean[["T"]])) {
      info$mean$T <- 
        mean_T(lambda, info, prior, verbose = verbose)
    }
    if (is.null(info$mean[["Q2"]])){
      info$mean$Q2 <- 
        mean_Q2(lambda, info, prior, 
                verbose = verbose)
    }
    
    # compute
    lambda$sigma <- 
      # compute E(sigma) * E(|Q_ij|^2)
      (info$mean$T * info$mean$Q2) %>% 
      # add over sample (j's)
      colSums
    
    # return
    list(lambda = lambda, info = info)
  }

update_MFVI_Q <- 
  function(lambda, prior, info, verbose = verbose){
    # update mean estimates
    if (is.null(info$mean[["Q1"]])){
      info <- 
        mean_Q1(lambda, info, prior, 
                verbose = verbose)      
    }
    if (is.null(info$mean[["T"]])) {
      info$mean$T <- 
        mean_T(lambda, info, prior, verbose = verbose)
    }
    if (is.null(info$mean[["sigma"]])) {
      info$mean$sigma <- 
        mean_sigma(lambda, info, prior, 
                   verbose = verbose)
    }
    
    # compute first parameter (mean for Q)
    lambda$Q1 <- 
      sapply(1:info$stat$J, 
             function(j) 
               prior$Q$mu_mplr[[j]] %*% 
               info$mean$Q1[-j, ]) %>% t
    
    # compute second parameter ('var' for |Q|_+^2)
    lambda$Q2 <- 
      matrix(info$mean$T, ncol = 1) %*% 
      info$mean$sigma
    
    # return
    list(lambda = lambda, info = info)
  }

update_MFVI_G <-
  function(lambda, prior, info, 
           verbose = verbose, plot = FALSE){
    # update mean estimates
    if (is.null(info$mean[["Q1"]])){
      info <- 
        mean_Q1(lambda, info, prior, 
                verbose = verbose)
    }
    if (is.null(info$mean[["Q2"]])){
      info$mean$Q2 <- 
        mean_Q2(lambda, info, prior, 
                verbose = verbose)      
    }    
    if (is.null(info$mean[["T"]])) {
      info$mean$T <- 
        mean_T(lambda, info, prior, verbose = verbose)
    }
    if (is.null(info$mean[["sigma"]])) {
      info$mean$sigma <- 
        mean_sigma(lambda, info, prior, 
                   verbose = verbose)
    }
    
    # calculate statistics
    info$stat$cond_mu <- 
      sapply(1:info$stat$J, 
             function(j) 
               prior$Q$mu_mplr[[j]] %*% 
               info$mean$Q1[-j, ]) %>% t
    
    info$stat$TxSigma <- 
      matrix(info$mean$T, ncol = 1) %*% 
      info$mean$sigma
    
    # MoM Estimation
    lambda_G2 <- 
      (lambda$G1/lambda$G1) * 
      (info$mean$Q1/(info$mean$Q2 - info$mean$Q1^2))
    lambda_G1 <- 
      info$mean$Q1 * lambda_G2
    
    if (plot){
      n_ij <- info$stat$n_ij
      s_j <- prior$Q$Sg_cond
      sd <- sqrt(s_j)
      lambda_Q1 <- lambda$Q1
      lambda_Q2 <- lambda$Q2
      
      mode_est <- info$plot$Q1$mode_est
      normalizer <- info$plot$Q1$normalizer
      
      # plot
      for (j in 1:info$stat$J){
        for (i in 1:(info$stat$I)){
          if (!is.na(lambda_G1[j, i])){
            # plot the mofo
            Q <- seq(mode_est[j, i] - 5*sd[j], 
                     mode_est[j, i] + 5*sd[j], 1e-3)
            dens_exct <- 
              Q_kernel(
                Q,
                lambda_Q1 = lambda_Q1[j, i],
                lambda_Q2 = lambda_Q2[j, i],
                n_ij = n_ij[j, i],
                s_j = s_j[j],
                const_adj = normalizer[j, i])
            
            dens_appr <- 
              dgamma(Q, lambda_G1[j, i], lambda_G2[j, i])
            
            dens_exct <- dens_exct/max(dens_exct)
            dens_appr <- dens_appr/max(dens_appr)
            
            plot(Q, dens_exct,
                 type = "l", ylab = "", 
                 main = 
                   paste0("N_", j, "_", i, "=", info$stat$n_ij[j, i]))
            lines(Q, dens_appr, col = 2)
          }
        }
      }
    }
    
    lambda$G1 <- lambda_G1
    lambda$G2 <- lambda_G2
    
    # return
    list(lambda = lambda, info = info)
  }


