default_val <- FALSE
if (default_val){
  Sigma = Sigma_obs
  prior = prior_sigma
  init = init_T 
  iter_max = 1e4
  iter_crit = 5e-5
  method = "MFVI"
  verbose = TRUE
  # debug parameters
  par_to_update = c("sigma", "Q", "G")
  early_termination = TRUE
  early_termination_crit = 1e-4
}

main_MFVI <- 
  function(N, Sigma,
           prior = NULL,
           init = NULL, 
           iter_max = 1e3, iter_crit = 1e-3, 
           method = "MFVI", 
           verbose = FALSE,
           # choose parameters to update
           par_to_update = c("T", "sigma", "Q", "G"),
           # shut down update scheme for converged parameter
           early_termination = FALSE, 
           early_termination_crit = 1e-4
  )
  {
    # the meta function where all things goes together
    
    # DOCUMENTATION:
    # function to perform variational inference for model 
    # outlined in Ren et al (2016) with following meta-steps
    # 1. Sample var parameters (T_j, sigma_i, Q_{ij})
    # 2. Sample copula parameteres
    
    # INPUT:
    # N:  [J x I matrix] data indicating OTU count 
    #     among I OTU categories over J samples
    # prior:      par values for parameter prior dist
    #             a length 3 list of (T, sigma, Q)
    # init:       init val for par variational parameter
    #             a length 3 list of (T, sigma, Q)
    # iter_max:   max number of iteration allowed
    # iter_crit:  convergence criteria
    # cop_method: method to calculate copula (currently only Gaussian)
    # cop_freq: frequency in updating tree structure. if 0 then never update
    
    # INTERNAL VARIABLE:
    # info: a list of stat used for gradient calculation 
    
    #### 0. Define Global Parameters ####
    J <- nrow(N) # sample count
    I <- ncol(N) # OTU count
    
    # assemble information container
    info <- NULL
    info$stat$I <- I
    info$stat$J <- J
    info$stat$n_j <- rowSums(N)
    info$stat$n_i <- colSums(N)
    info$stat$n_ij <- N
    if (is.null(info$oper$Q_change_thresh))
      # threshold for change in Q values in one iteration
      info$oper$Q_change_thresh <- 10
    
    # assemble prior container
    if (is.null(prior$sigma)){
      prior$sigma <- list(a = 0.25, b = 0.25)
    } else if (length(prior$sigma) == 1){
      prior$sigma <- list(a = prior$sigma/I, 
                          b = 0.5 - prior$sigma/I)
      if (!(prior$sigma$a > 0 & prior$sigma$a < 0.5)){
        print(paste0("sigma prior not in (0,", 0.5*I, ")"))
      }
    }
    
    if (is.null(prior$Q)){
      prior$Q$Sigma_inv <- solve(Sigma)
      prior$Q$mu_mplr <- 
        lapply(1:nrow(Sigma), 
               function(i){
                 solve(Sigma[-i, -i], Sigma[i, -i])
               })
      prior$Q$Sg_cond <-
        sapply(1:nrow(Sigma), 
               function(i){
                 Sigma[i, i] - 
                   sum(Sigma[i, -i] * prior$Q$mu_mplr[[i]])
               })
    }
    
    if (any(dim(prior$Q$Sigma_inv) != c(J, J))) 
      stop("dim(Sigma)=[", dim(Sigma)[1], ", ", 
           dim(Sigma)[2], "], expect [", J, ",", J, "]")
    
    #### 1. Initialization ####
    #### > 1.1 Parameter & Sample ====
    # T
    if (is.null(init$T)) 
      init$T <- info$stat$n_j * 2
    # sigma
    if (is.null(init$sigma)) 
      init$sigma <- runif(I)
    # Q, n_ij = 0
    if (is.null(init$Q1))
      init$Q1 <-
      rexp(I*J, 10) %>% matrix(nrow = J)
    if (is.null(init$Q2))
      init$Q2 <- 
      rexp(I*J, 10) %>% matrix(nrow = J)
    # Q, n_ij > 0
    if (is.null(init$G1))
      init$G1 <-
      rep(50, I*J) %>% matrix(nrow = J)
    if (is.null(init$G2))
      init$G2 <- 
      rep(100, I*J) %>% matrix(nrow = J)
    
    #init$Q1[N!=0] <- NaN
    #init$Q2[N!=0] <- NaN
    init$G1[N==0] <- NaN
    init$G2[N==0] <- NaN
    
    #### > 1.2 iter: Marginal Parameter History ====
    iter <- NULL
    
    iter$par$T <- array(NaN, dim = c(iter_max, J))
    iter$par$sigma <- array(NaN, dim = c(iter_max, I))
    iter$par$Q1 <- array(NaN, dim = c(iter_max, J, I)) 
    iter$par$Q2 <- array(NaN, dim = c(iter_max, J, I)) 
    iter$par$G1 <- array(NaN, dim = c(iter_max, J, I)) 
    iter$par$G2 <- array(NaN, dim = c(iter_max, J, I)) 
    
    iter$E$T <- array(NaN, dim = c(iter_max, J))
    iter$E$sigma <- array(NaN, dim = c(iter_max, I))
    iter$E$Q1 <- array(NaN, dim = c(iter_max, J, I)) 
    iter$E$Q2 <- array(NaN, dim = c(iter_max, J, I)) 

    iter$crit$ELBO <- array(NaN, dim = c(iter_max)) 
    iter$crit$par_dist <- array(NaN, dim = c(iter_max)) 
    iter$crit$prd_dist <- array(NaN, dim = c(iter_max)) 
    
    #### > 1.3 par_cur: Marginal Parameter Container ====
    lambda <- init
    
    #### 2. Sampling ####
    # parameters in each type to be updated
    # G: gamma par used to appr Q when n_ij > 0
    par_list <- par_to_update
    already_conv_par <- character(0)
    
    #pb <- txtProgressBar(1, iter_max, style = 3)
    time_start <- proc.time()
    for (i in 1:iter_max){
      #setTxtProgressBar(pb, i)
      
      #### > 2.1 Update meta info ====
      lambda_prev <- lambda
      mean_prev <- info$mean
      info$mean <- NULL
      info$oper$iter <- i
      
      # need a mean_MFVI function here?
      # info$mean$sigma <- 
      #   mean_sigma(lambda, info, prior, 
      #              verbose = verbose)
      # info <- 
      #   mean_Q2(lambda, info, prior, 
      #           verbose = verbose)
      # info$mean$T <- 
      #   mean_T(lambda, info, prior, 
      #          verbose = verbose)
      # info <- 
      #   mean_Q1(lambda, info, prior, 
      #           verbose = verbose)  
      
      #### > 2.2 Update Parameter ====
      if (verbose) 
        cat("\n >>> Iter", i, "Updating ")
      
      for (par_name in par_list){
        if (verbose) {
          cat("[", par_name, "]")
        }
        
        update_list <- 
          update_MFVI(
            lambda, prior, info, 
            par_name = par_name, 
            method = method,
            verbose = verbose)
        
        # update result and mean parameter
        lambda <- update_list$lambda
        info <- update_list$info
      }
      
      # makeup temporary mean parameters
      if (verbose) cat("\n")
      if (is.null(info$mean[["sigma"]])) {
        info$mean$sigma <- 
          mean_sigma(lambda, info, prior, 
                     verbose = verbose)
      }
      if (is.null(info$mean[["Q2"]])){
        info <- 
          mean_Q2(lambda, info, prior, 
                  verbose = verbose)
      }
      if (is.null(info$mean[["T"]])) {
        info$mean$T <- 
          mean_T(lambda, info, prior, verbose = verbose)
      }
      if (is.null(info$mean[["Q1"]])){
        info <- 
          mean_Q1(lambda, info, prior, 
                  verbose = verbose)      
      }
      
      #### > 2.3 Save History ====
      if (i > 1){ # after 1st update...
        # compute distance
        dist_list <- 
          lapply(names(info$mean), 
                 function(name){
                   (info$mean[[name]] - mean_prev[[name]]) %>% 
                     abs %>% mean(na.rm = TRUE)
                 }) %>% set_names(names(info$mean))
        
        pred_dist_total <- pred_dist(info, N)
        par_dist_total <- mean(abs(unlist(dist_list)))
        elbo_total <- ELBO_full(lambda, prior, info)
        
        if (verbose) {
          for(name in names(mean_prev)){
            cat(paste0(name, " = ", 
                       round(dist_list[[name]], 4), "; "))
          } 
          cat("Predist =", pred_dist_total, " ")
          cat("ELBO =", elbo_total, " ")
          cat("Pardist =", par_dist_total, "\n")
          
          iter$crit$ELBO[i] <- elbo_total
          iter$crit$prd_dist[i] <- pred_dist_total
          iter$crit$par_dist[i] <- par_dist_total
        }
        
        # if a parameter converged, don't update it anymore
        # TODO: need to fix for Q1 and Q2
        if (early_termination) {
          conv_par_name <- 
            names(dist_list)[unlist(dist_list) <= 
                               early_termination_crit]
          new_conv_par <- 
            !(conv_par_name %in% already_conv_par)
          
          if (sum(new_conv_par) > 0){
            cat("\n\n PARAMETER", conv_par_name[new_conv_par], "CONVERGED!!! ᕕ( ᐛ )ᕗ \n ")
            for (name in conv_par_name[new_conv_par])
              par_list <- par_list[par_list != name]
            
            already_conv_par <- 
              union(already_conv_par, conv_par_name)
            
            if (length(par_list) == 0)
              cat("\n\n ALL PARAMETER CONVERGED!!! ᕕ( ᐛ )ᕗ \n ")
          }
        }
        
        if (par_dist_total <= iter_crit){
          cat("\n\n TOTAL CONVERGENCE!!! ᕕ( ᐛ )ᕗ \n ")
          break
        } else if ((par_dist_total > iter_crit) & (i == iter_max)){
          cat("\n limit reached and no convergence... _(:3」∠)_  \n")
        }
      }
      # compute ELBO
      
      # store
      iter$par$T[i, ] <- lambda$T
      iter$par$sigma[i, ] <- lambda$sigma
      iter$par$Q1[i, , ] <- lambda$Q1
      iter$par$Q2[i, , ] <- lambda$Q2
      iter$par$G1[i, , ] <- lambda$G1
      iter$par$G2[i, , ] <- lambda$G2
      
      iter$E$T[i, ] <- info$mean$T
      iter$E$sigma[i, ] <- info$mean$sigma
      iter$E$Q1[i, , ] <- info$mean$Q1
      iter$E$Q2[i, , ] <- info$mean$Q2
    }
    time <- (proc.time() - time_start)/60
    
    #### 3. Compute Linear Response Covariance ####
    
    
    
    
    #### 4. Return ####
    list(lambda = lambda, info = info, iter = iter, time = time)
  }
