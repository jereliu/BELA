require(irlba)
require(Matrix)

default_val <- FALSE
if (default_val){
  prior = prior
  init = init 
  iter_max = 1e5
  iter_crit = 1e-3
  verbose = TRUE
  init_type = c("svd", "random")[1]
  # debug parameters
  method = "FISTA"
  par_to_update = c("Q")
  record = FALSE
  verbose_freq = 100
}

main_FISTA <- 
  function(N,
           prior = NULL,
           init = NULL, 
           iter_max = 1e4, iter_crit = 1e-3, 
           method = c("FISTA", "MALA")[1],
           verbose = FALSE,
           # choose parameters to update
           par_to_update = c("T", "sigma", "Q"), 
           record = FALSE, 
           verbose_freq = 1000,
           # method-specific 
           init_type = c("svd", "random")[1],
           method_mode = c("ISTA", "FISTA")[1],
           # mcmc parameters
           burn_in_ratio = 0.8, 
           thin_step = 100
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
    #             a length 5 list of (T, sigma, Q, lambda, rho)
    #             where: 
    #             lambda: penalty parameter, len 2 list of (sigma, Q) 
    #             rho: ADMM step size, scale default to 1      
    # init:       init val for parameters
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

    info$method_mode <- method_mode

    info$stat$I <- I
    info$stat$J <- J
    info$stat$n_j <- rowSums(N)
    info$stat$n_i <- colSums(N)
    info$stat$n_ij <- N
    
    info$stat$L$X <- rep(1, I)
    info$stat$L$Y <- rep(1, J)
    
    # assemble prior container
    if (is.null(prior$sigma)){
      prior$sigma <- list(a = 0.25, b = 0.25)
    } else if (length(prior$sigma) == 1){
      prior$sigma <- 
        list(a = prior$sigma/I, 
             b = 0.5 - prior$sigma/I)
      if (!(prior$sigma$a > 0 & prior$sigma$a < 0.5)){
        print(paste0("sigma prior not in (0,", 0.5*I, ")"))
      }
    }
    
    if (is.null(prior$Q$link)){
      prior$Q$link <- "exp"
    }
    
    if (is.null(prior$Q$link)){
      prior$Q$K <- min(I, J)
    }
    
    if (!is.null(prior$Q$Sigma)){
      prior$Q$Sigma$X <- prior$Q$Sigma$X %>% cov2cor
      prior$Q$Sigma$Y <- prior$Q$Sigma$Y %>% cov2cor
      prior$Q$Sigma$X_inv <- prior$Q$Sigma$X %>% pinv
      prior$Q$Sigma$Y_inv <- prior$Q$Sigma$Y %>% pinv
      # diagonal elements
      prior$Q$Sigma$X_diag <- 
        prior$lambda$X * diag(prior$Q$Sigma$X_inv)
      prior$Q$Sigma$Y_diag <- 
        prior$lambda$Y * diag(prior$Q$Sigma$Y_inv)
      # off diagonal elements
      prior$Q$Sigma$X_offd <- 
        prior$lambda$X * 
        (prior$Q$Sigma$X_inv - diag(diag(prior$Q$Sigma$X_inv)))
      prior$Q$Sigma$Y_offd <- 
        prior$lambda$Y *
        (prior$Q$Sigma$Y_inv - diag(diag(prior$Q$Sigma$Y_inv)))
    }
    
    if (!is.null(prior$Q$Sigma)){
      prior$Q$Gamma$X <- 
        prior$Q$Sigma$X %>% eigen %>%
        (function(x) 
          diag(1/sqrt(x$values[x$values > 1e-8])) %*%
           t(x$vector[, x$values > 1e-8])
        )
      prior$Q$Gamma$X <- 
        rbind(
          prior$Q$Gamma$X, 
          matrix(0, ncol = J, 
                 nrow = J - nrow(prior$Q$Gamma$X))
        )
      
      prior$Q$Gamma$Y <- 
        prior$Q$Sigma$Y %>% pinv %>% eigen %>%
        (function(x) 
          diag(1/sqrt(x$values[x$values > 1e-8])) %*%
           t(x$vector[, x$values > 1e-8])
        )
      prior$Q$Gamma$Y <- 
        rbind(
          prior$Q$Gamma$Y, 
          matrix(0, ncol = I, nrow = I - nrow(prior$Q$Gamma$Y))
        )
    }
    
    if (is.null(prior$lambda)){
      prior$lambda <- list(Q = 0.1, sigma = 0.1)
    }
    
    if (is.null(prior$eta)){
      # step inflation factor for fista
      prior$eta <- 1.1
    }
    
    #### 1. Initialization ####
    #### > 1.1 Parameter & Sample ====
    # T
    if (is.null(init$T)) 
      init$T <- info$stat$n_j
    # sigma
    if (is.null(init$sigma))
      init$sigma <- 
      info$stat$n_i/median(info$stat$n_i)
    # Q
    if (is.null(init$Q)){
      init$Q <- 
        ((N / (init$T %*% t(init$sigma))) + 
           .Machine$double.eps) %>% log
    }
    # Q, other (X and Y)
    if (init_type == "svd"){
      Q_svd <- irlba(init$Q, prior$Q$K)
      init$X <- Q_svd$u %*% diag(sqrt(Q_svd$d))
      init$Y <- Q_svd$v %*% diag(sqrt(Q_svd$d))
    } else if (init_type == "random"){
      init$X <- matrix(rnorm(J*prior$Q$K), nrow = J)/1e3
      init$Y <- matrix(rnorm(I*prior$Q$K), nrow = I)/1e3
    }
    
    #### > 1.2 iter: Marginal Parameter History ====
    iter <- NULL
    iter$crit$obj <- array(NaN, dim = c(iter_max))
    iter$est <- 
      array(NaN, 
            dim = 
              c(round(((1 - burn_in_ratio)*iter_max/thin_step)), p, n)
            )
    
    if (record){
      iter$par$T <- array(NaN, dim = c(iter_max, J))
      iter$par$sigma <- array(NaN, dim = c(iter_max, I))
      iter$par$Q <- array(NaN, dim = c(iter_max, J, I)) 
      
      iter$par$X <- array(NaN, dim = c(iter_max, J, K)) 
      iter$par$Y <- array(NaN, dim = c(iter_max, I, K)) 
      
      iter$crit$par_dist <- array(NaN, dim = c(iter_max)) 
      iter$crit$feas_gap <- array(NaN, dim = c(iter_max)) 
      
      iter$crit$obj_orig <- array(NaN, dim = c(iter_max))
      iter$crit$obj_prim <- array(NaN, dim = c(iter_max)) 
      iter$crit$obj_dual <- array(NaN, dim = c(iter_max)) 
    }
    
    #### 2. Sampling ####
    # parameters in each type to be updated
    
    #pb <- txtProgressBar(1, iter_max, style = 3)
    par_cur <- init
    mat_cur <- 
      (par_cur$X %*% t(par_cur$X)) %>% cov2cor
    
    par_list <- par_to_update
    time_start <- proc.time()
    
    cov_crit_cur <- Inf
    par_dist_total <- Inf
    obj_prev <- Inf
    obj_cur <- 0
    
    # MCMC record
    mcmc_rec_min <- 
      round(iter_max * burn_in_ratio)
    mcmc_rec_idx <- 1
    
    if (info$method_mode$optim == "FISTA"){
      info$oper$t_new <- 1
    }
    
    for (i in 1:iter_max){
      if (i > 2){
        if (cov_crit_cur <= iter_crit){
          break
        }
      }
      #### > 2.1 Update meta info ====
      mat_prev <- mat_cur
      par_prev <- par_cur
      obj_prev <- obj_cur
      info$oper$iter <- i
      
      #### > 2.2 Update Parameter ====
      if (verbose & (i %% verbose_freq == 0))
        cat("\n >>> Iter", i, "Updating ")
      
      for (par_name in par_list){
        if (verbose & (i %% verbose_freq == 0)) {
          cat("[", par_name, "]")
        }
        
        update_list <- 
          update_FISTA(
            par_cur, prior, info, 
            par_name = par_name, 
            method = method,
            verbose = verbose)
        
        # update result and mean parameter
        par_cur <- update_list$par
        info <- update_list$info
      }
      
      #### > 2.3 Save History ====
      obj_cur <- 
        objective_fista(par_cur, prior, info, 
                        type = "original")$total
      obj_part <- 
        objective_fista(par_cur, prior, info, 
                        type = "original")$scaler
      
      if (verbose & (i %% verbose_freq == 0)) {
        dist_list <- 
          lapply(names(par_cur), 
                 function(name){
                   (par_cur[[name]] - par_prev[[name]]) %>% 
                     abs %>% mean(na.rm = TRUE)
                 }) %>% set_names(names(par_cur))
        
        par_dist_total <- mean(abs(unlist(dist_list)))
        
        est_cur <- 
          (par_cur$Y %*% t(par_cur$X)) %>% cor
        mat_cur <- 
          (par_cur$X %*% t(par_cur$X)) %>% cov2cor
        cov_crit_cur <- norm(mat_cur - mat_prev)
        
        for(name in names(par_prev)){
          cat(paste0(name, " = ", 
                     round(dist_list[[name]], 4), "; "))
        }
        cat(", logLik =", obj_part[1])
        cat(", norm_X =", obj_part[2]/prior$lambda$X)
        cat(", norm_Y =", obj_part[3]/prior$lambda$Y)
        cat(", obj_diff =", cov_crit_cur)
        cat(", diff_est =", median(abs(est_cur - init$cor_tru_Q)))
        cat(", diff_prd =", median(abs(mat_cur - init$cor_tru_Q)))
      }
      
      if (i > 2){
        if (cov_crit_cur <= iter_crit){
          cat("\n\n TOTAL CONVERGENCE!!! ᕕ( ᐛ )ᕗ \n ")
          break
        } else if ((par_dist_total > iter_crit) & (i == iter_max)){
          cat("\n limit reached and no convergence... _(:3」∠)_  \n")
        }
      }
      
      iter$crit$obj[i] <- obj_cur
      iter$acc[i] <- info$oper$acc
      
      if (i > mcmc_rec_min){
        if (((i - mcmc_rec_min) %% thin_step) == 0){
          iter$est[mcmc_rec_idx, , ] <-
            (par_cur$Y %*% t(par_cur$X))
          mcmc_rec_idx <- mcmc_rec_idx + 1
        }
      }
      
      # store
      if (record){
        iter$par$T[i, ] <- par_cur$T
        iter$par$sigma[i, ] <- par_cur$sigma
        iter$par$Q[i, , ] <- par_cur$Q
        
        iter$crit$par_dist[i] <- par_dist_total
        iter$crit$feas_gap[i] <- feas_gap
        iter$crit$obj_orig[i] <- obj_orig
        iter$crit$obj_prim[i] <- obj_prim
        iter$crit$obj_dual[i] <- obj_dual
      }
    }
    
    time <- (proc.time() - time_start)/60
    
    #### 3. Return ####
    out_list <- 
      list(par = par_cur, 
           info = info, 
           iter = iter,
           time = time)

    out_list
  }
