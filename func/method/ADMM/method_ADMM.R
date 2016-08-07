require(irlba)

default_val <- FALSE
if (default_val){
  prior = prior
  init = init 
  iter_max = 1e4
  iter_crit = 1e-5
  verbose = TRUE
  # debug parameters
  method = "ADMM"
  par_to_update = c("Q")
}

main_ADMM <- 
  function(N,
           prior = NULL,
           init = NULL, 
           iter_max = 1e4, iter_crit = 1e-3, 
           method = c("ADMM", "SGLD"),
           verbose = FALSE,
           # choose parameters to update
           par_to_update = c("T", "sigma", "Q")
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
    info$stat$I <- I
    info$stat$J <- J
    info$stat$n_j <- rowSums(N)
    info$stat$n_i <- colSums(N)
    info$stat$n_ij <- N
    
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
      prior$Q$link <- "pos"
    }
    
    if (!is.null(prior$Q$Sigma)){
      prior$Q$Sigma$X <- prior$Q$Sigma$X %>% cov2cor
      prior$Q$Sigma$Y <- prior$Q$Sigma$Y %>% cov2cor
    }
    
    if (!is.null(prior$Q$Sigma)){
      prior$Q$Gamma$X <- 
        prior$Q$Sigma$X %>% pinv %>% eigen %>%
        (function(x) 
          diag(sqrt(x$values * (x$values > 0))) %*%
           t(x$vector)
        )
      prior$Q$Gamma$Y <- 
        prior$Q$Sigma$Y %>% pinv %>% eigen %>%
        (function(x) 
          diag(sqrt(x$values * (x$values > 0))) %*%
           t(x$vector)
        )
    }
    
    if (is.null(prior$lambda)){
      prior$lambda <- list(Q = 0.1, sigma = 0.1)
    }
    
    if (is.null(prior$rho)){
      prior$rho <- 1
    }
    
    #### 1. Initialization ####
    #### > 1.1 Parameter & Sample ====
    # T
    if (is.null(init$T)) 
      init$T <- info$stat$n_j/10
    # sigma
    if (is.null(init$sigma)) 
      init$sigma <- rep(1, I)
    # Q
    if (is.null(init$Q)){
      init$Q <- 
        N / (init$T %*% t(init$sigma))
      init$Q[init$Q <= 0] <- -1
    }
    # Q, other
    if (is.null(prior$Q$Sigma)){
      init$Z <- init$Q
      init$U <- matrix(0, nrow = J, ncol = I)
    } else {
      init$S <- init$Q
      init$W <- init$S %*% t(prior$Q$Gamma$Y)
      init$Z <- prior$Q$Gamma$X %*% init$W
      
      init$Us <- matrix(0, nrow = J, ncol = I)
      init$Uw <- matrix(0, nrow = J, ncol = I)
      init$Uz <- matrix(0, nrow = J, ncol = I)
    }
    #### > 1.2 iter: Marginal Parameter History ====
    iter <- NULL
    
    iter$par$T <- array(NaN, dim = c(iter_max, J))
    iter$par$sigma <- array(NaN, dim = c(iter_max, I))
    iter$par$Q <- array(NaN, dim = c(iter_max, J, I)) 
    iter$par$Z <- array(NaN, dim = c(iter_max, J, I)) 
    iter$par$U <- array(NaN, dim = c(iter_max, J, I)) 
    
    iter$crit$par_dist <- array(NaN, dim = c(iter_max)) 
    iter$crit$feas_gap <- array(NaN, dim = c(iter_max)) 
    
    iter$crit$obj_orig <- array(NaN, dim = c(iter_max))
    iter$crit$obj_prim <- array(NaN, dim = c(iter_max)) 
    iter$crit$obj_dual <- array(NaN, dim = c(iter_max)) 
    
    #### 2. Sampling ####
    # parameters in each type to be updated
    
    #pb <- txtProgressBar(1, iter_max, style = 3)
    par_cur <- init
    par_list <- par_to_update
    time_start <- proc.time()
    par_dist_total <- Inf
    
    for (i in 1:iter_max){
      if (i > 2){
        if (feas_gap <= iter_crit){
          break
        }
      }
      #### > 2.1 Update meta info ====
      par_prev <- par_cur
      info$oper$iter <- i
      
      #### > 2.2 Update Parameter ====
      if (verbose)
        cat("\n >>> Iter", i, "Updating ")
      
      for (par_name in par_list){
        if (verbose) {
          cat("[", par_name, "]")
        }
        
        update_list <- 
          update_ADMM(
            par_cur, prior, info, 
            par_name = par_name, 
            method = method,
            verbose = verbose)
        
        # update result and mean parameter
        par_cur <- update_list$par
        info <- update_list$info
        
        #### > 2.3 Save History ====
        dist_list <- 
          lapply(names(par_cur), 
                 function(name){
                   (par_cur[[name]] - par_prev[[name]]) %>% 
                     abs %>% mean(na.rm = TRUE)
                 }) %>% set_names(names(par_cur))
        
        par_dist_total <- mean(abs(unlist(dist_list)))
        
        obj_orig <- 
          objective(par_cur, prior, info, type = "original")$total
        obj_prim <- 
          objective(par_cur, prior, info, type = "primal")$total
        obj_dual <- 
          objective(par_cur, prior, info, type = "dual")$total
        
        feas_gap <- 
          max(abs(par_cur$Q - par_cur$Z))
        
        if (verbose) {
          for(name in names(par_prev)){
            cat(paste0(name, " = ", 
                       round(dist_list[[name]], 4), "; "))
          }
          cat("rank =", rankMatrix(par_cur$Q, tol = 1e-4)[1])
          cat(", norm =", sum(svd(par_cur$Q)$d))
          cat(", neg =", sum(par_cur$Q <= 0))
          cat(", obj =", obj_orig)
        }
        
        if (i > 1){
          if (feas_gap <= iter_crit){
            cat("\n\n TOTAL CONVERGENCE!!! ᕕ( ᐛ )ᕗ \n ")
            break
          } else if ((par_dist_total > iter_crit) & (i == iter_max)){
            cat("\n limit reached and no convergence... _(:3」∠)_  \n")
          }
        }
      }
      
      # store
      iter$par$T[i, ] <- par_cur$T
      iter$par$sigma[i, ] <- par_cur$sigma
      iter$par$Q[i, , ] <- par_cur$Q
      iter$par$Z[i, , ] <- par_cur$Z
      iter$par$U[i, , ] <- par_cur$U
      
      iter$crit$par_dist[i] <- par_dist_total
      iter$crit$feas_gap[i] <- feas_gap
      iter$crit$obj_orig[i] <- obj_orig
      iter$crit$obj_prim[i] <- obj_prim
      iter$crit$obj_dual[i] <- obj_dual
    }
    
    time <- (proc.time() - time_start)/60
    
    #### 3. Return ####
    list(par = par_cur, info = info, 
         iter = iter, time = time)
  }
