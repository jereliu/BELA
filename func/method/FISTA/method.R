require(irlba)
require(Matrix)

default_val <- FALSE
if (default_val){
  prior = prior
  init = init 
  iter_max = 1e4
  iter_crit = 1e-4
  verbose = TRUE
  # debug parameters
  method = "ADMM"
  par_to_update = c("Q")
  record = FALSE
}

main_FISTA <- 
  function(N,
           prior = NULL,
           init = NULL, 
           iter_max = 1e4, iter_crit = 1e-3, 
           method = c("ADMM", "SGLD"),
           verbose = FALSE,
           # choose parameters to update
           par_to_update = c("T", "sigma", "Q"), 
           record = FALSE
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
      prior$Q$link <- "pos"
    }
    
    if (is.null(prior$Q$link)){
      prior$Q$K <- min(I, J)
    }
    
    if (!is.null(prior$Q$Sigma)){
      prior$Q$Sigma$X <- prior$Q$Sigma$X %>% cov2cor
      prior$Q$Sigma$Y <- prior$Q$Sigma$Y %>% cov2cor
      prior$Q$Sigma$X_inv <- prior$Q$Sigma$X %>% pinv
      prior$Q$Sigma$Y_inv <- prior$Q$Sigma$Y %>% pinv
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
          matrix(0, ncol = I, nrow = I - nrow(prior$Q$Gamma$X))
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
          matrix(0, ncol = J, nrow = J - nrow(prior$Q$Gamma$Y))
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
      init$T <- info$stat$n_j/10
    # sigma
    if (is.null(init$sigma)) 
      init$sigma <- rep(1, I)
    # Q
    if (is.null(init$Q)){
      init$Q <- 
        Matrix(N / (init$T %*% t(init$sigma)), 
               sparse = TRUE)
    }
    # Q, other (X and Y)
    if (is.null(prior$Q$Sigma)){
      Q_svd <- irlba(init$Q, nv = prior$Q$K)
      init$X <- Q_svd$u %*% diag(sqrt(Q_svd$d))
      init$Y <- Q_svd$v %*% diag(sqrt(Q_svd$d))
    } 
    
    #### > 1.2 iter: Marginal Parameter History ====
    if (record){
      iter <- NULL
      
      iter$par$T <- array(NaN, dim = c(iter_max, J))
      iter$par$sigma <- array(NaN, dim = c(iter_max, I))
      iter$par$Q <- array(NaN, dim = c(iter_max, J, I)) 
      
      if (is.null(prior$Q$Sigma)){
        iter$par$Z <- array(NaN, dim = c(iter_max, J, I)) 
        iter$par$U <- array(NaN, dim = c(iter_max, J, I)) 
      } else {
        iter$par$S <- array(NaN, dim = c(iter_max, J, I)) 
        iter$par$W <- array(NaN, dim = c(iter_max, J, I)) 
        iter$par$Z <- array(NaN, dim = c(iter_max, J, I)) 
        
        iter$par$Us <- array(NaN, dim = c(iter_max, J, I)) 
        iter$par$Uw <- array(NaN, dim = c(iter_max, J, I)) 
        iter$par$Uz <- array(NaN, dim = c(iter_max, J, I)) 
      }
      
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
    par_list <- par_to_update
    time_start <- proc.time()
    par_dist_total <- Inf
    
    for (i in 1:iter_max){
      if (i > 2){
        if (par_dist_total <= iter_crit){
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
        
        if (is.null(prior$Q$Sigma)) {
          feas_gap <- (par_cur$Q - par_cur$Z) %>% 
            abs %>% max
        } else {
          feas_gap <- (par_cur$Q - par_cur$S) %>% 
            abs %>% max
        }
        
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
          if (par_dist_total <= iter_crit){
            cat("\n\n TOTAL CONVERGENCE!!! ᕕ( ᐛ )ᕗ \n ")
            break
          } else if ((par_dist_total > iter_crit) & (i == iter_max)){
            cat("\n limit reached and no convergence... _(:3」∠)_  \n")
          }
        }
      }
      
      
      # store
      if (record){
        iter$par$T[i, ] <- par_cur$T
        iter$par$sigma[i, ] <- par_cur$sigma
        iter$par$Q[i, , ] <- par_cur$Q
        
        if (is.null(prior$Q$Sigma)){
          iter$par$Z[i, , ] <- par_cur$Z
          iter$par$U[i, , ] <- par_cur$U
        } else {
          iter$par$S[i, , ] <- par_cur$S
          iter$par$W[i, , ] <- par_cur$W
          iter$par$Z[i, , ] <- par_cur$Z
          
          iter$par$Us[i, , ] <- par_cur$Us
          iter$par$Uw[i, , ] <- par_cur$Uw
          iter$par$Uz[i, , ] <- par_cur$Uz
        }
        
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
      list(par = par_cur, info = info, time = time)
    
    if (record) {
      out_list$iter <- iter
    }
    
    out_list
  }
