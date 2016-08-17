library(R.utils) # time out feature
library(Matrix)

update_ADMM <- 
  function(par, prior, info, 
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
      update_func(par, prior, info, verbose = verbose)
    
    #### 3. return ####
    update_list
  }

update_ADMM_Q <- 
  function(par, prior, info, verbose = FALSE){
    # parse function
    update_func <- 
      paste0("update_ADMM_Q_",
             ifelse(is.null(prior$Q$Sigma), 
                    "vanila", "sigma"))  %>% 
      parse(text = .) %>% eval
    
    # obtain result
    update_list <- 
      update_func(par, prior, info, verbose)
    
    update_list
  }

update_ADMM_Q_vanila <- 
  function(par, prior, info, verbose = FALSE){
    verbose = TRUE
    #### 1 assemble stats ####
    # TODO: Figure out what's going on!
    sigma_T <- (par$T %*% t(par$sigma))
    N <- info$stat$n_ij
    noisy_Z <- (par$Z - par$U/prior$rho)
    
    #### 2 update ####
    #### 2.1 get dual problem by min aug lagrangian ####
    # Q, get dual problem by minimizing aug lag
    if (verbose){
      cat("\n")
      objective(par, prior, info, type = "dual") %>%
        (function(x) c(x$scaler, x$total)) %>%
        round(3) %>% print    
    }
    par$Q <-
      update_Q(par$Q, prior$rho, sigma_T, N, noisy_Z, 
               link = prior$Q$link)
    
    if (verbose)
      objective(par, prior, info, type = "dual") %>%
      (function(x) c(x$scaler, x$total)) %>%
      round(3) %>% print
    
    # Z, 
    # soft thresholding, "Blow up" for stability
    # TODO: find reliable SVD method
    blow_factor <- 1e3
    Q_aug <- par$Q
    
    SVD_Q <- svd((Q_aug + par$U/prior$rho)*blow_factor)
    thres <- which(SVD_Q$d/blow_factor > prior$lambda$Q/prior$rho)
    eigen_thres <- 
      SVD_Q$d[thres]/blow_factor - prior$lambda$Q/prior$rho
    par$Z <- 
      (SVD_Q$u[, thres] %*% 
         diag(eigen_thres, nrow = length(eigen_thres)) %*% 
         t(SVD_Q$v[, thres]))
    
    if (verbose)
      objective(par, prior, info, type = "dual") %>%
      (function(x) c(x$scaler, x$total)) %>%
      round(3) %>% print
    
    #### 2.2 maximize dual problem ####
    # U
    par$U <- par$U + prior$rho * (Q_aug - par$Z)
    
    if (verbose)
      objective(par, prior, info, type = "dual") %>%
      (function(x) c(x$scaler, x$total)) %>%
      round(3) %>% print
    
    #### 3 return ####
    list(par = par, info = info)
  }

update_ADMM_Q_sigma <- 
  function(par, prior, info, verbose = FALSE){
    verbose = FALSE
    #### 1.1. Q update ####
    # TODO: Figure out what's going on!
    sigma_T <- (par$T %*% t(par$sigma))
    N <- info$stat$n_ij
    noisy_Q <- (par$Q + par$Us/prior$rho)
    
    if (verbose){
      cat("\n")
      objective(par, prior, info, type = "dual") %>%
        (function(x) c(x$scaler, x$total)) %>%
        round(3) %>%
        (function(x) c("Q", x)) %>% print    
    }
    
    par$Q <-
      update_Q(par$Q, prior$rho, sigma_T, N, noisy_Q, 
               link = prior$Q$link)
    
    if (verbose){
      cat("\n")
      objective(par, prior, info, type = "dual") %>%
        (function(x) c(x$scaler, x$total)) %>%
        round(3) %>%
        (function(x) c("Q", x)) %>% print   
    }
    
    #### 1.2. S and W update ####
    # least square solution
    par$S <- 
      (prior$rho * par$Q - par$Us + 
         (prior$rho * par$W + par$Uw) %*% 
         prior$Q$Gamma$Y) %*% 
      solve(prior$rho * diag(info$stat$I) + 
              prior$rho * prior$Q$Sigma$Y_inv)
    
    par$Us <- 
      par$Us + prior$rho * (par$S - par$Q)
    
    if (verbose){
      cat("\n")
      objective(par, prior, info, type = "dual") %>%
        (function(x) c(x$scaler, x$total)) %>%
        round(3) %>%
        (function(x) c("S", x)) %>% print   
    }
    
    par$W <- 
      solve(
        prior$rho * diag(info$stat$J) + 
          prior$rho * prior$Q$Sigma$X_inv
        ,
        (prior$rho * par$S %*% t(prior$Q$Gamma$Y) - 
           par$Uw + 
           t(prior$Q$Gamma$X) %*% 
           (prior$rho * par$Z + par$Uz)
        )
      ) 
    
    par$Uw <- 
      par$Uw + 
      prior$rho * 
      (par$W - par$S %*% t(prior$Q$Gamma$Y))

    if (verbose){
      cat("\n")
      objective(par, prior, info, type = "dual") %>%
        (function(x) c(x$scaler, x$total)) %>%
        round(3) %>%
        (function(x) c("W", x)) %>% print   
    }    
    
    #### 1.3. Z update ####
    Z_target <- 
      prior$Q$Gamma$X %*% par$W - 
      par$Uz/prior$rho
    par$Z <- update_Z(Z_target , prior)
    
    par$Uz <- 
      par$Uz + 
      prior$rho * 
      (par$Z - prior$Q$Gamma$X %*% par$W)
    
    if (verbose){
      cat("\n")
      objective(par, prior, info, type = "dual") %>%
        (function(x) c(x$scaler, x$total)) %>%
        round(3) %>%
        (function(x) c("Z", x)) %>% print   
    } 
    
    #### 2. return ####
    list(par = par, info = info)
  }

update_Q <- 
  function(Q, rho,
           sigma_T, N, noisy_Z,
           link = "pos")
  {
    if (link == "pos"){
      # close form when link_func(X) <- X * I(X>0)
      delta <- 
        (noisy_Z - sigma_T/rho)^2 + 4*N/rho
      Q[(Q > 0)] <-
        0.5 * 
        ((noisy_Z - sigma_T/rho) + sqrt(delta))[(Q > 0)]
      Q[Q <= 0] <-
        (noisy_Z[Q <= 0]) %>% (function(x) x*(x <= 0))
    } else if (link == "pos2"){
      # close form when link_func(X) <- X * I(X>0)
      delta <- 
        (noisy_Z - sigma_T/rho)^2 + 4*N/rho
      Q[(Q > 0)] <-
        0.5 * 
        ((noisy_Z - sigma_T/rho) + sqrt(delta))[(Q > 0)]
      Q[Q <= 0] <-
        (noisy_Z[Q <= 0]) %>% (function(x) x*(x <= 0))
    } else {
      link_func <- 
        paste0("link_", link) %>% 
        parse(text = .) %>% eval
      link_func_dir <- 
        paste0("link_", link, "_dir") %>% 
        parse(text = .) %>% eval
      # gradient descent when other link function
      # Q <- par$Q
      step <- 1e-4
      thres <- 1e-4
      grad_norm <- Inf
      eps <- .Machine$double.eps
      # obj_list <- vector("numeric", 1e5)
      # grad_list <- vector("numeric", 1e5)
      
      i <- 0
      while (grad_norm > thres){
        grad <- 
          (sigma_T - N/(link_func(Q) + eps)) * 
          link_func_dir(Q) + 
          rho * (Q - noisy_Z)
        Q <- Q - step * grad
        grad_norm <- mean(grad^2)
        
        i <- i + 1
        obj <-
          sum(sigma_T * link_func(Q) -
                N * log(link_func(Q) + eps)) +
          0.5 * rho * sum((Q - noisy_Z)^2)
        # obj_list[i] <- obj
        # grad_list[i] <- grad_norm
        
        if (i %% 1000 == 0) 
          cat(paste0(round(obj, 3), ".."))
      }
      cat("\n")
    }
    Q
  }

update_Z <- 
  function(Q, prior,
           blow_factor = 1e3){
    # soft thresholding operator for target Q
    SVD_Q <- svd(Q*blow_factor)
    thres <- which(SVD_Q$d/blow_factor > prior$lambda$Q/prior$rho)
    eigen_thres <- 
      SVD_Q$d[thres]/blow_factor - prior$lambda$Q/prior$rho
    Z <- 
      (SVD_Q$u[, thres] %*% 
         diag(eigen_thres, nrow = length(eigen_thres)) %*% 
         t(SVD_Q$v[, thres]))
    Z
  }