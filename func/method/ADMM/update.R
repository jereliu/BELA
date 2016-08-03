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

update_ADMM_Q_identity <- 
  function(par, prior, info, verbose = FALSE){
    verbose = FALSE
    #### 1 assemble stats ####
    # TODO: Figure out what's going on!
    sigma_T <- 
      (par$T %*% t(par$sigma))
    # prior$Q$Sigma$Y %*% (par$T %*% t(par$sigma)) %*%
    # prior$Q$Sigma$X
    N <- 
      info$stat$n_ij
    # prior$Q$Sigma$Y %*% (info$stat$n_ij) %*%
    # prior$Q$Sigma$X
    noisy_Z <- 
      (par$Z - par$U/prior$rho)
    # prior$Q$Sigma$Y %*% t(prior$Q$Gamma$Y) %*% 
    # (par$Z - par$U/prior$rho) %*%
    # prior$Q$Gamma$X %*% prior$Q$Sigma$X
    
    #### 2 update ####
    #### 2.1 get dual problem by min aug lagrangian ####
    # Q, get dual problem by minimizing aug lag
    if (verbose){
      cat("\n")
      objective(par, prior, info, type = "dual") %>%
        (function(x) c(x$scaler, x$total)) %>%
        round(3) %>% print    
    }
    
    delta <- 
      (noisy_Z - sigma_T/prior$rho)^2 + 4*N/prior$rho
    par$Q[(par$Q > 0)] <-
      0.5 * (
        (noisy_Z - 
           sigma_T/prior$rho)[(par$Q > 0)] + 
          sqrt(delta[(par$Q > 0)])
      )
    
    if (verbose)
      objective(par, prior, info, type = "dual") %>%
      (function(x) c(x$scaler, x$total)) %>%
      round(3) %>% print
    
    par$Q[par$Q <= 0] <-
      (noisy_Z[par$Q <= 0]) %>% (function(x) x*(x <= 0))
    
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

update_ADMM_Q_Sigma <- 
  function(par, prior, info, verbose = FALS){
    
  }
