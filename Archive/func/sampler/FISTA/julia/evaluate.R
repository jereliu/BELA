objective_fista <- 
  function(par, prior, info, type = "original"){
    link_func <- 
      paste0("link_", prior$Q$link) %>% 
      parse(text = .) %>% eval
    
    par$Q <- par$X %*% t(par$Y)
    f_Q <- link_func(par$Q)
    sigma_T_Q <- (par$T %*% t(par$sigma)) * f_Q
    sigma_T_Q[sigma_T_Q == 0] <- .Machine$double.eps
    
    loss <- 
      sigma_T_Q - info$stat$n_ij * log(sigma_T_Q)
    
    penalty <- 
      c(prior$lambda$X * norm(par$X, "F")^2, 
        prior$lambda$Y * norm(par$Y, "F")^2)
    
    list(total = sum(loss) + sum(penalty), 
         scaler = c(sum(loss), penalty), 
         loss_mat = loss)
  }


obj_ratio_fista <- function(par_cur, par_prev, prior, info){
  exp(objective_fista(par_cur, prior, info)$total - 
        objective_fista(par_prev, prior, info)$total)
}

# derivative of likelihood function
lik_dir <- 
  function(par, prior, info, type = "X"){
    #### 1. Prepare elements ####
    # latent factor related item
    other_fac <- 
      ifelse(type == "X", "Y", "X") %>% 
      paste0("par$", .) %>%
      parse(text = .) %>% eval
    
    # link_func-related item 
    link_func <- 
      paste0("link_", prior$Q$link) %>% 
      parse(text = .) %>% eval
    link_dir <- 
      paste0("link_", prior$Q$link, "_dir") %>% 
      parse(text = .) %>% eval
    
    Q_XY <- par$X %*% t(par$Y)
    
    f_Q <- link_func(par$Q)
    f_Q[f_Q == 0] <- 42 #numeric stability for 0 case
    
    f_Q_dir <- link_dir(par$Q)
    
    # other variables
    sigma_T <- (par$T %*% t(par$sigma))
    
    #### 2. Compute and collapse matrix ####
    comp1 <- 
      (sigma_T - info$stat$n_ij/f_Q) * f_Q_dir
    if (type == "Y") comp1 <- t(comp1)
    
    dir <- comp1 %*% other_fac
    
    #### 3. Return ####
    dir
  }
