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

