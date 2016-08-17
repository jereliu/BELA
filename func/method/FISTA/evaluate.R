objective_fista <- 
  function(par, prior, info, type = "original"){
    link_func <- 
      paste0("link_", prior$Q$link) %>% 
      parse(text = .) %>% eval
    f_Q <- link_func(par$Q)
    sigma_T_Q <- (par$T %*% t(par$sigma)) * f_Q
    sigma_T_Q[sigma_T_Q == 0] <- .Machine$double.eps
    
    loss <- 
      sigma_T_Q - info$stat$n_ij * log(sigma_T_Q)
    
    penalty <- 
      prior$Q$lambda * 
      c(norm(par$X, "F")^2, 
        norm(par$Y, "F")^2)
    
    
    
    list(total = sum(loss) + sum(penalty), 
         scaler = c(sum(loss), penalty), 
         loss_mat = loss)
  }

