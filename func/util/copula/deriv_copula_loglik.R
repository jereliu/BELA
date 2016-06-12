# derivative wrt dQ(u_i) F_copula(u_1, u_2)
deriv_copula_loglik <- 
  function (par, info) 
  {
    # currently null
    J <- info$stat$J
    
    deriv <- matrix(0, nrow = J, ncol = 1)
    deriv
  }


