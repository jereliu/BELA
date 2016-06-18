step_size <- function(t, type = "Robbins-Monro")
  {
    # Robbins-Monro schedule
    # 
    # DOCUMENTATION:
    # Calculate step size \rho_t, currently supported type are:
    # Robbins-Monro: \sum_t \rho_t = \infty, 
    #                \sum_t \rho_t^2 < \infty
  
    if (type == "Robbins-Monro") rho_t <- 1/t
    
    rho_t
  }