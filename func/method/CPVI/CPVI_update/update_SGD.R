update_SGD <- 
  function(par, prior, info, 
           par_name = "T",
           mode = c("marginal", "copula")[1], 
           verbose = FALSE)
  {
    # DOCUMENTATION
    # Variational SGD for T_j ~ Gamma(n^j, lambda_{T, j})
    # returns updated parameter par_name
    
    # INPUT
    # par_name: name of parameter to be updated
    # mode: decide parameter (marginal or copula) to update
    
    #### 1. calculate stepsize ####
    step_t <- step_size(info$oper$iter)
    
    #### 2. calculate gradient ####
    grad_t <- # a J * 1 vector
      gradient(par, prior, info, 
               par_name, mode = mode, 
               verbose)
    
    #### 3. update ####
    par[[mode]][[par_name]] <- 
      par[[mode]][[par_name]] + step_t * grad_t
    
    #### 4. return ####
    par
  }