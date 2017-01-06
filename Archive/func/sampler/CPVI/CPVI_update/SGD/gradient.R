gradient <- 
  function(par, prior, info, 
           par_name, mode, verbose){
    # A "middleman" function that assign gradient calculation 
    # to parameter-specfic gradient functions 
    # with name "gradient_[name]"
    
    # assign function based on parameter name
    gradient_func <- 
      # get name
      paste0("gradient_", par_name) %>% 
      # find function
      parse(text = .) %>% eval
    
    # evaluate n return 
    gradient_func(par, prior, info, mode, verbose)
  }