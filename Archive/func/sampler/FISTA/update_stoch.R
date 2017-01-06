# list of stochastic gradient functions
#### 1. method assignment ####
stochas_update <- 
  function(update_item, par, info){
    # methods: info$method_mode
    # PULA: P-Langevin, no adjustment
    # MALA: P-Langevin, MC adjustment
    # SGHD: Hamiltonian with Friction, no adjustment
    stochas_update_func <- 
      info$method_mode$stoch %>% 
      paste0("stoch_update_", .) %>% 
      parse(text = .) %>% eval
    
    out_list <- 
      stochas_update_func(update_item, par, info)
    
    out_list
  }

#### 2. specific methods ####
stoch_update_PULA <- 
  function(update_item, par, info){
    noise <-   
      matrix(rnorm(length(par[[update_item]])), 
             nrow = nrow(par[[update_item]])) %>%
      sweep(1, sqrt(info$oper$L[[update_item]]), "/")
    par[[update_item]] <- 
      par[[update_item]] + noise
    
    list(par = par, info = info)
  }

stoch_update_MALA <- 
  function(update_item, par, info){
    par_prop <- par
    # propose
    noise <-   
      matrix(rnorm(length(par[[update_item]])), 
             nrow = nrow(par[[update_item]])) %>%
      sweep(1, sqrt(info$oper$L[[update_item]]), "/")
    par_prop[[update_item]] <- 
      par[[update_item]] + noise
    # accept / reject
    p_prop <-
      min(1, obj_ratio_fista(par_prop, par, prior, info))
    info$oper$acc <- acc <- rbinom(1, 1, p_prop)
    if (acc) 
      par <- par_prop
    
    list(par = par, info = info)
  }

stoch_update_SGHD <- 
  function(update_item, par, info){
    
    list(par = par, info = info)
  }
