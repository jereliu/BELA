library(R.utils) # time out feature
library(Matrix)

update_FISTA <- 
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

update_FISTA_Q <- 
  function(par, prior, info, verbose = FALSE){
    update_item_list <- c("X", "Y")
    
    for (update_item in  update_item_list){
      # parse function
      update_func_Q <- 
        paste0("update_FISTA_Q_", update_item)  %>% 
        parse(text = .) %>% eval
      
      # obtain result
      update_list <- 
        update_func_Q(par, prior, info, verbose)
      
      par <- update_list$par
      info <- update_list$info
    }
    
    update_list
  }

update_FISTA_Q_X <- 
  function(par, prior, info, verbose = FALSE){
    if (prior$Q$link != "exp") 
      stop("only exp link function supported")
    
    #objective_fista(par, prior, info)$total
    # ISTA update
    T_sigma_fQ <- 
      par$T %*% t(par$sigma) * exp(par$X %*% t(par$Y))
    nabla_F <- 
      (T_sigma_fQ - info$stat$n_ij) %*% par$Y
    lip_F <- 
      (T_sigma_fQ %*% rowSums(par$Y^2)) %>% as.numeric
    
    D_1 <- 
      diag(
        1/(lip_F + prior$Q$Sigma$X_diag)
      )
    # D_2 <- 
    # diag(lip_F +
    # prior$lambda$X * diag(prior$Q$Sigma$X_inv)) -
    # prior$lambda$X * prior$Q$Sigma$X_inv
    D_2 <- diag(lip_F) - prior$Q$Sigma$X_offd
    
    par$X <- 
      (D_1 %*% (D_2 %*% par$X - nabla_F)) %>% 
      as.matrix
    
    #objective_fista(par, prior, info)$total
    
    list(par = par, info = info)
  }

update_FISTA_Q_Y <- 
  function(par, prior, info, verbose = FALSE){
    if (prior$Q$link != "exp")
      stop("only exp link function supported")
    # ISTA update
    T_sigma_fQ <- 
      par$T %*% t(par$sigma) * exp(par$X %*% t(par$Y))
    nabla_F <- 
      t(T_sigma_fQ - info$stat$n_ij) %*% par$X
    lip_F <- 
      (t(T_sigma_fQ) %*% rowSums(par$X^2)) %>% as.numeric
    
    D_1 <- 
      diag(
        1/(lip_F + prior$Q$Sigma$Y_diag)
      )
    # D_2 <-
    #   diag(lip_F +
    #   prior$lambda$Y * diag(prior$Q$Sigma$Y_inv)) -
    #   prior$lambda$Y * prior$Q$Sigma$Y_inv
    D_2 <- diag(lip_F)-prior$Q$Sigma$Y_offd
    
    par$Y <- 
      (D_1 %*% (D_2 %*% par$Y - nabla_F)) %>% 
      as.matrix
    
    #objective_fista(par, prior, info)$total
    list(par = par, info = info)    
  }

update_FISTA_Q_Q <- 
  function(par, prior, info, verbose = FALSE){
  }


# FISTA-specific core functions
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