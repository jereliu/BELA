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
    if ("FISTA" %in% info$method_mode){
      par_old <- par
      info$oper$t_pre <- info$oper$t_new
      info$oper$t_new <- 
        (1 + sqrt(1 + 4 * info$oper$t_pre^2))/2
    }
    
    for (update_item in update_item_list){
      # parse function
      update_func_Q <- 
        paste0("update_FISTA_Q_", update_item)  %>% 
        parse(text = .) %>% eval
      
      # obtain ISTA result
      update_list <- 
        update_func_Q(par, prior, info, verbose)
      
      par <- update_list$par
      info <- update_list$info
      
      # "fast" adjustment
      if ("FISTA" %in% info$method_mode){
        t_step <- 
          (info$oper$t_pre-1)/info$oper$t_new
        par[[update_item]] <- 
          par[[update_item]] + t_step * 
          (par[[update_item]] - par_old[[update_item]])
      }
      
      # stochastic update
      update_list <- 
        stochas_update(update_item, par, info)
      par <- update_list$par
      info <- update_list$info
    }
    
    list(par = par, info = info)
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
    
    info$oper$L$X <- lip_F
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
    
    
    info$oper$L$Y <- lip_F
    par$Y <- 
      (D_1 %*% (D_2 %*% par$Y - nabla_F)) %>% 
      as.matrix
    
    #objective_fista(par, prior, info)$total
    list(par = par, info = info)    
  }

update_FISTA_Q_Q <- 
  function(par, prior, info, verbose = FALSE){
    par$Q <- par$Y %*% t(par$X)
    
    list(par = par, info = info)    
  }


#

stochas_update <- 
  function(update_item, par, info){
    # methods: info$method_mode
    # PULA: P-Langevin, no adjustment
    # MALA: P-Langevin, MC adjustment
    # SGHD: Hamiltonian with Friction, no adjustment
    
    if ("PULA" %in% info$method_mode){
      noise <-   
        matrix(rnorm(length(par[[update_item]])), 
               nrow = nrow(par[[update_item]])) %>%
        sweep(1, sqrt(info$oper$L[[update_item]]), "/")
      par[[update_item]] <- 
        par[[update_item]] + noise
    } else if ("MALA" %in% info$method_mode){
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
    } else if ("SGHD" %in% info$method_mode){
      # resample momentum
      
    }
    
    
    list(par = par, info = info)
  }