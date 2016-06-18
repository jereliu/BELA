VI_biom <- 
  function(N, Sigma,
           prior = NULL,
           init = NULL, 
           iter_max = 1e3, iter_crit = 1e-3, 
           n_sample = 10, 
           cop_method = c("none", "gaussian", "select")[1],
           cop_iter = 0,
           verbose = FALSE
  )
  {
    # the meta function where all things goes together
    
    # DOCUMENTATION:
    # function to perform variational inference for model 
    # outlined in Ren et al (2016) with following meta-steps
    # 1. Sample var parameters (T_j, sigma_i, Q_{ij})
    # 2. Sample copula parameteres
    
    # INPUT:
    # N:  [J x I matrix] data indicating OTU count 
    #     among I OTU categories over J samples
    # prior:      par values for parameter prior dist
    #             a length 3 list of (T, sigma, Q)
    # init:       init val for par variational parameter
    #             a length 3 list of (T, sigma, Q)
    # iter_max:   max number of iteration allowed
    # iter_crit:  convergence criteria
    # n_sample: number of samples used for SGD
    # cop_method: method to calculate copula (currently only Gaussian)
    # cop_freq: frequency in updating tree structure. if 0 then never update
    
    # INTERNAL VARIABLE:
    # info: a list of stat used for gradient calculation 
    
    #### 0. Define Global Parameters ####
    I <- nrow(N) # sample count
    J <- ncol(N) # OTU count
    
    # assemble information container
    info <- NULL
    info$stat$I <- I
    info$stat$J <- J
    info$stat$n_j <- colSums(N)
    info$stat$n_i <- rowSums(N)
    info$stat$n_ij <- N
    info$oper$n_sample <- n_sample
    
    # assemble prior container
    if (is.null(prior$sigma)){
      prior$sigma <- list(a = 0.25, b = 0.25)
    } else if (length(prior$sigma) == 1){
      prior$sigma <- list(a = prior$sigma/I, 
                          b = 0.5 - prior$sigma/I)
      if (!(prior$sigma$a > 0 & prior$sigma$a < 0.5)){
        print(paste0("sigma prior not in (0,", 0.5*I, ")"))
      }
    }
    
    if (is.null(prior$Q))
      prior$Q <- list(Sigma_inv = solve(Sigma))
    
    #### 1. Initialization ####
    #### > 1.1 Parameter & Sample ====
    # par
    if (is.null(init$T)) 
      init$T <- info$stat$n_j * 10
    if (is.null(init$sigma)) 
      init$sigma <- runif(I)
    if (is.null(init$Q)) 
      init$Q <- 
      rexp(I*J, 0.1) %>% matrix(nrow = I)
    
    # sample
    sample <- NULL
    sample$T <- 
      genfunc_T(init$T, runif(J), info$stat$n_j)
    sample$sigma <- 
      genfunc_sigma(
        init$sigma, runif(I), 
        info$stat$n_i, prior$sigma$a, prior$sigma$b)
    sample$Q <- 
      sapply(1:I, 
             function(i) 
               genfunc_Q(init$Q[i, ], u = runif(J))) %>% t
    
    info$sample <- sample
    
    #### > 1.2 iter: Marginal Parameter History ====
    iter <- 
      vector("list", length = 3) %>% 
      set_names(c("T", "sigma", "Q"))
    
    iter$T <- array(NaN, dim = c(iter_max, J))
    iter$sigma <- array(NaN, dim = c(iter_max, I))
    iter$Q <- array(NaN, dim = c(iter_max, I, J)) 
    
    iter$T[1, ] <- init$T
    iter$sigma[1, ] <- init$sigma
    iter$Q[1, , ] <- init$Q
    
    #### > 1.3 par_cur: Marginal Parameter Container ====
    par_cur <- NULL
    par_cur$marginal <- init
    par_cur$copula$Q <- 
      # all pairwise copula parameters
      rep(0, choose(J, 2))
    
    #### 2. Sampling ####
    # parameter types to be updated
    mod_list <- c("marginal")
    # c("marginal", "copula")
    
    # parameters in each type to be updated
    par_list <- list(
      marginal = c("T", "sigma", "Q"), 
      copula = "Q"
    )
    
    pb <- txtProgressBar(1, iter_max, style = 3)
    
    for (i in 1:iter_max){
      setTxtProgressBar(pb, i)
      #### > 2.1 Update meta info ====
      info$oper$iter <- i
      
      # select tree structure
      if ((cop_iter > 0) && (i %% cop_iter == 1)){
        info$copula$Q <- 
          RVineStructureSelect(
            data = info$sample$Q/max(info$sample$Q), 
            familyset = 1)
      }
      
      #### > 2.2 Update Parameter ====
      if (verbose) 
        cat("\n Iter", i, " Updating ")
      
      for (mode in mod_list){
        for (par_name in par_list[[mode]]){
          par_cur <- 
            SGD_update(
              par_cur, prior, info, 
              par_name = par_name, 
              mode = mode, 
              verbose = verbose)
        }
      }
      
      if (verbose) 
        cat("\n")
      
      #### > 2.3 Save History ====
      # compute ELBO
    }
    
    #### 3. Return ####
    iter
  }