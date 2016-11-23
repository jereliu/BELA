# TODO: Dynamically adjust step size
library(rstan)
library(parallel)

glrm_sampler_hmc_stan_U <- 
  function(Y, lambda, family_name, 
           init, config, rec, info){
    # unpack family properties
    n <- info$n 
    p <- info$p
    k <- info$k
    true_theta <- info$true_par$theta
    
    family <- glmr_family(family_name)
    T_suff <- family$sufficient(Y)
    negloglik <- family$negloglik
    
    # unpack mcmc parameters
    record_freq <- config$record_freq
    time_max <- config$time_max
    
    iter_max <- config$sampler$iter_max
    step_size <- config$sampler$step_size
    frog_step <- config$sampler$frog_step  
    mmtm_freq <- config$sampler$mmtm_freq 
    rotn_freq <- config$sampler$rotn_freq  
    
    U_cur <- init$U
    V_cur <- init$V
    Ru_cur <- matrix(rnorm(n*k), ncol = k)
    Rv_cur <- matrix(rnorm(p*k), ncol = k) 
    
    # initiate sampler
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())    
    stan_addr <- "./func/sampler/stan/"
    model_name <- family_name
    stan_file <- paste0(stan_addr, model_name, "_U.stan")
    
    stan_data <- list(N = n, P = p, K = k, Y = Y, 
                      V = init$V,
                      lambda_u = lambda, 
                      lambda_v = lambda)
    init_func <- function() init
    
    # execute sampler
    cat("\n Compiling..")
    time0 <- Sys.time()
    model_out <- 
      stan(stan_file, model_name,
           chains = 1, data = stan_data,
           iter = iter_max, 
           init = init_func,
           seed = 100, 
           pars = NA, # want all parameters
           verbose = TRUE
      )
    # time0 <- as.POSIXlt(model_out@date, format = "%c")
    time_max <- difftime(Sys.time(), time0, units = "mins")
    cat("Done")
    
    
    # return
    rec_list <- model_out@sim$samples[[1]]
    rec_name <- sapply(names(rec_list), function(x) gsub("\\[.*", "", x)) 
    
    rec <- NULL

    for (name in unique(rec_name)){
      idx <- which(rec_name == name)
      if (length(idx) > 1){
        dim_2 <- 
          names(rec_list)[max(idx)] %>% gsub("^.*,|\\]", "", .) %>% as.numeric
        dim_1 <- length(idx)/dim_2
      } else {
        dim_1 = dim_2 = 1
      }
      
      rec[[name]] <- 
        rec_list[idx] %>% abind(along = 0) %>% t %>% 
        array(dim = c(iter_max, dim_1, dim_2))
    }
    
    #rec <- extract(model_out)
    # rec$Theta <- 
    #   lapply(1:dim(rec$U)[1], 
    #          function(i) rec$U[i, ,] %*% t(rec$V[i, , ])) %>% abind(along = 0)
    
    rec$obj <-
      sapply(1:iter_max, 
             function(i){
               negloglik(T_suff, rec$Theta[i, ,]) + 
                 (lambda/2) * (sum(init$V^2) + sum(rec$U[i, ,]^2))
             }
      )

    
    rec$time <- seq(0, time_max, length.out = iter_max)
    rec$hmc_param <- get_sampler_params(model_out, inc_warmup = FALSE)
    
    # return
    rec
  }