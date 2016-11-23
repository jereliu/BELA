# TODO: Dynamically adjust step size
library(rstan)
library(parallel)

glrm_sampler_hmc_stan_debug <- 
  function(Y, lambda, family_name, 
           init, config, rec, info){
    # unpack family properties
    n <- info$n 
    p <- info$p
    k <- info$k
    true_theta <- info$true_par$theta
    
    family <- glmr_family(family_name)
    dist_gen <- # obtain random sampler
      substr(family_name, 1, 5) %>% 
      paste0("_gen") %>% 
      parse(text = .) %>% eval
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
    stan_file <- paste0(stan_addr, model_name, ".stan")
    model_out <- stan_model(stan_file)
    file.create("./func/debug/stan_output.txt", overwrite = TRUE)
    
    sample_iter <- 10
    time0 <- proc.time()[3]
    
    cat(paste0("\n hmc_stan (k=", k, ", family=", family_name,") initiated..."))
    sink(file = "./func/debug/stan_output.txt", 
         type = c("output", "message"))
    
    pb <- txtProgressBar(min = 1, max = iter_max, style = 3)
    
    for (iter in 1:iter_max){
      setTxtProgressBar(pb, iter)
      ############################################################
      # execute sampler for U, V
      stan_data <- list(N = n, P = p, K = k, Y = Y, 
                        lambda_u = lambda, 
                        lambda_v = lambda)
      init_func <- function() list(U = U_cur, V = V_cur)
      
      sample_cur <- 
        sampling(
          model_out,
          chains = 1, data = stan_data,
          iter = sample_iter, 
          init = init_func,
          pars = NA, # want all parameters
          refresh= -1, 
          verbose = FALSE
        ) %>% 
        suppressWarnings %>% suppressMessages
      
      rec_cur <- rstan::extract(sample_cur)
      U_cur <- rec_cur$U[sample_iter/2, , ]
      V_cur <- rec_cur$V[sample_iter/2, , ]
      Theta_cur <-  U_cur %*% t(V_cur)
      
      ############################################################
      # sample for Y
      # Loops to Sample Y using exact likelihood
      par <- family$linkinv(Theta_cur)
      Y <- par
      Y[] <- vapply(Y, dist_gen, numeric(1))
      
      T_suff <- family$sufficient(Y)
      
      # record
      if (iter %% record_freq == 0){
        rec$U[iter/record_freq, , ] <- U_cur
        rec$V[iter/record_freq, , ] <- V_cur
        rec$Theta[iter/record_freq, , ] <- U_cur %*% t(V_cur)
        
        rec$error[iter/record_freq] <- 
          mean((U_cur %*% t(V_cur) - true_theta)^2) %>% sqrt
        rec$obj[iter/record_freq] <- 
          negloglik(T_suff, U_cur %*% t(V_cur)) + 
          (lambda/2) * (sum(V_cur^2) + sum(U_cur^2))
        
        rec$time[iter/record_freq] <- 
          (proc.time()[3] - time0)/60
      }
    }
    sink()
    cat("Done. \n")
    
    
    # return
    rec
  }