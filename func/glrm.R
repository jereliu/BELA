require(coda)

default <- FALSE

if (default){
  lambda = 1/phi_sd
  #samplr_name = c("gibbs", "hmc", "hmc_stan", "rotation")[3]
  #family_name = c("gaussian", "poisson", "poisson_softplus", "poisson_reluaapr", "binomial")[1]
  iter_max = c(1e5, 1e3)
  record_freq = 1
  time_max = 60 
  pred_num = 100
  step_size = c(0.001, 0.1)
  frog_step = 5
  rotn_freq = 10
  mmtm_freq = iter_max[2]/10
}

glrm <- 
  function(
    Y, X_r = NULL, X_c = NULL, lambda = 1,
    k = NULL, true_par = NULL, 
    init = NULL, init_MAP = FALSE,
    samplr_name = c("gibbs", "gibbs_debug", "hmc_stan", "hmc_stan_debug", "vi_stan"),
    family_name = c("gaussian", "poisson", "poisson_softplus", "binomial"),
    # sampler parameters: generic
    iter_max = c(1e5, 1e4), 
    record_freq = 10,
    time_max = 60,
    pred_num = 100,
    ess_num = 100,
    # sampler parameters: hmc
    step_size = c(0.001, 0.001), 
    frog_step = 5,
    rotn_freq = 10,
    mmtm_freq = iter_max[2]/10
  )
  {
    # time_max: maximum sampling time, in minutes
    samplr_name <- match.arg(samplr_name)
    family_name <- match.arg(family_name)
    
    #### 1. config ####
    config <- NULL
    config$record_freq <- record_freq
    config$time_max <- time_max
    
    config$optimiz$optim_tol <- 1e-4
    config$optimiz$iter_max <- iter_max[1]
    config$optimiz$step_size <- step_size[1]
    
    config$sampler$iter_max <- iter_max[2]
    config$sampler$step_size <- step_size[2]
    config$sampler$rotn_freq <- rotn_freq
    
    if (length(grep("hmc", samplr_name)) > 0){
      config$sampler$frog_step <- frog_step
      config$sampler$mmtm_freq <- mmtm_freq
    }
    
    #### 2. initiate ####
    info <- NULL
    info$n <- nrow(Y)
    info$p <- ncol(Y)
    if (is.null(k)){
      info$k <- min(n, p)
    } else info$k <- k
    info$true_par <- true_par
    true_theta <- info$true_par$theta
    
    if (is.null(init$U)) 
      init$U <- matrix(rnorm(n*k, sd = 1e-2), nrow = n)
    if (is.null(init$V)) 
      init$V <- matrix(rnorm(p*k, sd = 1e-2), nrow = p) 
    
    #### 3. Output Container ####
    # for sampling
    rec <- NULL
    rec$U <- array(NaN, c(iter_max[2]/record_freq, n, k))
    rec$V <- array(NaN, c(iter_max[2]/record_freq, p, k))
    rec$Theta <- array(NaN, c(iter_max[2]/record_freq, n, p))
    
    rec$acc <- array(NaN, c(iter_max[2]/record_freq, 2), 
                     dimnames = list(NULL, c("U", "V")))
    
    rec$error <- array(NaN, c(iter_max[2]/record_freq))
    rec$time <- array(NaN, c(iter_max[2]/record_freq))
    rec$obj <- array(NaN, c(iter_max[2]/record_freq)) 
    
    if (init_MAP){
      # for optimization
      rec_optim <- NULL
      rec_optim$U <- array(NaN, c(iter_max[1]/record_freq, n, k))
      rec_optim$V <- array(NaN, c(iter_max[1]/record_freq, p, k))
      
      rec_optim$error <- array(NaN, c(iter_max[1]/record_freq))
      rec_optim$time <- array(NaN, c(iter_max[1]/record_freq)) 
      rec_optim$obj <- array(NaN, c(iter_max[1]/record_freq)) 
    }
    #### 4. Put initiation at MAP ####
    if (init_MAP){
      rec_optim <- 
        glrm_optimizer(Y, lambda, family_name, 
                       init, config, rec_optim, info)
      
      init <- rec_optim$output
      
      # plot_idx <- which(!is.na(rec_optim$obj))
      # plot(rec_optim$obj[plot_idx], type = "l", ylab = "likelihood")
      # plot(rec_optim$error[plot_idx], type = "l", ylab = "error")
      
      if (any(sapply(init, is.na) %>% unlist)){
        stop("NA generated when optimizing for MAP. Consider changing step_size for optimizer")
      }
    }
    
    #### 5. Sampler with time stamp ####
    sampler_func <-
      paste0("glrm_sampler_", samplr_name) %>%
      parse(text = .) %>% eval()
    
    rec <- 
      sampler_func(Y, lambda, family_name, 
                   init, config, rec, info)
    
    # calculate effective sample size and prediction
    if (length(grep("debug", samplr_name)) == 0){
      rec$pred_error <- predMeanError(rec, data.sim$theta, pred_num)
      rec$ess <- essMatrix(rec, ess_num)
      rec$esdj <- ESJD(rec, config)
    }
    
    rec$init <- init
    rec$true_theta <- true_theta
    
    #### 6. Results ####
    rec
  }