default <- FALSE

if (default){
  samplr_name = c("gibbs", "hmc", "hmc_nuts", "rotation")[1]
  family_name = c("gaussian", "poisson", "binomial")[1]
  iter_max = 1e4
  record_freq = 10
  time_max = 60 
}

glrm <- 
  function(
    Y, X_r = NULL, X_c = NULL, 
    k = NULL, true_theta = NULL,
    samplr_name = c("gibbs", "hmc", "hmc_nuts", "rotation"),
    family_name = c("gaussian", "poisson", "binomial"),
    iter_max = 1e4, record_freq = 10,
    time_max = 60 
  )
  {
    # time_max: maximum sampling time, in minutes
    
    samplr_name <- match.arg(samplr_name)
    family_name <- match.arg(family_name)
    
    #### 1. initiate ####
    n <- nrow(Y)
    p <- ncol(Y)
    if (is.null(k)) k <- min(n, p)
    
    init <- NULL  
    init$U <- matrix(rnorm(n*k), nrow = n)
    init$V <- matrix(rnorm(n*p), nrow = p) 
    
    #### 2. Output Container ####
    rec_U <- array(NaN, c(n, k, iter_max/record_freq))
    rec_V <- array(NaN, c(n, k, iter_max/record_freq))
    
    #### 3. Sampler with time stamp ####
    sampler_func <- 
      paste0("glrm_sampler_", samplr_name) %>%
      parse(text = .) %>% eval()
      
      
    #### 4. results ####
  }