library(mvtnorm)

default <- FALSE
if (default){
  n <- 100
  p <- 20
  k <- 10
  X_r = NULL
  X_c = NULL
  Sig_u = NULL
  Sig_v = NULL
  family_name = 
    c("gaussian", "poisson", "poisson_softplus", "binomial")[3]
  phi_sd = 1  
  noise_sd = 0.01  
}

glrm_sample <- 
  function(
    n, p, k, 
    family_name = c("gaussian", "poisson", "poisson_softplus", 
                    "poisson_reluappr", "binomial"),
    Sig_type = rep("iid", 2),
    X_r = NULL, X_c = NULL,
    Sig_u = NULL, Sig_v = NULL,
    lambda = 1, noise_sd = 0.01
  )
  {
    family_name <- match.arg(family_name)
    
    #### 1. specify fixed parameter ####
    #### > 1.1 Sig ====
    sig_u_null <- is.null(Sig_u)
    sig_v_null <- is.null(Sig_v)
    
    if (is.null(Sig_u)){
      if (Sig_type[1] == "iid")
        Sig_u <- diag(n) / (lambda)
    } 
    if (is.null(Sig_v)){
      if (Sig_type[2] == "iid")
        Sig_v <- diag(p) / (lambda)
    }
    
    #### 2. generate random parameter ####
    #### > 2.1  r_i and c_j  ####
    if (is.null(X_r)){
      r_i <- rep(0, n)
    } else {
      beta_r <- rnorm(ncol(X_r))
      r_i <- X_r %*% beta_r
    }
    
    if (is.null(X_c)){
      c_j <- rep(0, p)
    } else {
      beta_c <- rnorm(ncol(X_c))
      c_j <- X_c %*% beta_c
    }
    
    #### > 2.2  U and V ####
    U <- rmvnorm(k, sigma = Sig_u) 
    V <- rmvnorm(k, sigma = Sig_v) 
    noise <- 0 * matrix(rnorm(n*p) * noise_sd, nrow = n)
    Q <- t(U) %*% V + noise
    
    #### > 2.3  eta ####
    theta <- outer(r_i, c_j, "+") + Q
    
    #### 3. generate observation ####
    # obtain distribution family and link function
    dist_family <- glrm_family(family_name)
    dist_gen <- # obtain random sampler
      substr(family_name, 1, 5) %>% 
      paste0("_gen") %>% 
      parse(text = .) %>% eval
    
    if (class(dist_family) != "family"){
      stop("distribution family name ('", 
           family_name, "') misspecified")
    }
    
    par <- dist_family$linkinv(theta)
    
    # generate observation
    Y <- par
    Y[] <- vapply(par, dist_gen, numeric(1))
    
    outlist <- 
      list(Y = Y, theta = theta, 
           U = U, V = V, r = r_i, c = c_j, 
           linkinv = dist_family$linkinv, 
           distgen = dist_gen)
    return(outlist)
  }
