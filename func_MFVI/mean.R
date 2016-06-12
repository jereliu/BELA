library(pracma) # numeric integration 
library(symmoments)

# T: Gamma dist
# sigma: 'exact' exponential family
# Q: (n=0) Gamma
#    (n>0) Mixture Normal: (same update as 'exact')

mean_T <- function(lambda, info, prior, verbose){
  if (verbose) 
    cat("..E(T)..")
  
  # simple Gamma mean
  lambda_T <- lambda$T
  n_j <- info$stat$n_j
  
  #
  n_j/lambda_T
}

mean_sigma <- function(lambda, info, prior, 
                       eps = 1e-10, verbose){
  # brute force numeric mean for log pdf kernel
  # lambda_sigma * sig + a * log sig + b * log (1-sig)
  if (verbose) 
    cat("..E(sigma)..")
  
  a <- prior$sigma$a
  b <- prior$sigma$b
  n_i <- info$stat$n_i
  I <- info$stat$I
  lambda_sigma <- lambda$sigma
  
  #
  nom <- sapply(
    1:I,
    function(i)
      integral(fun = sigma_kernel, 
               xmin = eps, xmax = 1 - eps, 
               lambda = lambda_sigma[i], 
               alpha = a + n_i[i] + 1, beta = b)
  ) %>% log
  
  denom <- sapply(
    1:I,
    function(i)
      integral(fun = sigma_kernel, 
               xmin = eps, xmax = 1 - eps,
               lambda = lambda_sigma[i], 
               alpha = a + n_i[i], beta = b)
  ) %>% log
  
  # return expectation
  res <- exp(nom - denom)
  
  # numeric sanity check
  # sig <- seq(0, 1, 1e-3)
  # plot(sig, 
  #      sigma_kernel(sig, lambda_sigma[i], a + n_i[i], b), 
  #      type = "l")
  
  res
}

#### 2. posterior kernel function ####

sigma_kernel <- 
  function(sigma, lambda, alpha, beta){
    exp(-lambda*sigma) * 
      sigma^(alpha-1) * (1-sigma)^(beta-1)
    # # numeric sanity check
    # sig <- seq(0, 1, 1e-3)
    # plot(sig, sigma_kernel(sig, 500, a + n_i[i] + 1, b), type = "l")
  }

mean_Q1_wt_n <- 
  function(lambda_G1, lambda_G2){
    # function to calculate E(|Q|_+) when n_ij > 0
    # appr with gamma(lambda_Q1, lambda_Q2)
    lambda_G1/lambda_G2
  }

mean_Q2_wt_n <- 
  function(lambda_G1, lambda_G2){
    # function to calculate E(|Q|_+^2) when n_ij > 0
    # appr with gamma(lambda_Q1, lambda_Q2)
    (lambda_G1 + lambda_G1^2)/(lambda_G2^2)
  }

mean_Q1_no_n <- 
  function(lambda_Q1, lambda_Q2, s_j){
    # function to calculate E(|Q|_+) when n_ij = 0
    # appr with mixture normal    
    mu_1 <- lambda_Q1
    mu_2 <- lambda_Q1/(1 + 2 * lambda_Q2 * s_j)
    sd_1 <- sqrt(s_j)
    sd_2 <- sd_1/sqrt(1 + 2 * lambda_Q2 * s_j)
    
    #
    m_1 <- 
      mu_1 - sd_1 * 
      dnorm(-mu_1/sd_1)/pnorm(-mu_1/sd_1)
    
    m_2 <- 
      mu_2 + sd_2 * 
      dnorm(-mu_2/sd_2)/(1-pnorm(-mu_2/sd_2))
    
    #   
    p_denom <- 
      dnorm(0, mu_1, sd_1) * (1-pnorm(0, mu_2, sd_2)) + 
      dnorm(0, mu_2, sd_2) * pnorm(0, mu_1, sd_1)
    
    p_1 <- 
      dnorm(0, mu_2, sd_2) * pnorm(0, mu_1, sd_1)/
      p_denom
    p_2 <- 
      dnorm(0, mu_1, sd_1) * (1-pnorm(0, mu_2, sd_2))/
      p_denom
    
    #
    (p_1 * m_1 + p_2 * m_2)/(p_1 + p_2)
  }

mean_Q2_no_n <- 
  function(lambda_Q1, lambda_Q2, s_j){
    # function to calculate E(|Q|_+^2) when n_ij = 0
    # appr with mixture normal
    
    mu_1 <- lambda_Q1
    mu_2 <- lambda_Q1/(1 + 2 * lambda_Q2 * s_j)
    sd_1 <- sqrt(s_j)
    sd_2 <- sqrt(s_j)/sqrt(1 + 2 * lambda_Q2 * s_j)
    
    # calculate truncated mean and var
    alpha <- -mu_2/sd_2
    lambda <- dnorm(-alpha)/(1-pnorm(-alpha))
    delta <- lambda * (lambda - alpha)
    
    m_2 <- mu_2 + sd_2 * lambda
    var_2 <- sd_2^2 * (1 - delta)
    
    # calculate mixing proportion
    p_denom <- 
      dnorm(0, mu_1, sd_1) * (1-pnorm(0, mu_2, sd_2)) + 
      dnorm(0, mu_2, sd_2) * pnorm(0, mu_1, sd_1)
    p_2 <- 
      dnorm(0, mu_1, sd_1) * (1-pnorm(0, mu_2, sd_2))/
      p_denom
    
    #   
    p_2 * (var_2 + m_2^2)
  }


test_mean_Q1 <- FALSE
test_mean_Q2 <- FALSE

if (test_mean_Q1){
  mean_Q1_no_n(lambda_Q1 = lambda_Q1[j, i],
               lambda_Q2 = lambda_Q2[j, i],
               s_j = s_j[j])
  
  integral(fun = Q_kernelxQ,
           xmin = -10,
           xmax = 10,
           lambda_Q1 = lambda_Q1[j, i],
           lambda_Q2 = lambda_Q2[j, i],
           n_ij = n_ij[j, i],
           s_j = s_j[j],
           const_adj = mode_est[j, i]^(2*n_ij[j, i]))/
    integral(fun = Q_kernel,
             xmin = -10,
             xmax = 10,
             lambda_Q1 = lambda_Q1[j, i],
             lambda_Q2 = lambda_Q2[j, i],
             n_ij = n_ij[j, i],
             s_j = s_j[j],
             const_adj = mode_est[j, i]^(2*n_ij[j, i]))
}

if (test_mean_Q2) {
  mean_Q2_no_n(lambda_Q1 = lambda_Q1[j, i],
               lambda_Q2 = lambda_Q2[j, i],
               s_j = s_j[j])
  
  # numeric ground truth
  integral(fun = Q_kernel, 
           xmin = mode_est[j, i] - 10*sd[j], 
           xmax = mode_est[j, i] + 10*sd[j], 
           lambda_Q1 = lambda_Q1[j, i], 
           lambda_Q2 = lambda_Q2[j, i],
           n_ij = n_ij[j, i] + 1, 
           s_j = s_j[j],
           const_adj = normalizer[j, i])/
    integral(fun = Q_kernel, 
             xmin = mode_est[j, i] - 10*sd[j], 
             xmax = mode_est[j, i] + 10*sd[j], 
             lambda_Q1 = lambda_Q1[j, i], 
             lambda_Q2 = lambda_Q2[j, i],
             n_ij = n_ij[j, i], 
             s_j = s_j[j], 
             const_adj = normalizer[j, i])
}