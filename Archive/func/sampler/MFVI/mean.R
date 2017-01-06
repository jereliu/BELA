library(pracma) # numeric integration 
library(matrixStats)

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

mean_sigma <- function(
  lambda, info, prior, eps = 1e-15, 
  verbose, plot = FALSE){
  # brute force numeric mean for log pdf kernel
  # lambda_sigma * sig + a * log sig + b * log (1-sig)
  if (verbose) 
    cat("..E(sigma)..")
  
  a <- prior$sigma$a
  b <- prior$sigma$b
  n_i <- info$stat$n_i
  I <- info$stat$I
  lambda_sigma <- lambda$sigma
  
  # guess mode
  b_md <- (b - a - n_i)/lambda_sigma -1
  c_md <- (a + n_i - 1)/lambda_sigma
  delta_md <- (b_md^2 - 4*c_md)
  
  mode <- (0.5 * (- b_md - sqrt(delta_md))) %>% 
    pmax(0, .) %>% pmin(1)
  
  -lambda_sigma * mode*(1-mode) + 
    (a + n_i - 1) * (1 - mode) + 
    (b - 1) * mode
  
  # compute log normalizer
  log_normalizer <- 
    sapply(1:I, 
           function(i){
             sigma_val <- seq(0, 1, 1e-5)
             kernel_val <- 
               sigma_kernel(
                 sigma_val, 
                 lambda = lambda_sigma[i], 
                 alpha = a + n_i[i], beta = b,
                 log = TRUE)
             # plot(sigma_val, exp(kernel_val))
             finite_id <- is.finite(kernel_val)
             max(kernel_val[finite_id])
           }
    )
  
  x_min <- 
    sapply(1:I, 
           function(i){
             sigma_val <- seq(0, 1, 1e-5)
             kernel_val <- 
               sigma_kernel(
                 sigma_val, 
                 lambda = lambda_sigma[i], 
                 alpha = a + n_i[i], beta = b,
                 log_normalizer = log_normalizer[i])
             # plot(sigma_val, exp(kernel_val))
             finite_id <- 
               (is.finite(kernel_val) & 
                  (kernel_val > 1e-20))
             min(sigma_val[finite_id])
           }
    )
  x_max <- 
    sapply(1:I, 
           function(i){
             sigma_val <- seq(0, 1, 1e-5)
             kernel_val <- 
               sigma_kernel(
                 sigma_val, 
                 lambda = lambda_sigma[i], 
                 alpha = a + n_i[i], beta = b,
                 log_normalizer = log_normalizer[i])
             # plot(sigma_val, exp(kernel_val))
             finite_id <- 
               (is.finite(kernel_val) & 
                  (kernel_val > 1e-20))
             max(sigma_val[finite_id])
           }
    )
  
  #
  nom <- sapply(
    1:I,
    function(i)
      integral(fun = sigma_kernel,
               xmin = x_min[i], 
               xmax = x_max[i],
               lambda = lambda_sigma[i],
               alpha = a + n_i[i] + 1, beta = b,
               log_normalizer = log_normalizer[i])
  ) %>% log
  
  denom <- # something wrong with denom
    sapply(1:I,
           function(i){
             if (nom[i] == -Inf){
               1
             } else {
               integral(fun = sigma_kernel,
                        xmin = x_min[i], 
                        xmax = x_max[i],
                        lambda = lambda_sigma[i],
                        alpha = a + n_i[i], beta = b,
                        log_normalizer = log_normalizer[i])
             }
           }
    ) %>% log
  
  # return expectation
  res <- exp(nom - denom)
  
  # visual sanity check
  if (plot) {
    for (i in 1:I){
      sigma_val <- seq(0, 1, 1e-4)
      kern_val <- 
        sigma_kernel(
          sigma_val,
          lambda = lambda_sigma[i], 
          alpha = a + n_i[i], 
          beta = b, 
          log_normalizer = log_normalizer[i]
        )
      plot(sigma_val, kern_val, 
           type = "l", 
           main = 
             paste0("Iter", info$oper$iter, ": sigma ", paste0(i))
      )
      abline(v = res[i], col = 2)
    }
  }
  
  res
}

#### 2. posterior kernel function ####

sigma_kernel <- 
  function(sigma, lambda, alpha, beta, 
           log_normalizer = 0, log = FALSE){
    log_out <- 
      -lambda * sigma + 
      (alpha-1) * log(sigma) + 
      (beta-1) * log(1-sigma) - 
      log_normalizer
    
    log_out[sigma == 0] <- -Inf
    log_out[sigma == 1] <- -Inf
    
    # # numeric sanity check
    # sig <- seq(0, 1, 1e-3)
    # plot(sig, sigma_kernel(sig, 500, a + n_i[i] + 1, b), type = "l")
    
    if (log) {
      return(log_out)
    } else {
      return(exp(log_out))
    }
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
      mu_1 - 
      exp(log(sd_1) + 
            dnorm(-mu_1/sd_1, log = TRUE) - 
            pnorm(-mu_1/sd_1, log.p = TRUE)
      )
    
    m_2 <- 
      mu_2 + 
      exp(log(sd_2) + 
            dnorm(-mu_2/sd_2, log = TRUE) - 
            pnorm(-mu_2/sd_2, lower.tail = FALSE, log.p = TRUE)
      )
    
    #   
    p_denom_log <- 
      logSumExp(
        c(
          dnorm(0, mu_1, sd_1, log = TRUE) + 
            pnorm(0, mu_2, sd_2, log.p = TRUE, 
                  lower.tail = FALSE),
          dnorm(0, mu_2, sd_2, log = TRUE) + 
            pnorm(0, mu_1, sd_1, log.p = TRUE)
        )
      )
    
    p_1 <- 
      (dnorm(0, mu_2, sd_2, log = TRUE) + 
         pnorm(0, mu_1, sd_1, log.p = TRUE) -
         p_denom_log) %>% exp
    p_2 <- 
      (dnorm(0, mu_1, sd_1, log = TRUE) + 
         pnorm(0, mu_2, sd_2, lower.tail = FALSE, log.p = TRUE) -
         p_denom_log) %>% exp
    
    #
    if (p_1 == 0) {
      # handle 0*Inf = NaN
      return(m_2)
    } else if (p_2 == 0) {
      return(m_1)
    } else 
      return(p_1 * m_1 + p_2 * m_2)
  }

var_Q1_no_n <- 
  function(lambda_Q1, lambda_Q2, s_j){
    # function to calculate E(|Q|_+) when n_ij = 0
    # appr with mixture normal    
    mu_1 <- lambda_Q1
    mu_2 <- lambda_Q1/(1 + 2 * lambda_Q2 * s_j)
    sd_1 <- sqrt(s_j)
    sd_2 <- sd_1/sqrt(1 + 2 * lambda_Q2 * s_j)
    
    #
    beta <- -mu_1/sd_1
    lambda_1 <- 
      (dnorm(beta, log = TRUE) - 
         pnorm(beta, log.p = TRUE)) %>% exp
    delta_1 <- lambda_1 * (lambda_1 + beta)
    
    m_1 <- mu_1 - sd_2 * delta_1
    var_1 <- sd_2^2 * (1 - delta_1)
    
    #
    alpha <- -mu_2/sd_2
    lambda_2 <- 
      (dnorm(alpha, log = TRUE) - 
         pnorm(alpha, lower.tail = FALSE, log.p = TRUE)) %>% 
      exp
    delta_2 <- lambda_2 * (lambda_2 - alpha)
    
    m_2 <- mu_2 + sd_2 * lambda_2
    var_2 <- sd_2^2 * (1 - delta_2)    
    
    #   
    p_denom_log <- 
      logSumExp(
        c(
          dnorm(0, mu_1, sd_1, log = TRUE) + 
            pnorm(0, mu_2, sd_2, log.p = TRUE, 
                  lower.tail = FALSE),
          dnorm(0, mu_2, sd_2, log = TRUE) + 
            pnorm(0, mu_1, sd_1, log.p = TRUE)
        )
      )
    
    p_1 <- 
      (dnorm(0, mu_2, sd_2, log = TRUE) + 
         pnorm(0, mu_1, sd_1, log.p = TRUE) -
         p_denom_log) %>% exp
    p_2 <- 
      (dnorm(0, mu_1, sd_1, log = TRUE) + 
         pnorm(0, mu_2, sd_2, lower.tail = FALSE, log.p = TRUE) -
         p_denom_log) %>% exp
    
    #
    if (p_1 == 0) {
      # handle 0*Inf = NaN
      return(v_2)
    } else if (p_2 == 0) {
      return(v_1)
    } else 
      return(p_1^2 * var_1 + p_2^2 * var_2)
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
    lambda <- 
      (dnorm(alpha, log = TRUE) - 
         pnorm(alpha, lower.tail = FALSE, log.p = TRUE)) %>% 
      exp
    delta <- lambda * (lambda - alpha)
    
    m_2 <- mu_2 + sd_2 * lambda
    var_2 <- sd_2^2 * (1 - delta)
    
    # calculate mixing proportion
    p_denom_log <- 
      logSumExp(
        c(
          dnorm(0, mu_1, sd_1, log = TRUE) + 
            pnorm(0, mu_2, sd_2, log.p = TRUE, 
                  lower.tail = FALSE),
          dnorm(0, mu_2, sd_2, log = TRUE) + 
            pnorm(0, mu_1, sd_1, log.p = TRUE)
        )
      )
    p_2 <- 
      (
        dnorm(0, mu_1, sd_1, log = TRUE) + 
          pnorm(0, mu_2, sd_2, lower.tail = FALSE, log.p = TRUE) -
          p_denom_log
      ) %>% exp
    
    #   
    if (p_2 == 0) {
      # handle 0*Inf = NaN
      return(0)
    } else {
      return(p_2 * (var_2 + m_2^2))
    }
  }


test_mean_Q1 <- FALSE
test_mean_Q2 <- FALSE

if (test_mean_Q1){
  mean_Q1_no_n(lambda_Q1 = lambda$Q1[j, i],
               lambda_Q2 = lambda$Q2[j, i],
               s_j = sqrt(prior$Q$Sg_cond)[j])
  
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
              s_j = sqrt(prior$Q$Sg_cond)[j])
  
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