require(magrittr)
require(dplyr)

source("./func/util/source_Dir.R")
sourceDir("./func")

#### 1. Data Generation ####
n <- 10
p <- 20
alpha <- 0.01 * p

sigma.value <- seq(0.001,0.999,0.001)
n_sigval <- length(sigma.value)

sigma.prob <- 
  # generate quantile 
  c(0, pbeta(sigma.value, alpha/p, 1/2-alpha/p)) %>% 
  # compute interval length
  (function(tmp) (tmp[-1] - tmp[-length(tmp)])) %>%
  # force to sum to 1
  extract(-n_sigval) %>% 
  (function(tmp) c(tmp, 1-sum(tmp)))


data <- 
  boyu_sample(
    lscounts = 100, # num measure per sample
    n = n, # num of samples 
    p = p, # num of categories
    m = 3, # factor dimension
    K = 2, # population groups
    a.er = 1, b.er = 0.3, #
    sigma.value = sigma.value, 
    sigma.prob = sigma.prob, 
    strength = 3)

N <- data$data[[1]] %>% t
Sigma_true <- t(data$Y.tru) %*% data$Y.tru 
Sigma_obs <- cov(data$Q)


#### 2. Brute-Force MF Variational Inference ####
prior_sigma <- 
  list(sigma = list(a = alpha/p, b = 1/2-alpha/p))

res <- main_MFVI(
  N, Sigma_obs, 
  prior = prior_sigma,
  init = NULL, 
  iter_max = 1e3,
  iter_crit = 1e-3,
  method = "MFVI",
  verbose = TRUE)


#### 3. Reconstruction result ####
lambda <- res$lambda
info <- res$info

P_j_raw <- 
  (info$mean$sigma * t(info$mean$Q2)) %>% t
P_j <- P_j_raw/rowSums(P_j_raw)

P_emp <- N/rowSums(N)

mean(as.matrix(abs(P_j - P_emp)))

#### 2. Investigate shape of Q_kernel ####
x <- seq(-10, 10, 1e-3)
density <- 
  Q_kernel(Q_ij = x, 
           lambda_Q1 = 0.1, lambda_Q2 = 100, 
           n_ij = 0, s_j = 1)
plot(x, density, type = "l")

