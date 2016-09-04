library(magrittr)
library(dplyr)
library(ggplot2)
library(abind)
library(gplots)
source("./func/util/source_Dir.R")
sourceDir("./func")

#### 1. Data Generation ####
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


n_list <- c(50, 100, 500, 1000)
p_list <- c(100, 500, 1000)
rep_num <- 1
record <- array(
  NA, dim = c(length(n_list), length(p_list), rep_num),
  dimnames = list(n = n_list, p = p_list, rep = 1:rep_num)
)


for (i in 1:length(n_list)){
  for (j in 1:length(p_list)){
    for (k in 1:rep_num){
      n <- n_list[i]
      p <- p_list[j]
      alpha <- 5
      
      print(paste0("n = ", n, ", p = ", p, ", rep = ", k))

      data.sim <- 
        boyu_sample(
          lscounts = 1e5, # num measure per sample
          n = n, # num of samples 
          p = p, # num of categories
          m = round(n/5), # factor dimension
          K = 2, # population groups
          a.er = 1, b.er = 0.3, #
          sigma.value = sigma.value, 
          sigma.weight = sigma.prob)
      
      #### 2. ADMM Optimization ####
      file_addr_ADMM <- "./temp_res/ADMM/"
      
      N <- data.sim$data %>% t
      
      prior <- 
        list(sigma = list(a = alpha/p, b = 1/2-alpha/p))
      prior$lambda <- list(Q = 1, sigma = 0.1)
      
      # prior$Q$Sigma <-
      #   list(X = t(data.sim$Y.tru) %*% data.sim$Y.tru,
      #        Y = diag(ncol(N)))
      prior$Q$link = "pos"
      
      init <- NULL
      #init$sigma <- data.sim$sigma
      res_ADMM <- 
        main_ADMM(
          N,
          prior = prior,
          init = init, 
          iter_max = 1e4,
          iter_crit = 1e-4,
          verbose = FALSE,
          # debug parameters
          method = "ADMM",
          par_to_update = c("Q")
        )
      
      #save(res_ADMM, file = paste0(file_addr_ADMM, "res_ADMM.RData"))
      #load(paste0(file_addr_ADMM, "res_ADMM.RData"))
      
      par <- res_ADMM$par
      info <- res_ADMM$info
      #res_ADMM$iter <- NULL
      
      record[i, j, k] <- res_ADMM$time[3]
    }
  }
}

record_old <- record

#### 3.1. Reconstruction result ####
# covariance estimate
cor_tru <- 
  (t(data.sim$Y.tru) %*% data.sim$Y.tru + 
     diag(rep(data.sim$er, ncol(data.sim$Y.tru)))) %>%
  cov2cor

cor_mvi <- data.sim$Q %>% cor
cor_est <- par$Q %>% t %>% cor
cor_est2 <- par$Q %>% svd %>% 
  (function(x) x$u %*% diag(x$d) %*% t(x$u))

rv_coef(cor_est, cor_mvi)
rv_coef(cor_est2, cor_mvi)

heatmap_2(cor_tru)
heatmap_2(cor_mvi)
heatmap_2(cor_est)
heatmap_2(cor_est2)

heatmap_2(abs(cor_mvi - cor_est)/abs(cor_mvi))



sigma_T <- 
  prior$Q$Sigma$Y %*% (par$T %*% t(par$sigma)) %*%
  prior$Q$Sigma$X
N <- 
  prior$Q$Sigma$Y %*% (info$stat$n_ij) %*%
  prior$Q$Sigma$X
noisy_Z <- 
  prior$Q$Sigma$Y %*% t(prior$Q$Gamma$Y) %*% 
  (par$Z - par$U/prior$rho) %*%
  prior$Q$Gamma$X %*% prior$Q$Sigma$X

par_Q <- par$Q * (par$Q > 0)
delta <- 
  (noisy_Z - sigma_T/prior$rho)^2 + 4*N/prior$rho
par_Q <-
  0.5 * ((noisy_Z - 
            sigma_T/prior$rho) + sqrt(delta))
par_Q[par_Q == 0] <-  .Machine$double.eps

sum(sigma_T * par_Q  - N * log(par_Q)) + 
  sum((par_Q - noisy_Z)^2)
