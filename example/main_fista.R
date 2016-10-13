library(magrittr)
library(dplyr)
library(ggplot2)
library(abind)
library(gplots)
source("./func/util/source_Dir.R")
sourceDir("./func")

#### 1. Data Generation ####
n <- 20
p <- 1000
m <- 5
K <- 3
alpha <- 2

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

data.sim <- 
  boyu_sample(
    lscounts = 1e4, # num measure per sample
    n = n, # num of samples 
    p = p, # num of categories
    m = m, # factor dimension
    K = K, # population groups
    a.er = 1, b.er = 0.3, #
    sigma.value = sigma.value, 
    sigma.weight = sigma.prob, 
    type = "mult",
    link = "exp")

#### 2. FISTA Optimization ####
file_addr_FISTA <- "./temp_res/FISTA/"

N <- data.sim$data %>% t
write.csv(N, file = "./func/method/FISTA/julia/N.csv", 
          row.names = FALSE)

cor_tru_S <- 
  (t(data.sim$Y.tru) %*% data.sim$Y.tru  + 
     diag(rep(data.sim$er, ncol(data.sim$Y.tru)))) %>%
  cov2cor

prior <- 
  list(sigma = list(a = alpha/p, b = 1/2-alpha/p))
prior$lambda <- list(X = 1, Y = 1, sigma = 0.1)

prior$Q$K <- m


prior$Q$Sigma <-
  list(X = #cor_tru_S,
        diag(nrow(N)), 
       Y = diag(ncol(N)))
prior$Q$link <- "exp"

init <- NULL
init$T <- rowSums(N)
init$sigma <- #rep(1, ncol(N))
  data.sim$sigma

init$cor_tru_Q <- data.sim$Q %>% cor

N_temp <- (init$sigma * exp(data.sim$Q)) %>% t

method_mode <- 
  list(optim = "FISTA", stoch = "MALA")
# init$sigma <- data.sim$sigma
# TODO:watch out for data scale!

res_FISTA <- 
  main_FISTA(
    N,
    prior = prior,
    init = init, 
    iter_max = 1e5,
    iter_crit = 1e-3,
    verbose = TRUE,
    # debug parameters
    method = "FISTA",
    par_to_update = c("Q"),
    record = FALSE,
    verbose_freq = 1e3,
    # 
    init_type = "random",
    method_mode = method_mode,
    # mcmc parameters
    burn_in_ratio = 0.8, 
    thin_step = 10
    )

plot(res_FISTA$iter$crit$obj[4e4:1e5], type = "l")

#save(res_ADMM, 
#     file = paste0(file_addr_ADMM, "res_ADMM.RData"))
#load(paste0(file_addr_ADMM, "res_ADMM.RData"))

par <- res_FISTA$par
info <- res_FISTA$info

par_true <- par
par_true$X <- data.sim$Y.tru %>% t
par_true$Y <- data.sim$X.tru %>% t
par_true$sigma <- data.sim$sigma
par_true$T <- rowSums(N)

# heatmap
cor_tru <- 
  (t(data.sim$Y.tru) %*% data.sim$Y.tru + 
     diag(rep(data.sim$er, ncol(data.sim$Y.tru)))) %>%
  cov2cor

cor_mvi <- data.sim$Q %>% cor
cor_est <- (par$Y %*% t(par$X)) %>% cor

cor_est <-
  res_FISTA$iter$est %>% 
  apply(1, function(x) cor(x)) %>%
  apply(1, mean) %>% 
  matrix(n, n)

cor_est <-
  res_FISTA$iter$est %>% 
  apply(1, mean) %>% 
  matrix(n, p) %>% t %>% cor

rv_coef(cor_est, cor_mvi)

heatmap_2(cor_tru)
heatmap_2(cor_mvi)
heatmap_2(cor_est)
