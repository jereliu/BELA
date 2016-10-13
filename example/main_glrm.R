library(magrittr)
library(dplyr)
library(ggplot2)
library(abind)
library(gplots)
source("./func/util/source_Dir.R")
sourceDir("./func")

#### 1. Data Generation ####
n <- 20
p <- 100
k <- 20

data.sim <- 
  glrm_sample(
    n = n, # num of samples 
    p = p, # num of categories
    k = k, # factor dimension
    family_name = "gaussian"
  )

#### 2. FISTA Optimization ####
Y <- data.sim$Y


#save(res_ADMM, 
#     file = paste0(file_addr_ADMM, "res_ADMM.RData"))
#load(paste0(file_addr_ADMM, "res_ADMM.RData"))

par <- res_FISTA$par
info <- res_FISTA$info

obj_FISTA <- res_FISTA$iter$crit$obj
obj_ISTA <- res_FISTA$iter$crit$obj

rng <- 2e4:1e5
plot(log(-obj_ISTA[rng]), type = "l")
lines(log(-obj_FISTA[rng]), col = 2)
#res_ADMM$iter <- NULL

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

rv_coef(cor_est, cor_mvi)

heatmap_2(cor_tru)
heatmap_2(cor_mvi)
heatmap_2(cor_est)
