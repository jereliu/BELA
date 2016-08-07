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
alpha <- 10

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
    lscounts = 1000, # num measure per sample
    n = n, # num of samples 
    p = p, # num of categories
    m = 20, # factor dimension
    K = 2, # population groups
    a.er = 1, b.er = 0.3, #
    sigma.value = sigma.value, 
    sigma.prob = sigma.prob, 
    strength = 3, 
    link = "soft")


#### 2. Gold Standard MCMC ####
do_MCMC <- FALSE
if (do_MCMC){
  file_addr_MCMC <- "./temp_res/MCMC/"
  if (exists("data.sim")){
    save(data.sim, file = paste0(file_addr_MCMC, "data.RData"))
  } else {
    load(paste0(file_addr_MCMC, "data.RData"))
  }
  
  hyper <- 
    list(nv = 3, a.er = 1, b.er = 0.3, 
         a1 = 3, a2 = 4, m = n )
  
  for( i in 1:length(data.sim$data) ){
    time.start <- proc.time()
    print(i)
    i = 1
    start = list( sigma = sample( sigma.value, size = p, replace = T, prob = sigma.prior ),
                  T.aug = colSums( data.sim$data[[i]] ),
                  Q = matrix( 0.5, nrow = p, ncol = n ),
                  X = matrix( rnorm( p*hyper$m ), nrow = hyper$m ),
                  Y = matrix( rnorm( n*hyper$m ), nrow = hyper$m ),
                  er = 1/rgamma( 1, shape = hyper$a.er, rate = hyper$b.er ),
                  delta = c( rgamma( 1, shape = 3, rate = 1 ), 
                             rgamma( hyper$m-1, shape = 4, rate = 1 ) ),
                  phi = matrix( rgamma( n*hyper$m, shape = 3/2, rate = 3/2 ), nrow = n ) )
    main.mcmc.shrink( data.sim$data[[i]], start, hyper, sigma.value, sigma.prior, 
                      paste0(file_addr_MCMC, "res/"), 
                      save.obj = c("sigma","Q", "Y","er"), step = 50000, thin=5 )
    time.mcmc <- proc.time() - time.start
    all.res = lapply(list.files(paste0(file_addr_MCMC, "res/"),
                                pattern = "_", full.names = TRUE), readRDS )
    all.corr = lapply( all.res, function(x) cov2cor( t(x$Y)%*%x$Y + diag( rep( x$er, ncol(x$Y) ) ) ) )
    all.corr.use = all.corr[sample(1:length(all.corr),size = 1000, replace = F)]
    
    statis.res = statis( all.corr.use, 2 )
    
    pic = plot_statis( statis.res$coord, n, 2, label = T )
    ggsave( paste(file_addr_MCMC, "plot_", i, ".pdf", sep = "" ), pic )
  }
  
  # extract MCMC prediction
  all.res <-
    lapply(list.files(
      paste0(file_addr_MCMC, "res/"),
      pattern = "_", full.names = TRUE), 
      readRDS )
  
  Q_iter <- 
    lapply(all.res, function(obj) obj$Q) %>% 
    abind(along = 3)
  mean_MCMC_Q <- apply(Q_iter, c(1, 2), mean)
  
  pred_iter <- 
    lapply(all.res, 
           function(obj){
             Q2 <- pmax(obj$Q, 0)^2
             t(obj$sigma * Q2) %>% 
               apply(1, function(x) x/sum(x)) %>% t
           }
    ) %>% abind(along = 3)
  
  mean_MCMC_pred <- apply(pred_iter, c(1, 2), mean)
  orig_mat <- 
    apply(t(data.sim$data[[1]]), 1, function(x) x/sum(x)) %>% t
  orig_mat_norm <- 
    orig_mat %>% abs %>% mean
  
  (orig_mat-mean_MCMC_pred) %>% 
    divide_by(orig_mat_norm) %>% abs %>% mean
}

#### 3. ADMM Optimization ####
file_addr_ADMM <- "./temp_res/ADMM/"

N <- data.sim$data[[1]] %>% t

prior <- 
  list(sigma = list(a = alpha/p, b = 1/2-alpha/p))
prior$lambda <- list(Q = 1, sigma = 0.1)
# prior$Q$Sigma <- 
#   list(X = diag(ncol(N)), 
#        Y = diag(nrow(N)))
         #t(data.sim$Y.tru) %*% data.sim$Y.tru)

init <- NULL
#init$sigma <- data.sim$sigma
res_ADMM <- 
  main_ADMM(
    N,
    prior = prior,
    init = init, 
    iter_max = 1e3,
    iter_crit = 1e-5,
    verbose = TRUE,
    # debug parameters
    method = "ADMM",
    par_to_update = c("Q")
  )

res_ADMM$iter <- NULL

save(res_ADMM, file = paste0(file_addr_ADMM, "res_ADMM.RData"))
load(paste0(file_addr_ADMM, "res_ADMM.RData"))

par <- res_ADMM$par
info <- res_ADMM$info
res_ADMM$iter <- NULL

#### 3.1. Reconstruction result ####
heatmap_2 <- 
  function(obj, ...){
    heatmap.2(obj, 
              Rowv = FALSE, Colv = FALSE,
              dendrogram = "none", trace = "none",
              ...)
  }
# covariance estimate
cor_tru <- 
  (t(data.sim$Y.tru) %*% data.sim$Y.tru + 
     diag(rep(data.sim$er, ncol(data.sim$Y.tru)))) %>%
  cov2cor

cor_mvi <- data.sim$Q %>% cor
cor_est <- par$Q %>% t %>% cor

heatmap_2(cor_tru)
heatmap_2(cor_mvi)
heatmap_2(cor_est)

heatmap(abs(cor_mvi - cor_est)/abs(cor_mvi))



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
