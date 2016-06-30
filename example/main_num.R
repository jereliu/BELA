library(magrittr)
library(dplyr)
library(ggplot2)
library(abind)
source("./func/util/source_Dir.R")
sourceDir("./func")

#### 1. Data Generation ####
n <- 22
p <- 64
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
    lscounts = 500, # num measure per sample
    n = n, # num of samples 
    p = p, # num of categories
    m = 3, # factor dimension
    K = 2, # population groups
    a.er = 1, b.er = 0.3, #
    sigma.value = sigma.value, 
    sigma.prob = sigma.prob, 
    strength = 3)


#### 2. Gold Standard MCMC ####
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

#### 3. Brute-Force MF Variational Inference ####
file_addr_MFVI <- "./temp_res/MFVI/"

N <- data.sim$data[[1]] %>% t
Sigma_true <- t(data.sim$Y.tru) %*% data.sim$Y.tru 
Sigma_obs <- cov(data.sim$Q)

P <- 
  (data.sim$sigma * 
     pmax(data.sim$Q, 0)^2) %>% t %>% 
  apply(1, function(x) x/sum(x)) %>% t

prior_sigma <- 
  list(sigma = list(a = alpha/p, b = 1/2-alpha/p))
init_T <- 
  list(T = rowSums(N)/100)

res_MFVI <- 
  main_MFVI(
    N, Sigma_obs, 
    prior = prior_sigma,
    init = init_T, 
    iter_max = 1e4,
    iter_crit = 1e-4,
    method = "MFVI",
    verbose = TRUE, 
    # debug parameters
    par_to_update = c("T", "sigma", "Q", "G"),
    early_termination = TRUE, 
    early_termination_crit = 1e-4)

save(res_MFVI, file = paste0(file_addr_MFVI, "res_MFVI.RData"))
load(paste0(file_addr_MFVI, "res_MFVI.RData"))

lambda <- res_MFVI$lambda
info <- res_MFVI$info
iter <- res_MFVI$iter

plot_diag <- TRUE
if (plot_diag){
  iter_idx <- 1:(info$oper$iter-1)
  for (res_type in c("par", "E")){
    for (res_name in c("Q1", "Q2", "T", "sigma")){
      pdf(paste0(file_addr_MFVI, res_type, "_", res_name, ".pdf"))
      # graphically examing change in mean_Q2 over iteration
      data_plot_raw <- iter[[res_type]][[res_name]]
      
      if (length(dim(data_plot_raw)) == 3){
        data_plot <- 
          data_plot_raw[iter_idx, , ]
      } else if (length(dim(data_plot_raw)) == 2) {
        data_plot <- 
          data_plot_raw[iter_idx, ]
      }
      
      iter_cur <- dim(data_plot)[1]
      plot(0, type = "n", 
           xlim = c(1, iter_cur), 
           ylim = range(data_plot), 
           main = paste(res_type, res_name)
      )
      
      if (length(dim(data_plot)) == 3){
        for (j in 1:info$stat$J){
          for (i in 1:info$stat$I){
            lines(1:iter_cur, 
                  data_plot[, j, i], 
                  col = 
                    c(rgb(0,0,0,0.3), 
                      rgb(1,0,0,0.3)
                    )[2 - (N[j, i]>0)]
            )
          }
        }
      } else if (length(dim(data_plot)) == 2) {
        for (j in 1:ncol(data_plot)){
          lines(1:iter_cur, data_plot[, j], 
                col = rgb(0,0,0,0.5))
        }
      }
      dev.off()
    }
  }
}


#### 3.1. Reconstruction result ####
pred_dist(info, N)

P_j_raw <- 
  (info$mean$sigma * t(info$mean$Q2)) %>% t
P_j <- P_j_raw/rowSums(P_j_raw)

P_emp <- N/rowSums(N)

mean(as.matrix(abs(P_j - P_emp)))


# covariance estimate
cov_est <- (t(data.sim$Y.tru) %*% data.sim$Y.tru + 
              diag(rep(data.sim$er, ncol(data.sim$Y.tru)))) 
cov_mvi <- #iter$E$Q1[459, , ] %>% t %>% cor
  info$mean$Q1 %>% t %>% cor
  
heatmap(cov2cor(cov_est))
heatmap(info$mean$Q1 %>% t %>% cor)


heatmap((cov2cor(cov_est) - cov_mvi)/cov2cor(cov_est))


# ELBO assessment
elbo_plot <- iter$crit$ELBO[
  is.finite(iter$crit$ELBO)]
plot(elbo_plot, type = "l", 
     xlab = "Iterations", ylab = "Unnormalized ELBO")
