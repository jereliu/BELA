source("../DirichletFactor/MCMC_fns.R")
source("../DirichletFactor/utilities.R")
source("../DirichletFactor/simulation_study/contour_sim.R")

require(magrittr)
require(dplyr)

source("./func/util/source_Dir.R")
sourceDir("./func")
sourceDir("./func_MFVI")



#### 1. Brute-Force MF Variational Inference ####
alpha = 10
n = 22
p = 68
sigma.value = seq(0.001,0.999,0.001)
tmp = c(0, pbeta( sigma.value, alpha/p, 1/2-alpha/p ) )
sigma.prior = sapply( 1:999, function(x) tmp[x+1]-tmp[x] )
sigma.prior[length(sigma.prior)] = sigma.prior[length(sigma.prior)] + 1-sum(sigma.prior)

data <- 
  sim.for.contour(1000, n = 20, p = 10, 
                  m = 3, 
                  K = 2, 
                  a.er = 1, 
                  b.er = 0.3, 
                  sigma.value,
                  sigma.prior, 
                  strength = 3 )

data <- 
  biom_sample(n_sample = 20, n_OTU = 50, 
              n_block = 2, n_obs = 50)

N <- data$N
Sigma <- data$Sigma

res <- main_MFVI(
  N, Sigma, 
  prior = NULL,
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
