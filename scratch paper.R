#### 1. calculate mean for each Q_ij, Q_ij^2 ####

mode_est <- info$plot$Q1$mode_est
normalizer <- info$plot$Q1$normalizer
n_ij <- info$stat$n_ij
s_j <- prior$Q$Sg_cond
sd <- sqrt(s_j)
lambda_Q1 <- lambda_prev$Q1
lambda_Q2 <- lambda_prev$Q2

mean <- info$mean$Q1
file_pref <- "./temp_plot/meanQ2 "

round(info$mean$Q2, 3)
apply(data$Q, 1, function(r) pmax(0, r))^2 %>% round(3)

for (i in 1:I){
  for (j in 1:J){
    pdf(paste0(file_pref, " ", j, "-", i, 
               " N = ", n_ij[j, i], ".pdf"))
    if (N[j, i] == 0){
      Q <- seq(info$mean$Q1[j, i] - 5 *sd[j], 
               info$mean$Q1[j, i] + 5 *sd[j], 1e-3)
      kernel_val <- 
        Q_kernel_v0(
          Q,
          lambda_Q1 = lambda$Q1[j, i],
          lambda_Q2 = lambda$Q2[j, i],
          n_ij = n_ij[j, i],
          s_j = s_j[j],
          const_adj = 1)
      idx <- is.finite(kernel_val)
      
      plot(pmax(0, Q)[idx]^2,
           kernel_val[idx]/max(kernel_val[idx]),
           type = "l", ylab = "", 
           main = paste0("Q: ", j, "-", i, 
                         " N = ", n_ij[j, i]))
    } else {
      Q <- seq(1e-10, 
               mode_est[j, i] + 10*sd[j], 1e-3)
      plot(Q, #pmax(0, Q)^2,
           Q_kernel(
             pmax(0, Q)^2,
             lambda_Q1 = lambda_Q1[j, i],
             lambda_Q2 = lambda_Q2[j, i],
             n_ij = n_ij[j, i],
             s_j = s_j[j],
             log_const_adj = normalizer[j, i]),
           type = "l", ylab = "", 
           main = paste0("Q: ", j, "-", i, 
                         " N = ", n_ij[j, i]))
    }
    abline(v = mean[j, i], col = 2)
    dev.off()
  }
}

#### 2. sigma mean ####
eps <- 1e-10
a <- 2
b <- 10
lambda <- 15

E_sig <- function(power_x, power_1_x, a, b, lambda){
  nom <- 
    integral(fun = sigma_kernel, 
             xmin = eps, xmax = 1 - eps, 
             lambda = lambda, 
             alpha = a + power_x, 
             beta = b + power_1_x)
  denom <- 
    integral(fun = sigma_kernel, 
             xmin = eps, xmax = 1 - eps, 
             lambda = lambda, 
             alpha = a, beta = b
    )
  nom/denom
}
  
E_x_1_x <- E_sig(1, 1, a, b, lambda)
E_1_x <- E_sig(0, 1, a, b, lambda)
E_x <- E_sig(1, 0, a, b, lambda)
E_x2 <- E_sig(2, 0, a, b, lambda)
denom <- 
  integral(fun = sigma_kernel, 
           xmin = eps, xmax = 1 - eps, 
           lambda = lambda, 
           alpha = a - 1, beta = b - 1
  )

# integration by part formula
- a * E_1_x + b * E_x + lambda * E_x_1_x

- a + (a + b + lambda) * E_x  - lambda * E_x2

# -a * Phi(lambda, a, b) + 
# (a + b + lambda) * Phi(lambda, a+1, b) - 
# lambda * Phi(lambda, a+2, b) = 0


cov_est <- (t(data.sim$Y.tru) %*% data.sim$Y.tru + 
              diag(rep(data.sim$er, ncol(data.sim$Y.tru)))) 
cov2cor(cov_est) %>% image

info$mean$Q1 %>% t %>% cor %>% image
