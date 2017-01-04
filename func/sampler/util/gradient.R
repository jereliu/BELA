# conditional gradient for HMC

grad_hmc_U <- function(T_suff, U, V, lambda, d1){
  (-T_suff + d1(U %*% t(V))) %*% V + lambda * U
}

grad_hmc_V <- function(T_suff, U, V, lambda, d1){
  t(-T_suff + d1(U %*% t(V))) %*% U + lambda * V
}

lik_target <- function(U, V, T_suff, lambda, dist_family){
  dist_family$negloglik(T_suff, U %*% t(V)) + 
    (lambda/2) * (sum(U^2) + sum(V^2))
}

lik_target_S <- function(S, T_suff, lambda, dist_family){
  U <- S[1:nrow(T_suff), ]
  V <- S[(nrow(T_suff)+1):nrow(S), ]
  
  dist_family$negloglik(T_suff, U %*% t(V)) + 
    (lambda/2) * (sum(U^2) + sum(V^2))
}

# Total gradient wrt negative log likelihood and wrt KSD, for SVGD
grad_neglik_S <- 
  function(S, T_suff, lambda, dist_family){
    # S: (n + p) x k vector
    # T_suff: n x p vector
    n <- nrow(T_suff)
    p <- ncol(T_suff)
    k <- ncol(S)
    
    I_U <- rbind(diag(n), matrix(0, nrow = p, ncol = n))
    I_V <- rbind(matrix(0, nrow = n, ncol = p), diag(p)) 
    U <- S[1:n, ]
    V <- S[(n+1):(n+p), ]
    
    #
    T_mat <- I_U %*% T_suff %*% t(I_V)
    comp1 <- - T_mat - t(T_mat)
    
    # 
    dAS_mat <- dist_family$partition$d(U %*% t(V))
    comp2 <- 
      rbind(
        cbind(matrix(0, nrow = n, ncol = n), dAS_mat),
        cbind(t(dAS_mat), matrix(0, nrow = p, ncol = p))
      )
    
    # 
    (comp1 + comp2 + lambda * diag(n+p)) %*% S
  }

# gradient of KL-divergence (which is KSD)
grad_ksd <- 
  function(S, alpha = -1, 
           T_suff, lambda, dist_family){
    # S <- array(rnorm(n_sample * (n + p) * k),
    #            dim = c(n_sample, (n + p)*k))
    n_sample <- dim(S)[1]
    n <- nrow(T_suff)
    p <- ncol(T_suff)
    
    # produce kernel matrix
    SS <- S %*% t(S)
    S2 <- array(rowSums(S^2), dim = c(n_sample, 1))
    S2e <- repmat(S2, 1, n_sample)
    H <- (S2e + t(S2e) - 2 * SS)
    if (alpha == -1){
      alpha <- log(n)/(median(H)^2)
    }
    K <- exp(-alpha * H/2)
    
    # calculate negative log likelihood
    neg_grad_S <- 
      apply(S, 1, 
            function(S_i)
              matrix(S_i, nrow = n+p) %>%
              grad_neglik_S(T_suff, lambda, dist_family) %>%
              as.vector
      ) %>% t
    
    # calculate ksd gradient
    grad_out <-
      (K %*% (-neg_grad_S - alpha * S) + 
      alpha * rowSums(K) * S)/n_sample
    
    # output 
    list(grad = grad_out, alpha = alpha)
  }