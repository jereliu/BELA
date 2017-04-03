# iter <- 1e2
# Theta <- rec$Theta[iter, , ]
# U <- rec$U[iter, , ]
# V <- rec$V[iter, , ] 

# TODO: find derivative with respect to Theta
score_svd <- # derivative of likelihood wrt singular valuez
  function(Theta, U, V, T_suff, lambda, 
           dist_family, n_eigen = NULL, 
           log = TRUE
  ){
    # extract svd component
    Theta_svd <- svd(Theta)
    if (is.null(n_eigen)) n_eigen <- sum(Theta_svd$d > 1e-5)
    comp_svd <- # e^l * t(e^r)
      lapply(1:n_eigen, function(i) 
        Theta_svd$u[, i] %*% t(Theta_svd$v[, i]) )
    
    # compute model component
    A_deriv <- dist_family$partition$d
    comp_model <- 
      -T_suff + A_deriv(Theta) #+ 
      #2*lambda* (U %*% t(1/V) + (1/U) %*% t(V))
    
    # compute score output, by default use log
    sapply(1:n_eigen, function(i) sum(comp_svd[[i]] * comp_model))
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
      alpha <- (median(H))/log(n_sample+1)
    }
    K <- exp(-H/alpha)
    
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

# support function for KSD
repmat <- 
  function (X, a = 1, b = 1){
    rows <- dim(X)[1]
    cols <- dim(X)[2]
    if (is.null(cols)) 
      cols <- 1
    rowRep <- matrix(rep(t(X), a), ncol = cols, byrow = TRUE)
    newX <- matrix(rep(rowRep, b), ncol = cols * b)
    return(newX)
  }

