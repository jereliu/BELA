# method to calculate E(Q) and E(|Q|^2) based on numeric integreation
# used when use 'exact' variational family that's in the same family as original p

mean_Q1 <- function(
  lambda, info, prior, 
  adj_power = 2, plot = FALSE, 
  verbose, finer_mode = TRUE)
{
  if (verbose) {
    cat("..E(Q)..")
    plot <- TRUE
  }
  
  # brute force numeric mean for Q
  # adj_power: number used to tune the integration procedure
  I <- info$stat$I
  J <- info$stat$J 
  n_ij <- info$stat$n_ij
  s_j <- prior$Q$Sg_cond
  sd <- sqrt(s_j)
  
  lambda_Q1 <- lambda$Q1
  lambda_Q2 <- lambda$Q2
  
  # estimate mode/normalizer to 'tame' kernel function
  S_j_mat <- matrix(rep(s_j, I), nrow = J)
  lambda_S_1 <- 
    1 + 2 * lambda_Q2 * S_j_mat
  
  info$plot$Q1$mode_est <- mode_est <- 
    (lambda_Q1 + sqrt(lambda_Q1^2 + 16*n_ij*S_j_mat*lambda_S_1))/
    (4 * lambda_S_1)
  
  if (finer_mode){
    # find out value & position of mode
    info$plot$Q1$normalizer <- normalizer <- 
      # grid search for mode and mode value
      outer(1:J, 1:I,
            Vectorize(
              function(j, i){
                Q <- seq(
                  mode_est[j, i] - 10*sd[j], 
                  mode_est[j, i] + 10*sd[j], 1e-2)
                
                Q_resp <- 
                  Q_kernel(
                    Q,
                    lambda_Q1 = lambda_Q1[j, i],
                    lambda_Q2 = lambda_Q2[j, i],
                    n_ij = n_ij[j, i],
                    s_j = s_j[j])
                
                #Q[which.max(Q_resp)]
                max(Q_resp)
              }
            ))
    
    info$plot$Q1$mode_est <- mode_est <-
      outer(1:J, 1:I,
            Vectorize(
              function(j, i){
                Q <- seq(
                  mode_est[j, i] - 10*sd[j], 
                  mode_est[j, i] + 10*sd[j], 1e-2)
                
                Q_resp <- 
                  Q_kernel(
                    Q,
                    lambda_Q1 = lambda_Q1[j, i],
                    lambda_Q2 = lambda_Q2[j, i],
                    n_ij = n_ij[j, i],
                    s_j = s_j[j])
                
                Q[which.max(Q_resp)]
              }
            ))
  }
  
  if (plot){
    # find max count element
    coord <- which(n_ij == max(n_ij), arr.ind = TRUE)
    j <- coord[1]
    i <- coord[2]
    
    # plot the mofo
    Q <- seq(mode_est[j, i] - 10*sd[j], 
             mode_est[j, i] + 10*sd[j], 1e-3)
    
    plot(Q,
         Q_kernel(
           Q,
           lambda_Q1 = lambda_Q1[j, i],
           lambda_Q2 = lambda_Q2[j, i],
           n_ij = n_ij[j, i],
           s_j = s_j[j],
           const_adj = normalizer[j, i]),
         type = "l", ylab = "", 
         main = paste0("Q"))
  }
  
  nom <-
    outer(1:J, 1:I,
          Vectorize(
            function(j, i){
              if (n_ij[j, i] == 0){
                mean_Q1_no_n(lambda_Q1 = lambda_Q1[j, i],
                            lambda_Q2 = lambda_Q2[j, i],
                            s_j = s_j[j])
              } else {
                integral(fun = Q_kernelxQ,
                         xmin =
                           (1-(n_ij[j, i] == 0))*
                           (mode_est[j, i] - 10*sd[j]),
                         xmax = mode_est[j, i] + 10*sd[j],
                         lambda_Q1 = lambda_Q1[j, i],
                         lambda_Q2 = lambda_Q2[j, i],
                         n_ij = n_ij[j, i],
                         s_j = s_j[j],
                         const_adj = normalizer[j, i])
              }
            }
          )
    )
  
  denom <- 
    outer(1:J, 1:I,
          Vectorize(
            function(j, i){
              if (n_ij[j, i] == 0){
                1 
              } else {
                integral(fun = Q_kernel, 
                         xmin = 
                           (1-(n_ij[j, i] == 0))* 
                           (mode_est[j, i] - 10*sd[j]), 
                         xmax = mode_est[j, i] + 10*sd[j], 
                         lambda_Q1 = lambda_Q1[j, i], 
                         lambda_Q2 = lambda_Q2[j, i],
                         n_ij = n_ij[j, i], 
                         s_j = s_j[j], 
                         const_adj = normalizer[j, i])
              }
            }
          )
    )
  
  # return expectation
  res <- nom/denom
  if (plot) 
    abline(v = res[j, i], col = 2)
  
  # return
  info$mean$Q1 <- res
  info
}

mean_Q2 <- 
  function(lambda, info, prior, 
           adj_power = 2,
           plot = FALSE, 
           finer_mode = TRUE,
           verbose){
    if (verbose) {
      cat("..E(|Q|+^2)..")
      plot <- TRUE
    }
    
    # brute force numeric mean for |Q|+^2
    I <- info$stat$I
    J <- info$stat$J     
    n_ij <- info$stat$n_ij
    s_j <- prior$Q$Sg_cond
    sd <- sqrt(s_j)
    
    lambda_Q1 <- lambda$Q1
    lambda_Q2 <- lambda$Q2
    
    # estimate mode
    S_j_mat <- matrix(rep(s_j, I), nrow = J)
    lambda_S_1 <- 
      1 + 2 * lambda_Q2 * S_j_mat
    
    mode_est <- 
      (lambda_Q1 + sqrt(lambda_Q1^2 + 8*n_ij*S_j_mat*lambda_S_1))/
      (2 * lambda_S_1)
    
    if (finer_mode){
      # find out value & position of mode
      normalizer <- 
        # grid search for mode and mode value
        outer(1:J, 1:I,
              Vectorize(
                function(j, i){
                  Q <- seq(
                    mode_est[j, i] - 10*sd[j], 
                    mode_est[j, i] + 10*sd[j], 1e-2)
                  
                  Q_resp <- 
                    Q_kernel(
                      Q,
                      lambda_Q1 = lambda_Q1[j, i],
                      lambda_Q2 = lambda_Q2[j, i],
                      n_ij = n_ij[j, i],
                      s_j = s_j[j])
                  
                  #Q[which.max(Q_resp)]
                  max(Q_resp)
                }
              ))
      mode_est <- 
        outer(1:J, 1:I,
              Vectorize(
                function(j, i){
                  Q <- seq(
                    mode_est[j, i] - 10*sd[j], 
                    mode_est[j, i] + 10*sd[j], 1e-2)
                  
                  Q_resp <- 
                    Q_kernel(
                      Q,
                      lambda_Q1 = lambda_Q1[j, i],
                      lambda_Q2 = lambda_Q2[j, i],
                      n_ij = n_ij[j, i],
                      s_j = s_j[j])
                  
                  Q[which.max(Q_resp)]
                }
              ))
    }
    
    # # # # # numeric sanity check
    if (plot){
      # find max count element
      coord <- which(n_ij == max(n_ij), arr.ind = TRUE)
      j <- coord[1]
      i <- coord[2]
      
      # plot the mofo
      Q <- seq(0, mode_est[j, i]+ 10*sd[j], 1e-3)
      plot(pmax(0, Q)^2,
           Q_kernel(
             Q,
             lambda_Q1 = lambda_Q1[j, i],
             lambda_Q2 = lambda_Q2[j, i],
             n_ij = n_ij[j, i],
             s_j = s_j[j],
             const_adj = normalizer[j, i]
           ),
           type = "l", ylab = "", 
           main = paste0("|Q|_+^2"))
    }
    
    #### Calculation
    #
    nom <- 
      outer(1:J, 1:I,
            Vectorize(
              function(j, i){
                if (n_ij[j, i] == 0){
                  mean_Q2_no_n(lambda_Q1 = lambda_Q1[j, i],
                               lambda_Q2 = lambda_Q2[j, i],
                               s_j = s_j[j])
                } else {                
                  integral(fun = Q_kernel, 
                           xmin = 0, 
                           xmax = mode_est[j, i] + 10*sd[j], 
                           lambda_Q1 = lambda_Q1[j, i], 
                           lambda_Q2 = lambda_Q2[j, i],
                           n_ij = n_ij[j, i] + 1, 
                           s_j = s_j[j],
                           const_adj = normalizer[j, i])
                }
              }
            )
      ) 
    
    denom <- 
      outer(1:J, 1:I,
            Vectorize(
              function(j, i){
                if (n_ij[j, i] == 0){
                  1 
                } else {
                  integral(fun = Q_kernel, 
                           xmin = mode_est[j, i] - 10*sd[j], 
                           xmax = mode_est[j, i] + 10*sd[j], 
                           lambda_Q1 = lambda_Q1[j, i], 
                           lambda_Q2 = lambda_Q2[j, i],
                           n_ij = n_ij[j, i], 
                           s_j = s_j[j], 
                           const_adj = normalizer[j, i])
                }
              }
            )
      )
    
    
    # return expectation
    res <- nom/denom
    if (plot)       
      abline(v = res[j, i], col = 2)
    
    # return 
    res
  }



#### 2. posterior kernel function ####
Q_kernel <- 
  function(Q_ij, lambda_Q1, lambda_Q2, 
           n_ij, s_j, const_adj = 1){
    
    (pmax(Q_ij, 0))^(2*n_ij) /const_adj * 
      exp(-lambda_Q2 * pmax(Q_ij, 0)^2 - 
            0.5 * Q_ij^2 /s_j + lambda_Q1 * Q_ij /s_j)
  }

Q_kernelxQ <- 
  function(Q_ij, lambda_Q1, lambda_Q2, 
           n_ij, s_j, const_adj = 1){
    
    Q_ij * (pmax(Q_ij, 0))^(2*n_ij) /const_adj * 
      exp(-lambda_Q2 * pmax(Q_ij, 0)^2 - 
            0.5 * Q_ij^2 /s_j + lambda_Q1 * Q_ij/s_j)
  }