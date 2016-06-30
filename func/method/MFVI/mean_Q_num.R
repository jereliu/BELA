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
  
  ratio_thresh <- info$oper$Q_change_thresh
  # estimate mode/normalizer to 'tame' kernel function
  S_j_mat <- matrix(rep(s_j, I), nrow = J)
  lambda_S_1 <- 
    1 + 2 * lambda_Q2 * S_j_mat
  
  info$plot$Q1$mode_est <- mode_est <- 
    (lambda_Q1 + 
       sqrt(lambda_Q1^4 + 
              64*n_ij*(S_j_mat^4)*lambda_S_1 + 
              32 * n *(S_j_mat^2)))/
    (8 * lambda_S_1 * (S_j_mat^2) + 4)
  
  if (finer_mode){
    # find out value & position of mode
    info$plot$Q1$mode_est <- mode_est <-
      # grid search for mode and mode value
      outer(1:J, 1:I,
            Vectorize(
              function(j, i){
                # TODO: adaptive Q range
                Q <- seq(0, mode_est[j, i] + 50 * sd[j], 1e-2)                
                Q_resp <- 
                  Q_kernel(
                    Q,
                    lambda_Q1 = lambda_Q1[j, i],
                    lambda_Q2 = lambda_Q2[j, i],
                    n_ij = n_ij[j, i],
                    s_j = s_j[j], 
                    log = TRUE
                  )
                
                finite_idx <- is.finite(Q_resp)
                Q <- Q[finite_idx]
                Q_resp <- Q_resp[finite_idx]
                
                Q[which.max(Q_resp)]
              }
            ))
    
    # find out value & position of mode
    info$plot$Q1$normalizer <- normalizer <-
      # grid search for mode and mode value
      outer(1:J, 1:I,
            Vectorize(
              function(j, i){
                Q_kernel(
                  mode_est[j, i],
                  lambda_Q1 = lambda_Q1[j, i],
                  lambda_Q2 = lambda_Q2[j, i],
                  n_ij = n_ij[j, i],
                  s_j = s_j[j], 
                  log = TRUE
                )
              }
            ))
    
    info$plot$Q1$range_min <- range_min <-
      # grid search for mode and mode value
      outer(1:J, 1:I,
            Vectorize(
              function(j, i){
                # TODO: adaptive Q range
                Q <- 
                  seq(0, mode_est[j, i] + 50 * sd[j],
                      1e-2)
                
                Q_resp <- 
                  Q_kernel(
                    Q,
                    lambda_Q1 = lambda_Q1[j, i],
                    lambda_Q2 = lambda_Q2[j, i],
                    n_ij = n_ij[j, i],
                    s_j = s_j[j], 
                    log_const_adj = normalizer[j, i]
                  )
                
                range_idx <- which(Q_resp > 1e-10)
                finite_idx <- which(is.finite(Q_resp))
                final_idx <- intersect(range_idx, finite_idx)
                
                Q <- Q[final_idx]
                Q_resp <- Q_resp[final_idx]
                
                min(Q) - sd[j]
              }
            ))
    
    info$plot$Q1$range_max <- range_max <-
      # grid search for mode and mode value
      outer(1:J, 1:I,
            Vectorize(
              function(j, i){
                # TODO: adaptive Q range
                Q <- 
                  seq(0, mode_est[j, i] + 50 * sd[j],
                      1e-2)
                
                Q_resp <- 
                  Q_kernel(
                    Q,
                    lambda_Q1 = lambda_Q1[j, i],
                    lambda_Q2 = lambda_Q2[j, i],
                    n_ij = n_ij[j, i],
                    s_j = s_j[j], 
                    log_const_adj = normalizer[j, i]
                  )
                
                range_idx <- which(Q_resp > 1e-10)
                finite_idx <- which(is.finite(Q_resp))
                final_idx <- intersect(range_idx, finite_idx)
                
                Q <- Q[final_idx]
                Q_resp <- Q_resp[final_idx]
                
                max(Q) + sd[j]
              }
            ))
  }
  
  if (plot){
    # find max count element
    coord <- which(n_ij == max(n_ij), arr.ind = TRUE)
    j <- coord[1]
    i <- coord[2]
    
    # plot the mofo
    Q <- seq(range_min[j, i], range_max[j, i], 1e-3)
    
    plot(Q,
         Q_kernel(
           Q,
           lambda_Q1 = lambda_Q1[j, i],
           lambda_Q2 = lambda_Q2[j, i],
           n_ij = n_ij[j, i],
           s_j = s_j[j],
           log_const_adj = normalizer[j, i]),
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
                         xmin = 1e-10,
                         xmax = range_max[j, i],
                         lambda_Q1 = lambda_Q1[j, i],
                         lambda_Q2 = lambda_Q2[j, i],
                         n_ij = n_ij[j, i],
                         s_j = s_j[j],
                         log_const_adj = normalizer[j, i])
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
                         xmin = range_min[j, i], 
                         xmax = range_max[j, i], 
                         lambda_Q1 = lambda_Q1[j, i], 
                         lambda_Q2 = lambda_Q2[j, i],
                         n_ij = n_ij[j, i], 
                         s_j = s_j[j], 
                         log_const_adj = normalizer[j, i])
              }
            }
          )
    )
  
  # return expectation
  res <- nom/denom
  if (plot) 
    abline(v = res[j, i], col = 2)
  
  # check if result makes sense
  if (!is.null(info$mean$Q1)){
    ratio_change <- 
    #  abs((info$mean$Q1 - res)/(info$mean$Q1)) %>% max(na.rm = TRUE)
      abs((info$mean$Q1 - res)) %>% max(na.rm = TRUE)
    if(ratio_change > ratio_thresh){
      stop("some Q1 value grow ", round(ratio_change, 3), " times (>threshold = ", ratio_thresh, ") in this iteration")
    }
  }
  # return
  info$mean$Q1 <- res
  info
}

mean_Q2 <- 
  function(lambda, info, prior, 
           finer_mode = TRUE,
           plot = FALSE, 
           plot_diag = FALSE,
           verbose)
  {
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
    ratio_thresh <- info$oper$Q_change_thresh
    
    lambda_Q1 <- lambda$Q1
    lambda_Q2 <- lambda$Q2
    
    # estimate mode
    S_j_mat <- matrix(rep(s_j, I), nrow = J)
    lambda_S_1 <- 
      1 + 2 * lambda_Q2 * S_j_mat
    
    info$plot$Q1$mode_est <- mode_est <- 
      (lambda_Q1 + 
         sqrt(lambda_Q1^4 + 
                64*n_ij*(S_j_mat^4)*lambda_S_1 + 
                32 * n *(S_j_mat^2)))/
      (8 * lambda_S_1 * (S_j_mat^2) + 4)
    
    if (finer_mode){
      # find out value & position of mode
      info$plot$Q1$mode_est <- mode_est <-
        # grid search for mode and mode value
        outer(1:J, 1:I,
              Vectorize(
                function(j, i){
                  # TODO: adaptive Q range
                  Q <- 
                    seq(0, mode_est[j, i] + 50 * sd[j],
                        1e-2)
                  
                  Q_resp <- 
                    Q_kernel(
                      Q,
                      lambda_Q1 = lambda_Q1[j, i],
                      lambda_Q2 = lambda_Q2[j, i],
                      n_ij = n_ij[j, i],
                      s_j = s_j[j], 
                      log = TRUE
                    )
                  
                  finite_idx <- is.finite(Q_resp)
                  Q <- Q[finite_idx]
                  Q_resp <- Q_resp[finite_idx]
                  
                  Q[which.max(Q_resp)]
                }
              ))
      
      # find out value & position of mode
      info$plot$Q1$normalizer <- normalizer <-
        # grid search for mode and mode value
        outer(1:J, 1:I,
              Vectorize(
                function(j, i){
                  Q_kernel(
                    mode_est[j, i],
                    lambda_Q1 = lambda_Q1[j, i],
                    lambda_Q2 = lambda_Q2[j, i],
                    n_ij = n_ij[j, i],
                    s_j = s_j[j], 
                    log = TRUE
                  )
                }
              ))
      
      info$plot$Q1$range_min <- range_min <-
        # grid search for mode and mode value
        outer(1:J, 1:I,
              Vectorize(
                function(j, i){
                  # TODO: adaptive Q range
                  Q <- 
                    seq(0, mode_est[j, i] + 50 * sd[j],
                        1e-2)
                  
                  Q_resp <- 
                    Q_kernel(
                      Q,
                      lambda_Q1 = lambda_Q1[j, i],
                      lambda_Q2 = lambda_Q2[j, i],
                      n_ij = n_ij[j, i],
                      s_j = s_j[j], 
                      log_const_adj = normalizer[j, i]
                    )
                  
                  range_idx <- which(Q_resp > 1e-10)
                  finite_idx <- which(is.finite(Q_resp))
                  final_idx <- intersect(range_idx, finite_idx)
                  
                  Q <- Q[final_idx]
                  Q_resp <- Q_resp[final_idx]
                  
                  min(Q) - sd[j]
                }
              ))
      
      # find out value & position of mode
      info$plot$Q1$range_max <- range_max <-
        # grid search for mode and mode value
        outer(1:J, 1:I,
              Vectorize(
                function(j, i){
                  # TODO: adaptive Q range
                  Q <- 
                    seq(0, mode_est[j, i] + 50 * sd[j],
                        1e-2)
                  
                  Q_resp <- 
                    Q_kernel(
                      Q,
                      lambda_Q1 = lambda_Q1[j, i],
                      lambda_Q2 = lambda_Q2[j, i],
                      n_ij = n_ij[j, i],
                      s_j = s_j[j], 
                      log_const_adj = normalizer[j, i]
                    )
                  
                  range_idx <- which(Q_resp > 1e-10)
                  finite_idx <- which(is.finite(Q_resp))
                  final_idx <- intersect(range_idx, finite_idx)
                  
                  Q <- Q[final_idx]
                  Q_resp <- Q_resp[final_idx]
                  
                  max(Q) + sd[j]
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
      Q <- seq(range_min[j, i], range_max[j, i], 1e-3)
      plot(pmax(0, Q)^2,
           Q_kernel(
             Q,
             lambda_Q1 = lambda_Q1[j, i],
             lambda_Q2 = lambda_Q2[j, i],
             n_ij = n_ij[j, i],
             s_j = s_j[j],
             log_const_adj = normalizer[j, i]
           ),
           type = "l", ylab = "", 
           main = paste0("|Q|_+^2 ", j, "-", i))
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
                           xmin = range_min[j, i], 
                           xmax = range_max[j, i], 
                           lambda_Q1 = lambda_Q1[j, i], 
                           lambda_Q2 = lambda_Q2[j, i],
                           n_ij = n_ij[j, i] + 1, 
                           s_j = s_j[j],
                           log_const_adj = normalizer[j, i])
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
                           xmin = range_min[j, i], 
                           xmax = range_max[j, i], 
                           lambda_Q1 = lambda_Q1[j, i], 
                           lambda_Q2 = lambda_Q2[j, i],
                           n_ij = n_ij[j, i], 
                           s_j = s_j[j], 
                           log_const_adj = normalizer[j, i]
                  )
                }
              }
            )
      )
    
    # return expectation
    res <- nom/denom
    if (plot)       
      abline(v = res[j, i], col = 2)
    
    if (plot_diag){
      file_dir <- "./temp_plot/"
      # graphically examing mean for Q2 est 
      for (j in 1:J){
        for (i in 1:I){
          pdf(paste0(file_dir, "|Q|_+^2 ", j, "-", i, ", N = ", N[j, i], ".pdf"))
          Q <- seq(0, mode_est[j, i]+ 50*sd[j], 1e-3)
          kern_val <- 
            Q_kernel(
              Q,
              lambda_Q1 = lambda_Q1[j, i],
              lambda_Q2 = lambda_Q2[j, i],
              n_ij = n_ij[j, i],
              s_j = s_j[j],
              log_const_adj = normalizer[j, i]
            )
          idx <- which(kern_val > 1e-10)
          plot(pmax(0, Q)[idx]^2,
               kern_val[idx],
               type = "l", ylab = "", 
               main = paste0("|Q|_+^2 ", j, "-", i, ", N = ", N[j, i]))
          abline(v = res[j, i], col = 2)
          dev.off()
        }
      }
    }
    
    # check if result makes sense
    if (!is.null(info$mean$Q2)){
      ratio_change <- 
        abs((info$mean$Q2 - res)/(info$mean$Q2)) %>% max(na.rm = TRUE)
      if(ratio_change > ratio_thresh){
        stop("some Q2 value grow ", round(ratio_change, 3), " times (>threshold = ", ratio_thresh, ") in this iteration")
      }
    }
    
    # return
    info$mean$Q2 <- res
    info
  }



#### 2. posterior kernel function ####
Q_kernel <- 
  function(Q_ij, lambda_Q1, lambda_Q2, 
           n_ij, s_j, log_const_adj = 0, 
           log = FALSE){
    
    log_out <- 
      (2*n_ij) * log(pmax(Q_ij, 0)) -
      lambda_Q2 * pmax(Q_ij, 0)^2 - 
      0.5 * Q_ij^2 /s_j + lambda_Q1 * Q_ij /s_j - 
      log_const_adj
    
    if (log) {
      return(log_out)
    } else {
      return(exp(log_out))
    }
  }

Q_kernelxQ <- 
  function(Q_ij, lambda_Q1, lambda_Q2, 
           n_ij, s_j, log_const_adj = 0, 
           log = FALSE){
    log_out <- 
      log(Q_ij) + (2*n_ij) * log(pmax(Q_ij, 0)) -
      lambda_Q2 * pmax(Q_ij, 0)^2 - 
      0.5 * Q_ij^2 /s_j + lambda_Q1 * Q_ij /s_j - 
      log_const_adj
    
    if (log) {
      return(log_out)
    } else {
      return(exp(log_out))
    }
  }

Q_kernelxQ_n <- 
  function(Q_ij, lambda_Q1, lambda_Q2, 
           n_ij, s_j, log_const_adj = 0, 
           n = 2, log = FALSE){
    log_out <- 
      n * log(Q_ij) + (2*n_ij) * log(pmax(Q_ij, 0)) -
      lambda_Q2 * pmax(Q_ij, 0)^2 - 
      0.5 * Q_ij^2 /s_j + lambda_Q1 * Q_ij /s_j - 
      log_const_adj
    
    if (log) {
      return(log_out)
    } else {
      return(exp(log_out))
    }
  }
