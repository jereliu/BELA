pred_dist <- 
  function(info, N){
    pred_mat <- 
      (info$mean$sigma * t(info$mean$Q2)) %>% 
      t %>% apply(1, function(x) x/sum(x))
    orig_mat <- 
      apply(N, 1, function(x) x/sum(x))
    orig_mat_norm <- 
      orig_mat %>% abs %>% mean
    
    (orig_mat-pred_mat) %>% 
      divide_by(orig_mat_norm) %>% abs %>% mean
  }

ELBO_full <- 
  function(lambda, prior, info){
    # computing overall ELBO in MFVI, 
    # check appendix for documentation
    
    # calc var(Q) under variational distribution
    I <- info$stat$I
    J <- info$stat$J 
    n_ij <- info$stat$n_ij
    s_j <- prior$Q$Sg_cond
    sd <- sqrt(s_j)
    
    info$stat$Q_var <- 
      outer(1:J, 1:I,
            Vectorize(
              function(j, i){
                if (n_ij[j, i] == 0){
                  var_Q1_no_n(lambda_Q1 = lambda_Q1[j, i],
                              lambda_Q2 = lambda_Q2[j, i],
                              s_j = s_j[j])
                } else {                
                  info$mean$Q2[j, i] - info$mean$Q1[j, i]^2
                }
              }
            )
      )
    
    #
    term_1 <- 
      2 * (info$stat$TxSigma * info$mean$Q2) %>% sum
    
    term_2 <- 
      lapply(1:I, 
             function(i)
               sum(
                 as.vector(prior$Q$Sigma_inv) * 
                   as.vector(diag(info$stat$Q_var[, i]))
               )
      ) %>% do.call(sum, .)
    
    term_3 <- 
      ((0.5 * info$mean$Q2 - info$mean$Q1^2)/s_j) %>%
      sum
    
    term_1 - term_2 + term_3
  }


ELBO_cond <- 
  function(lambda, prior, info, 
           par_name = "T",verbose = FALSE){
    # compute conditional ELBO
    if (par_name == "T"){
      par_q <- lambda$T
      par_p <- 
        # compute E(sigma) * E(|Q_ij|^2)
        (info$mean$sigma * t(info$mean$Q2)) %>% 
        # add over OTU
        colSums
      n_j <- info$stat$n_j
      E_T <- info$mean$T
      
      ELBO_val <- 
        n_j * (log(par_p) - log(par_q)) - 
        (par_p - par_q) * E_T
      
    } else if (par_name == "sigma"){
      par_p <- 
        (info$mean$T * info$mean$Q2) %>% 
        # add over sample (j's)
        colSums
      par_q <- lambda$sigma
      E_sigma <- info$mean$sigma
      
      ELBO_val <- 
        (par_q - par_p) * E_sigma
    } else if (par_name == "Q"){
      par_p <- 
        matrix(info$mean$T, ncol = 1) %*% 
        info$mean$sigma
      par_q <- lambda$Q2
      E_Q2 <- info$mean$Q2
      
      ELBO_val <- 
        (par_q - par_p) * E_Q2
    }
    ELBO_val
  }