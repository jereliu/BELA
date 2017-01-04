# # acceptance probability for Gibbs sampler
# acc_prob_U_ind <- 
#   function(U_cur, U_old, V_cur, V_old, 
#            lambda, family, T_suff){
#     dim1 <- nrow(U_cur)
#     dim2 <- ncol(U_cur)
#     
#     negloglik <- family$negloglik
#     prob_ij <- function(i, j){
#       U_prop <- U_old
#       U_prop[i, j] <- U_cur[i, j]
#       
#       # probability
#       (+(-negloglik(T_suff, U_prop %*% t(V_cur)) + 
#            negloglik(T_suff, U_old %*% t(V_cur)) -
#            (lambda/2) * 
#            (sum(U_prop^2) - sum(U_old^2))
#       )) %>% min(0, .) %>% exp
#     }
#     
#     acc_prob <- 
#       outer(1:dim1, 1:dim2, Vectorize(prob_ij))
#     acc_prob
#   }
# 
# acc_prob_V_ind <- 
#   function(U_cur, U_old, V_cur, V_old, 
#            lambda, family, T_suff){
#     dim1 <- nrow(V_cur)
#     dim2 <- ncol(V_cur)
#     
#     negloglik <- family$negloglik
#     prob_ij <- function(i, j){
#       V_prop <- V_old
#       V_prop[i, j] <- V_cur[i, j]
#       
#       # probability
#       (+(-negloglik(T_suff, U_cur %*% t(V_prop)) + 
#            negloglik(T_suff, U_cur %*% t(V_old)) -
#            (lambda/2) * 
#            (sum(V_prop^2) - sum(V_old^2))
#       )) %>% min(0, .) %>% exp
#     }
#     
#     acc_prob <- 
#       outer(1:dim1, 1:dim2, Vectorize(prob_ij))
#     acc_prob
#   }

# acceptance probability for Gibbs sampler
acc_prob_U <- 
  function(U_cur, U_old, V_cur, V_old, i, 
           lambda, family, T_suff){
    dim1 <- nrow(U_cur)
    dim2 <- ncol(U_cur)
    
    negloglik <- family$negloglik
    d1 <- family$partition$d
    d2 <- family$partition$d2
    
    #### calculate proposal components
    # new
    Theta_new <- U_cur %*% t(V_old)
    A_d1_new <- d1(Theta_new)
    A_d2_new <- d2(Theta_new)
    B_new <- T_suff - A_d1_new + A_d2_new * Theta_new
    lnr_coef_u_new <- B_new %*% V_old    # n x k, each row lnr coef for u_i
    sigma_new <-
      solve(
        t(V_old) %*% diag(A_d2_new[i, ]) %*% V_old +
          lambda * diag(dim2)
      )
    mu_new <-  sigma_new %*% lnr_coef_u_new[i, ]    
    
    # old
    Theta_old <- U_old %*% t(V_old)
    A_d1_old <- d1(Theta_old)
    A_d2_old <- d2(Theta_old)
    B_old <- T_suff - A_d1_old + A_d2_old * Theta_old
    lnr_coef_u_old <- B_old %*% V_old    # n x k, each row lnr coef for u_i
    sigma_old <-
      solve(
        t(V_old) %*% diag(A_d2_old[i, ]) %*% V_old +
          lambda * diag(dim2)
      )
    mu_old <- sigma_old %*% lnr_coef_u_old[i, ]
    
    #### calculate proposal distribution
    p_ratio <- 
      ((-negloglik(T_suff, U_cur %*% t(V_cur)) + 
          negloglik(T_suff, U_old %*% t(V_cur)) -
          (lambda/2) * 
          (sum(U_cur^2) - sum(U_old^2))
      ))
    
    q_ratio <-
      dmvnorm(U_cur[i, ], mean = mu_old, 
              sigma = sigma_old, log = TRUE) - 
      dmvnorm(U_old[i, ], mean = mu_new, 
              sigma = sigma_new, log = TRUE)
    
    acc_prob <- p_ratio - q_ratio
    acc_prob <- acc_prob %>% min(0, .) %>% exp
    
    acc_prob
  }

acc_prob_V <- 
  function(U_cur, U_old, V_cur, V_old, j,
           lambda, family, T_suff){
    dim1 <- nrow(V_cur)
    dim2 <- ncol(V_cur)
    
    negloglik <- family$negloglik
    d1 <- family$partition$d
    d2 <- family$partition$d2
    
    #### calculate proposal components
    # new
    Theta_new <- U_cur %*% t(V_cur)
    A_d1_new <- d1(Theta_new)
    A_d2_new <- d2(Theta_new)
    B_new <- T_suff - A_d1_new + A_d2_new * Theta_new
    lnr_coef_v_new <- t(B_new) %*% U_cur # p x k, each row lnr coef for v_j
    sigma_new <- solve(
      t(U_cur) %*% diag(A_d2_new[, j]) %*% U_cur +
        lambda * diag(dim2))
    mu_new <- sigma_new %*% lnr_coef_v_new[j, ]
    
    # old
    Theta_old <- U_cur %*% t(V_old)
    A_d1_old <- d1(Theta_old)
    A_d2_old <- d2(Theta_old)
    B_old <- T_suff - A_d1_old + A_d2_old * Theta_old
    lnr_coef_v_old <- t(B_old) %*% U_cur # p x k, each row lnr coef for v_j
    sigma_old <- solve(
      t(U_cur) %*% diag(A_d2_old[, j]) %*% U_cur +
        lambda * diag(dim2))
    mu_old <- sigma_old %*% lnr_coef_v_old[j, ]
    
    #### calculate proposal distribution
    p_ratio <- 
      ((-negloglik(T_suff, U_cur %*% t(V_cur)) + 
          negloglik(T_suff, U_cur %*% t(V_old)) -
          (lambda/2) * 
          (sum(V_cur^2) - sum(V_old^2))
      ))
    
    
    q_ratio <-
      dmvnorm(V_cur[j, ], mean = mu_old, 
              sigma = sigma_old, log = TRUE) - 
      dmvnorm(V_old[j, ], mean = mu_new, 
              sigma = sigma_new, log = TRUE)
    
    acc_prob <- p_ratio - q_ratio
    acc_prob <- acc_prob %>% min(0, .) %>% exp
    
    acc_prob
  }


# acceptance probability for rotation sampler
acc_prob_rotation <- 
  function(U_cur, U_old, V_cur, V_old, 
           lambda, family){
    
  }

