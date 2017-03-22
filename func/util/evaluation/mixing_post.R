library(KSD)

default <- FALSE
if (default){
  n_iter <- dim(rec$Theta)[1]
  idx_range <- round(n_iter/2):n_iter
  ty <-
    mixing_stein_svd(
      eig_list = rec$eig_list[idx_range, ],
      # parameters to produce score function
      Theta = rec$Theta[idx_range, , ],
      U = rec$U[idx_range, , ],
      V = rec$V[idx_range, , ],
      T_suff = rec$Y,
      lambda = rec$info$lambda,
      dist_family = rec$info$family,
      n_eigen = 1
    )
}

# KSD on leading ksd values of matrix
mixing_stein_svd <- 
  function(
    eig_list, score_list = NULL, n_eigen = NULL, 
    # parameters to produce score function
    Theta, U, V, T_suff, lambda, family_name,
    # hyperparameter for Gaussian kernel
    width = -1, nboot = 1000, ...){
    # define global variable
    n_sample <- dim(Theta)[1]
    if (is.null(n_eigen)) 
      # if no input then examine all singular values
      n_eigen <- ncol(eig_list)
    
    if (is.null(score_list)){
      # produce score value evaluated at singular values
      score_list <- 
        score_svd_batch(
          n_eigen, lambda, family_name,
          Theta, U, V, T_suff
        )
    }
    
    if (length(score_list) == n_sample)
      score_list <- t(score_list)
    
    # produce singular value for KSD()
    eig_list <- eig_list[, 1:n_eigen]
    score_list <- score_list[, 1:n_eigen]
    
    # evaluate KSD
    KSD(eig_list, score_list, 
        width = width, nboot = nboot)
  }


score_svd_batch <- 
  function(n_eigen, 
           lambda, family_name,
           Theta, U, V, T_suff){
    n_sample <- dim(Theta)[1]
    dist_family <- glrm_family(family_name)
    
    score_list <- 
      sapply(1:n_sample, 
             function(i)
               score_svd(
                 Theta[i, , ], 
                 U = U[i, , ], V = V[i, , ],
                 T_suff, lambda, dist_family, 
                 n_eigen = n_eigen)
      ) %>% t
    
    if(nrow(score_list) == 1)
      score_list <- t(score_list)
    
    score_list
  }

# KSD on distance between U and V parameters
mixing_stein_param <-
  function(
    negrad_func, theta_list, theta_row = NULL, 
    # hyperparameter for Gaussian kernel
    width = -1, nboot = 1000, ...)
  {
    # gaussian kernel with bandwidth equal to 1/alpha
    neggrad_list <- 
      apply(theta_list, 1, 
            function(theta){
              negrad_func(matrix(theta, nrow = theta_row), ...)
            }
      ) %>% t 
    
    KSD(theta_list, -neggrad_list, 
        width = width, nboot = nboot)
  }

