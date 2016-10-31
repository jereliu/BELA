# acceptance probability for Gibbs sampler
acc_prob_U <- 
  function(U_cur, U_old, V_cur, V_old, 
           lambda, family){
    dim1 <- nrow(U_cur)
    dim2 <- ncol(U_cur)
    
    negloglik <- family$negloglik
    prob_ij <- function(i, j){
      U_prop <- U_old
      U_prop[i, j] <- U_cur[i, j]
      
      # probability
      (-negloglik(T, U_prop %*% t(V_cur)) + 
        negloglik(T, U_old %*% t(V_cur)) -
        lambda * 
        (sum(U_prop^2) - sum(U_old^2))
      ) %>% min(0, .) %>% exp
    }
    
    acc_prob <- 
      outer(1:dim1, 1:dim2, Vectorize(prob_ij))
    acc_prob
  }

acc_prob_V <- 
  function(U_cur, U_old, V_cur, V_old, 
           lambda, family){
    dim1 <- nrow(V_cur)
    dim2 <- ncol(V_cur)
    
    negloglik <- family$negloglik
    prob_ij <- function(i, j){
      V_prop <- V_old
      V_prop[i, j] <- V_cur[i, j]
      
      # probability
      (-negloglik(T, U_cur %*% t(V_prop)) + 
        negloglik(T, U_cur %*% t(V_old)) -
        lambda * 
        (sum(V_prop^2) - sum(V_old^2))
      ) %>% min(0, .) %>% exp
    }
    
    acc_prob <- 
      outer(1:dim1, 1:dim2, Vectorize(prob_ij))
    acc_prob
  }

# acceptance probability for rotation sampler
acc_prob_rotation <- 
  function(U_cur, U_old, V_cur, V_old, 
           lambda, family){
    
  }

