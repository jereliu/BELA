objective <- 
  function(par, prior, info,
           type = c("original", "primal", "dual")){
    link_func <- 
      paste0("link_", prior$Q$link) %>% 
      parse(text = .) %>% eval
    
    f_Q <- link_func(par$Q)
    sigma_T_Q <- (par$T %*% t(par$sigma)) * f_Q
    sigma_T_Q[sigma_T_Q == 0] <- .Machine$double.eps
    
    loss <- 
      sigma_T_Q - info$stat$n_ij * log(sigma_T_Q)
    
    if (is.null(prior$Q$Sigma)){
      # objective function with 1 parameter ADMM
      if (type == "original"){
        penalty <- 
          prior$lambda$Q * nuclear(par$Q)
      } else if (type == "primal"){
        penalty <- 
          c(prior$lambda$Q * nuclear(par$Z),
            prior$rho * sum((par$Z - par$Q)^2)/2)
      } else if (type == "dual"){
        penalty <- 
          c(prior$lambda$Q * nuclear(par$Z),
            sum(par$U * (par$Q - par$Z)),
            prior$rho * sum((par$Z - par$Q)^2)/2
          )
      }
    } else {
      # objective function for 4 parameter ADMM
      if (type == "original"){
        penalty <- 
          prior$lambda$Q * nuclear(par$Q)
      } else if (type == "primal"){
        penalty <- 
          c(prior$lambda$Q * nuclear(par$Z),
            prior$rho * 
              sum((par$Z - prior$Q$Gamma$X %*% par$W)^2)/2,
            prior$rho * 
              sum((par$W - par$S %*% t(prior$Q$Gamma$Y))^2)/2,
            prior$rho * sum((par$S - par$Q)^2)/2
          )
      } else if (type == "dual"){
        penalty <- 
          c(prior$lambda$Q * nuclear(par$Z),
            prior$rho * 
              sum((par$Z - prior$Q$Gamma$X %*% par$W)^2)/2,
            prior$rho * 
              sum((par$W - par$S %*% t(prior$Q$Gamma$Y))^2)/2,
            prior$rho * sum((par$S - par$Q)^2)/2,
            sum(par$Uz * (par$Z - prior$Q$Gamma$X %*% par$W)),
            sum(par$Uw * (par$W - par$S %*% t(prior$Q$Gamma$Y))),
            sum(par$Us * (par$S - par$Q))
          )
      }
    }
    
    list(total = sum(loss) + sum(penalty), 
         scaler = c(sum(loss), penalty), 
         loss_mat = loss)
  }


pred_dist_abs <- 
  function(par, prior, info){
    f_Q <- par$Q * (par$Q > 0)
    pred <- (par$T %*% t(par$sigma)) * f_Q
    info$stat$n_ij - pred
  }

rv_coef <- function(X, Y){
  # center column vector
  X <- t(t(X) - colMeans(X))
  Y <- t(t(Y) - colMeans(Y))
  
  S_xy <- t(X) %*% Y
  S_xx <- t(X) %*% (X)
  S_yy <- t(Y) %*% (Y) 
  
  COVV <- Trace(S_xy %*% t(S_xy))
  VAV_X <- Trace(S_xx %*% S_xx) 
  VAV_Y <- Trace(S_yy %*% S_yy) 
  
  COVV/sqrt(VAV_X * VAV_Y)
}

#### extra functions ####
nuclear <- function(X) sum(svd(X)$d)
