library(numDeriv)
library(VGAM) # for Nakagami distribution

# DOCUMENTATION
# functions to calculate derivative (wrt VARIATIONAL PARAMETER!)
# of sample-generation functions (aka inverse CDF)

# INPUT
# mode:
# method: method used to calculate derivative, either 
#         "numeric": numeric diff. using numDeriv package
#         "implicit" (TODO) close-form deriv using implicit function theorem

#### 2. Derivative Functions ####

genfunc_deriv_T <- 
  function(var_par, info, prior,
           mode = "marginal", method = "numeric")
  {
    if (mode != "marginal"){
      stop("only 'marginal' supported for T")
    }
    
    if (method == "numeric"){
      deriv <- 
        sapply(1:J, 
               function(j) 
                 grad(genfunc_T, 
                      var_par[j], 
                      u = info$sample$u[j], 
                      n_j = info$stat$n_j[j])
        )
    } else if (method == "implicit") {
      stop("method 'implicit' not yet implemented")
      # ## TO CORRECT
      # deriv <- 
      #   info$sample$T *
      #   (pgamma(info$sample$T,
      #           shape = info$stat$n_j + 1,
      #           rate = var_par) -
      #      pgamma(info$sample$T,
      #             shape = info$stat$n_j,
      #             rate = var_par)
      #   )/dgamma(info$sample$T,
      #            shape = info$stat$n_j,
      #            rate = var_par)
    }
    
    deriv
  }



genfunc_deriv_sigma <- 
  function(var_par, info, prior,
           mode = "marginal", method = "numeric")
  {
    if (mode != "marginal"){
      stop("only 'marginal' supported for sigma")
    }
    
    if (method == "numeric"){
      deriv <- 
        sapply(1:I, 
               function(i) 
                 grad(genfunc_sigma, 
                      var_par[i], 
                      u = info$sample$u[i], 
                      n_i = info$stat$n_i[i], 
                      a = prior$sigma$a,
                      b = prior$sigma$b)
        )
    } else if (method == "implicit") {
      stop("method 'implicit' not yet implemented")
      # ## TO CREATE/CORRECT
    }
    deriv
  }

genfunc_deriv_Q <- 
  function(var_par, info, prior,
           mode = "marginal", method = "numeric"){
    if (mode == "marginal"){
      if (method == "numeric"){
        deriv <- 
          sapply(1:J, 
                 function(j) 
                   grad(
                     genfunc_Q, 
                     var_par[, j], 
                     u = info$sample$u[, j])
          )
      } else if (method == "implicit") {
        stop("method 'implicit' not yet implemented")
        # ## TO CREATE/CORRECT
      }    
    } else if (mode == "copula") {
      warning("mode 'copula' not implemented")
      deriv <- 0
    }
    deriv
  }

#### 2. Generation Function ####
genfunc_T <- 
  function(lambda, u, n_j) 
    qgamma(u, n_j, rate = lambda)

genfunc_sigma <-
  function(lambda, u, n_i, a, b) 
    qbeta(u, a + n_i + lambda, b)

genfunc_Q <- # nakagami distribution
  function(lambda, u){
    qnaka(u, scale = 1.5 * lambda, shape = 1.5) 
  }