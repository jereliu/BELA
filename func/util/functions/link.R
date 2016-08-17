# 
link_identity <- function(x) x
link_identity_dir <- function(x) 
  matrix(1, nrow = nrow(x), ncol = ncol(x))

#
link_soft <- function(x, k = 10){
  x / (1 - exp(-k*x))
}

link_soft_dir <- function(x, k = 10){
  (1 - exp(-k*x))^(-1) -
    k * x /((1 - exp(-k*x)) * (exp(k*x)-1))
}

#
link_pos <- function(x) x*I(x>0)
link_pos_dir <- function(x) I(x>0)
  
# 
link_pos2 <- function(x) x^2 *I(x>0)
link_pos2_dir <- function(x) 2 * x * I(x>0)

# 
link_exp <- function(x) exp(x)
link_exp_dir <- function(x) exp(x)
