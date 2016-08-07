link_soft <- function(x, k = 10){
  x / (1 - exp(-k*x))
}

link_soft_dir <- function(x, k = 10){
  (1 - exp(-k*x))^(-1) -
    k * x /((1 - exp(-k*x)) * (exp(k*x)-1))
}

link_pos <- function(x) x*I(x>0)
