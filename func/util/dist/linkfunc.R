softplus <- function(theta){
  log(1 + exp(theta))
}

reluappr <- function(theta){
  theta/(1 - exp(-100 * theta))
}
