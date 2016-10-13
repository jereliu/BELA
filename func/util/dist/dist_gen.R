library(gamlss)

binom_gen <- function(p) rbinom(1, 1, p)
poiss_gen <- function(lambda) rpois(1, lambda)
gauss_gen <- function(mu) rnorm(1, mu, sd = 1)
logno_gen <- function(mu) exp(rnorm(1, mu, 1))
