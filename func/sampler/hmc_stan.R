# call glrm_sampler_stan
library(rstan)
library(parallel)

glrm_sampler_hmc_stan <- 
  function(Y, lambda, 
           family_name = c("gaussian", "poisson"), 
           prior_name = c("gaussian", "dirichlet"),
           init, config, rec, info)
  {
    glrm_sampler_stan(
      Y, lambda, 
      family_name, prior_name,
      init, config, rec, info,
      stan_algorithm = "hmc"
    )
  }