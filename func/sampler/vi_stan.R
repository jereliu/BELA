# call glrm_sampler_stan
library(rstan)
library(parallel)

glrm_sampler_vi_stan <- 
  function(Y, lambda, 
           family_name = c("gaussian", "poisson"), 
           prior_name = c("gaussian", "dirichlet"),
           init, config, rec, info)
  {
    glrm_sampler_stan(
      Y, lambda, 
      family_name, prior_name,
      init, config, rec, info,
      stan_algorithm = "vi"
    )
  }

# glrm_sampler_vi_stan <- 
#   function(Y, lambda, 
#            family_name = c("gaussian", "poisson"), 
#            prior_name = c("gaussian", "dirichlet"),
#            init, config, rec, info)
#   {
#     family_name <- match.arg(family_name)
#     prior_name <- match.arg(prior_name)
#     
#     # unpack family properties
#     n <- info$n 
#     p <- info$p
#     k <- info$k
#     true_theta <- info$true_par$theta
#     
#     # unpack mcmc parameters
#     record_freq <- config$record_freq
#     time_max <- config$time_max
#     
#     iter_max <- config$sampler$iter_max
#     step_size <- config$sampler$step_size
#     frog_step <- config$sampler$frog_step  
#     mmtm_freq <- config$sampler$mmtm_freq 
#     rotn_freq <- config$sampler$rotn_freq  
#     
#     U_cur <- init$U
#     V_cur <- init$V
#     Ru_cur <- matrix(rnorm(n*k), ncol = k)
#     Rv_cur <- matrix(rnorm(p*k), ncol = k) 
#     
#     # initiate sampler
#     rstan_options(auto_write = TRUE)
#     options(mc.cores = parallel::detectCores())    
#     stan_addr <- "./func/sampler/stan/"
#     model_name <- family_name
#     
#     if (prior_name == "dirichlet") {# TODO
#       model_name <- prior_name
#     }
#     
#     stan_file <- paste0(stan_addr, model_name, ".stan")
#     
#     stan_data <- list(N = n, P = p, K = k, Y = Y, 
#                       lambda_u = lambda, 
#                       lambda_v = lambda)
#     init_func <- function() init
#     
#     # execute sampler
#     cat("\n Compiling..")
#     time0 <- Sys.time()
#     model_out <- 
#       stan_model(stan_file, model_name) %>%
#       vb(data = stan_data, init = init_func,
#          seed = 100, 
#          output_samples = iter_max)
#     # time0 <- as.POSIXlt(model_out@date, format = "%c")
#     time_max <- difftime(Sys.time(), time0, units = "mins")
#     cat("Done")
#     
#     # return
#     rec <- rstan::extract(model_out)
#     rec$Theta <- 
#       lapply(1:dim(rec$U)[1], 
#              function(i) rec$U[i, ,] %*% t(rec$V[i, , ])) %>% abind(along = 0)
#     
#     rec$obj <- rec$lp__
#     rec$time <- seq(0, time_max, length.out = iter_max)
#     
#     # return
#     rec
#   }