library(magrittr)
library(dplyr)
library(ggplot2)
library(abind)
library(gplots)

source("./func/util/source_Dir.R")
sourceDir("./func")

default <- FALSE
if (default){
  data_seed = cfig$data_seed
  trial_seed = cfig$trial_seed
  n = N
  p = P
  K_true = cfig$K_true
  K_model = cfig$K_model 
  SNR = cfig$SNR
  LAMBDA = cfig$LAMBDA
  FAMILY = cfig$FAMILY
  PRIOR = cfig$PRIOR
  SAMPLR = cfig$SAMPLR
  iter_max = c(1e5, 1e2)
  parm_updt = c("U", "V")
}

glrm_worker <- 
  function(
    data_seed, 
    trial_seed,
    n = 50, p = 100, 
    K_true = c(2, 15, 20), 
    K_model = c(2, 15, 20), 
    SNR = 100,
    LAMBDA = 2,
    FAMILY = c("gaussian", "poisson"),
    PRIOR = c("gaussian", "sparse"),
    SAMPLR = c("gibbs", "hmc_stan", "vi_stan"), 
    record_freq = 10,
    parm_updt = c("U", "V"),
    iter_max = c(1e5, 5e3)
  ){
    rand_seeds <- 
      list(data = data_seed, samplr = trial_seed)
    
    for (family_name in FAMILY){
      for (prior_name in PRIOR){
        for (snr in SNR){
          for (k_true in K_true){
            for (k_model in K_model){
              for (lambda in LAMBDA){
                
                phi_sd = 1/sqrt(lambda)
                step_optim = 0.001
                step_sampl = 0.01
                
                set.seed(rand_seeds$data)
                data.sim <-
                  glrm_sample(
                    n = n, # num of samples 
                    p = p, # num of categories
                    k = k_true, # factor dimension
                    family_name = family_name,
                    lambda = lambda,
                    noise_sd = 1/(sqrt(lambda)*snr)
                  )
                
                Y <- data.sim$Y
                true_par <- data.sim[c("U", "V", "theta")]
                true_par$sd <- lambda
                #true_par$theta <- theta_true
                
                rec <- NULL
                #rec$init <- init_hmc
                for (samplr_name in SAMPLR){
                  
                  if (TRUE){
                    # if (is.null(rec$init)){
                    init_MAP <- FALSE
                    init <- NULL
                    init$V <- t(data.sim$V)
                    init$U <- t(data.sim$U)
                  } else {
                    init_MAP <- FALSE
                    init <- rec$init
                    # init$V <- rec$V[dim(rec$U)[1], , ]
                    # init$U <- rec$U[dim(rec$U)[1], , ]
                  }
                  
                  #### 2. Sampling ####
                  set.seed(rand_seeds$samplr)
                  rec <-
                    glrm(Y, lambda = lambda, k = k_model, 
                         true_par = true_par,
                         init = init, init_MAP = init_MAP,
                         samplr_name = samplr_name,
                         family_name = family_name,
                         prior_name = prior_name,
                         iter_max = iter_max, 
                         record_freq = record_freq,
                         time_max = 60, 
                         step_size = c(step_optim, step_sampl), 
                         frog_step = 5,
                         rotn_freq = iter_max[2],
                         mmtm_freq = iter_max[2], 
                         samp_seed = rand_seeds$samplr, 
                         parm_updt = parm_updt
                    )
                  
                  rec$data.sim <- data.sim
                  
                  #### 3. Evaluation ####
                  # plot
                  task_title <- 
                    paste0(
                      family_name, " ", samplr_name, 
                      " k_true=", k_true, " k_model=", k_model, 
                      " snr = ", snr)
                  
                  # save (disabled for cluster)
                  save_eigen <- TRUE
                  save_score <- TRUE
                  
                  var_name <- 
                    paste0(family_name, "_", prior_name,
                           "_ktr", k_true, "_kmd", k_model,
                           "_snr", snr, "_", samplr_name)
                  paste0("rec_", var_name, " <- rec") %>%
                    parse(text = .) %>% eval()
                  
                  file_name <- 
                    paste0("./result/mixing_stat/rec_", 
                           var_name, "_", data_seed, "_", 
                           trial_seed, ".csv")
                  
                  if (save_eigen){
                    # save Theta 1. 1
                    # rec_target <- 
                    #   matrix(c(data.sim$theta[1, 1], 
                    #            rec$Theta[, 1, 1]), nrow = 1) 
                    
                    n_eigen <- 5
                    # save 2nd largest eigen value
                    eig0 <- svd(data.sim$theta)$d[1:n_eigen]
                    eig_list <- 
                      apply(rec$Theta, 1, 
                            function(theta) 
                              svd(theta)$d[1:n_eigen])
                    rec_target <- cbind(eig0, eig_list)
                    rec$eigen <- rec_target
                    
                    if (!file.exists(file_name)){
                      colnames(rec_target) <-
                        c("S0", 1:(ncol(rec_target)-1))
                    }
                    
                    rec_target <- 
                      rbind( c(0, 0, rec$time), rec_target)
                    
                    write.table(
                      rec_target,
                      append = FALSE, #file.exists(file_name),
                      quote = FALSE, sep = ",",
                      qmethod = "double",
                      file = file_name,
                      row.names = FALSE,
                      col.names = !file.exists(file_name)
                    )
                  } 
                  if (save_score) {
                    # compute likelihood score information for ksd
                    # add to existing file
                    score_list <-
                      score_svd_batch(
                        n_eigen = 5, 
                        lambda, family_name,
                        rec$Theta, rec$U, rec$V, 
                        T_suff = rec$Y
                      ) 
                    
                    write.table(
                      cbind(0, t(score_list)/(n*p*k_model)), # add 0 in place of initial values
                      append = TRUE, #file.exists(file_name),
                      quote = FALSE, sep = ",",
                      qmethod = "double",
                      file = file_name,
                      row.names = FALSE,
                      col.names = !file.exists(file_name)
                    )
                    # save(rec, 
                    #      file = 
                    #        paste0("./result/mixing_stat/rec_", 
                    #               var_name, "_", data_seed, "_", 
                    #               trial_seed, ".RData")
                    # )
                  }
                }
              }
            }
          }
        }
      }
    }
    
    return(0)
  }



plot(0,0, type = "n",
     xlim = c(1,iter_max[2]), ylim = c(0,0.1))
for (i in 1:k_model){
  lines(rec$gamma_rank_cumprod[,,i], col = i)
}



sapply(1:iter_max[2],
       function(i)
         sum(rec$gamma_rank_cumprod[i, , ] > 1e-2)
) %>% plot(type = "l", ylim = c(0, k_model))
