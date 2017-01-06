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
  K = cfig$K
  SNR = cfig$SNR
  LAMBDA = cfig$LAMBDA
  FAMILY = cfig$FAMILY
  SAMPLR = cfig$SAMPLR
}

glrm_worker <- 
  function(
    data_seed, 
    trial_seed,
    n = 10, p = 50, 
    K = c(2, 15, 20), 
    SNR = 100,
    LAMBDA = 2,
    FAMILY = c("gaussian", "poisson"),
    SAMPLR = c("gibbs", "hmc_stan", "vi_stan"), 
    record_freq = 10
  ){
    rand_seeds <- 
      list(data = data_seed, 
           samplr = trial_seed)
    
    save_rec_target <- TRUE
    
    for (family_name in FAMILY){
      for (snr in SNR){
        for (k in K){
          for (lambda in LAMBDA){
            phi_sd = 1/sqrt(lambda)
            step_optim = 0.001
            step_sampl = 0.01
            
            set.seed(rand_seeds$data)
            data.sim <- 
              glrm_sample(
                n = n, # num of samples 
                p = p, # num of categories
                k = k, # factor dimension
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
              # choose iter based on 
              iter_max <- c(1e5, 5e3)
              
              if (TRUE){
                # if (is.null(rec$init)){
                init_MAP <- FALSE
                init <- NULL
                # init$V <- t(data.sim$V)
                # init$U <- t(data.sim$U)
              } else {
                init_MAP <- FALSE
                init <- rec$init
                # init$V <- rec$V[dim(rec$U)[1], , ]
                # init$U <- rec$U[dim(rec$U)[1], , ]
              }
              
              #### 2. Sampling ####
              set.seed(rand_seeds$samplr)
              rec <- 
                glrm(Y, lambda = lambda, k = k, 
                     true_par = true_par,
                     init = init, init_MAP = init_MAP,
                     samplr_name = samplr_name,
                     family_name = family_name,
                     iter_max = iter_max, 
                     record_freq = record_freq,
                     time_max = 60, 
                     step_size = c(step_optim, step_sampl), 
                     frog_step = 5,
                     rotn_freq = iter_max[2],
                     mmtm_freq = iter_max[2], 
                     samp_seed = rand_seeds$samplr)
              
              rec$data.sim <- data.sim
              
              #### 3. Evaluation ####
              # plot
              task_title <- 
                paste0(
                  family_name, " ", 
                  samplr_name, " k=", k, " snr = ", snr)
              
              # save (disabled for cluster)
              var_name <- paste0(family_name,
                                 "_k", k, "_snr", snr,
                                 "_", samplr_name)
              paste0("rec_", var_name, " <- rec") %>%
                parse(text = .) %>% eval()
              
              if (save_rec_target){
              rec_target <- 
                matrix(c(data.sim$theta[1, 1], rec$Theta[, 1, 1]), nrow = 1) 
              }
              
              # save file 
              if (save_rec_target){
                file_name <- 
                  paste0("./result/mixing_stat/rec_", 
                         var_name,  "_", data_seed, ".csv")
                
                if (!file.exists(file_name)){
                  colnames(rec_target) <-
                    c("S0", 1:(length(rec_target)-1))
                }
                
                write.table(
                  rec_target,
                  append = file.exists(file_name),
                  quote = FALSE, sep = ",",
                  qmethod = "double",
                  file = file_name,
                  row.names = FALSE,
                  col.names = !file.exists(file_name))
                
                
                # save(rec_target, 
                #      file = 
                #        paste0("./result/mixing_stat/rec_", 
                #               var_name,  "_", data_seed, "_", 
                #               trial_seed, ".RData"))
              } else {
                rec2 <- rec
                rec <- NULL
                rec_idx <- 
                  seq(1, dim(rec2$U)[1], 
                      length.out = iter_max[2]/record_freq) %>% round
                  
                rec$U <- rec2$U[rec_idx, , ]
                rec$V <- rec2$V[rec_idx, , ]
                rec$Y <- rec2$Y
                  
                var_name <- 
                  paste0(family_name, "_k", k, "_snr", snr, "_", samplr_name)
                
                save(rec, file = 
                       paste0("./result/rec_", var_name, "_", data_seed, "_", trial_seed, ".RData"))
              }
              
              # # save test statistic
              # mixing_stat <- 
              #   sapply(
              #     1:length(rec$U[, 1, 1]),
              #     function(iter)
              #       mean(rec$data.sim$U[1, 1] < rec$U[1:iter, 1, 1])
              #   )
              # write.csv(mixing_stat, 
              #           file = 
              #             paste0("./result/mixing_stat/rec_", 
              #                    var_name, "_mixtest_", data_seed, ".csv"), 
              #           row.names = FALSE)
            }
          }
        }
      }
    }
    
    return(0)
  }