rm(list = ls()) #WoG

library(magrittr)
library(dplyr)
library(ggplot2)
library(abind)
library(gplots)
source("./func/util/source_Dir.R")
sourceDir("./func")

#### 1. Data Generation ####
n <- 20
p <- 100
k <- 20
family_name <- c("gaussian", "poisson")[1]
samplr_name <- "hmc"
snr <- 10


for (family_name in c("gaussian", "poisson")){
  #for (snr in c(100, 10, 1, 0.5)){
  for (k in c(2, 5, 10, 15, 20)){
    phi_sd = 0.5
    step_optim = 0.001
    step_sampl = 0.01

    set.seed(100)
    data.sim <- 
      glrm_sample(
        n = n, # num of samples 
        p = p, # num of categories
        k = k, # factor dimension
        family_name = family_name,
        phi_sd = phi_sd,
        noise_sd = phi_sd/snr
      )
    
    Y <- data.sim$Y
    true_theta <- data.sim$theta
    init <- NULL
    init$V <- t(data.sim$V)
    
    for (samplr_name in c("gibbs", "hmc")){
      if (samplr_name == "gibbs") {
        iter_max <- c(1e5, 1e3)
      } else if (samplr_name == "hmc"){
        iter_max <- c(1e5, 2e5)
      }
      #### 2. Sampling ####
      set.seed(1000)
      rec <- 
        glrm(Y, lambda = 1, 
             k = k, true_theta = true_theta,
             init = init,
             samplr_name = samplr_name,
             family_name = family_name,
             iter_max = iter_max, 
             record_freq = 10,
             time_max = 60, 
             step_size = c(step_optim, step_sampl), 
             frog_step = 5,
             rotn_freq = iter_max[2],
             mmtm_freq = iter_max[2])
      
      #### 3. Evaluation ####
      # plot
      task_title <- 
        paste0(
          family_name, " ", 
          samplr_name, " k=", k, " snr = ", snr)
      plot(rec$time, rec$error, type = "l", 
           main = paste0("sample error ", task_title))
      plot(rec$pred_error, type = "l", 
           main = paste0("prediction error ", task_title))
      plot(rec$time, rec$obj, type = "l", 
           main = paste0("obj ", task_title))
      plot(rec$time[-1], rec$esdj, type = "l", 
           main = paste0("obj ", task_title))
      
      # save 
      var_name <- paste0("rec_", family_name, 
                         "_k", k, "_snr", snr,
                         "_", samplr_name)
      paste0(var_name, " <- rec") %>% 
        parse(text = .) %>% eval()
      save(rec, file = paste0("./result/", var_name, ".RData"))
    }
  }
}
#}