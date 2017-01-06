rm(list = ls()) #WoG

library(magrittr)
library(dplyr)
library(ggplot2)
library(abind)
library(gplots)
library(data.table)
library(psych)
library(infotheo)
library(kernlab)
source("./func/util/source_Dir.R")
sourceDir("./func")

raw_dir <- "./result/mixing_stat/"
tar_dir <- "./result/mixing_model/"

#### 0.specify model categories ####
FAMILY <- c("gaussian", "poisson")[1]
K <- c(2, 10, 15) #c(2, 10, 15)
SNR <- 100
SAMPLR <- c("gibbs", "hmc_stan")
#c("gibbs", "hmc_stan", "vi_stan")

#### 1. restructure raw file into array #### 
file_list <- list.files(raw_dir)
cfig_list <- read.csv("cfigList.csv")
data_id <- unique(cfig_list$data_seed)
trial_id <- unique(cfig_list$trial_seed)

for (family_name in FAMILY){
  for (k in K[2]){
    for (snr in SNR){
      for (sampler in SAMPLR){
        # capture relevant file index
        sett_handle <-
          paste0(family_name, "_k", k, "_snr", snr, "_", sampler)
        data_list <- 
          grep(sett_handle, file_list, value = TRUE)
        #initiate container for all quantiles 
        theta_container <- NULL
        
        cat(paste0("\n ", sett_handle, " in process...\n"))
        pb <- 
          txtProgressBar(min = 1, max = length(data_list), style = 3)
        for (d_id in 1:length(data_list)){
          setTxtProgressBar(pb, d_id)
          # for each dataset in setting, read-in all reps
          data_file_names <- data_list[d_id]
          
          #### 1.1 harvest sample to compute quantile ----
          # read in initial file, find chain length, build file container
          # also build quantile container for current setting if not exist
          metric_KL <- FALSE
          metric_stein <- FALSE
          if (metric_stein){
            # store Y, U_1 and V_1
            file_name <- paste0(raw_dir, data_file_names)
            load(file_name)
            if (is.null(theta_container)){
              iter_max <- dim(rec$U)[1] #obtain number of iteration
              iter_idx <- 1:50
              # iter_idx <- 
              #   c(seq(1, round(iter_max/2), length.out = 24),
              #     seq(round(iter_max/2)+100, iter_max, length.out = 16)) %>%
              #   round %>% unique
              theta_container$Y <- rec$Y
              theta_container$iter_idx <- iter_idx
              theta_container$theta <- 
                array(NaN, 
                      dim = 
                        c(length(iter_idx), # iteration
                          length(data_list), # repetition
                          dim(rec$U)[3]*
                            sum(dim(rec$U)[2] + dim(rec$V)[2])) # number of variable
                )
            }
          } else if (metric_KL) {
            # store Y, U_1 and V_1
            file_name <- paste0(raw_dir, data_file_names)
            load(file_name)
            if (is.null(theta_container)){
              theta_container <- 
                array(NaN, 
                      dim = 
                        c(length(data_list), # repetition
                          dim(rec$U)[1], # iteration
                          dim(rec$U)[3]*2 + 1) # number of variable
                )
            }
          } else {
            data_file <- 
              read.csv(paste0(raw_dir, data_file_names))
            
            if (is.null(theta_container)){
              theta_container <- 
                matrix(NaN, 
                       nrow = length(data_list), 
                       ncol = ncol(data_file))
            }
          }
          #### 1.2. compute obtain quantities then store ----
          
          # add quantile computed from data to 
          # container for current setting
          if (metric_stein) {
            theta_container$theta[, d_id, ] <- 
              sapply(iter_idx, 
                     function(id)
                       as.vector(rbind(rec$U[id, , ], 
                                       rec$V[id, , ]))
              ) %>% t
          } else if (metric_KL){
            theta_container[d_id, , ] <- 
              cbind(
                # sufficient statistics
                rep(glrm_family(family_name)$sufficient(
                  rec$data.sim$Y[1, 1]), dim(rec$U)[1]), 
                rec$U[, 1, ], rec$V[, 1, ])
          } else {
            theta_container[d_id, ] <- 
              as.matrix(data_file)[1, ]
          }
          # # visualize the chain
          # plot(file_container[1, ], 
          #      col = rgb(0,0,0,0.2),
          #      type = "l")
          # for (i in 1:nrow(file_container)){
          #   lines(file_container[i, ], 
          #         col = rgb(0,0,0,0.3))
          # }
        }
        
        # 
        if (metric_stein){
          postfix <- "_stein.RData"
        } else if (metric_KL){
          postfix <- "_KL.RData"
        } else {
          postfix <- "_geweke.RData"
        }
        
        save(theta_container,
             file = paste0(tar_dir, sett_handle, postfix))
      }
    }
  }
}

#### 2. analyze array using ecdf ####
calc_metric <- TRUE

if (calc_metric){
  array_list <- list.files(tar_dir)
  array_list <- array_list[grep("geweke.RData", array_list)]
  mixing_tvdist_list <- list()
  
  # stein test p-value
  for (family_name in FAMILY){
    lambda <- 2
    family_obj <- glrm_family(family_name)
    loglikfunc <- function(u, v, T_suff){
      family_obj$negloglik(sum(u*v), T_suff) +
        lambda/2 * (sum(u^2) + sum(v^2))
    }
    
    for (k in K[2]){
      for (snr in SNR){
        for (sampler in SAMPLR){
          # capture relevant file index
          file_handle <-
            paste0(family_name, "_k", k, "_snr", snr, "_", sampler)
          # load theta_container, row for data, col for iteration
          print(paste0("loading ", file_handle, ".."))
          
          load(paste0(tar_dir, 
                      grep(file_handle, array_list, value = TRUE)))
          print("Done! Calculating Test Statistic..")
          
          metric_stein <- FALSE 
          metric_KL <- FALSE

          if (metric_stein){
            ksd_list <- rep(0, dim(theta_container$theta)[1])
            p_list <- rep(0, dim(theta_container$theta)[1])
            pb <-
              txtProgressBar(
                min = 1,
                max = dim(theta_container$theta)[1],
                style = 3)
            
            for (i in 1:dim(theta_container$theta)[1]){
              setTxtProgressBar(pb, i)
              S_cur <- theta_container$theta[i, , ]
              T_suff <- theta_container$Y
              
              result <- 
                mixing_stein(
                  grad_neglik_S,
                  S_cur, 
                  theta_row = sum(dim(T_suff)),
                  width = -1,
                  nboot = 1000,
                  T_suff = T_suff,
                  lambda = lambda,
                  dist_family = glrm_family(family_name)
                )
              ksd_list[i] <- result$info$ksd_V
              p_list[i] <- result$p
            }
            print(p_list)
            mixing_tvdist_list[[file_handle]] <- ksd_list
            
          } else if (metric_KL) {
            # if collects whole uv_1, then check KL-divergence
            KL_div <-
              rep(NaN, length = dim(theta_container)[2])
            
            pb <-
              txtProgressBar(
                min = 1, max = dim(theta_container)[2],
                style = 3)
            
            for (i in 1:dim(theta_container)[2]){
              setTxtProgressBar(pb, i)
              KL_div[i] <-
                mixing_KL_dist(
                  theta_container[, i, ],
                  loglikfunc
                )
            }
            
            mixing_tvdist_list[[file_handle]] <- KL_div
            
          } else {
            iter_idx <- 1:200
              # seq(1, (ncol(theta_container)-1), 10)
            komo_stat <-
              matrix(NaN,
                nrow = length(iter_idx), ncol = 3)
            pb <-
              txtProgressBar(
                min = 1,
                max = ncol(theta_container)-1,
                style = 3)
            
            for (i in iter_idx){
              setTxtProgressBar(pb, i)
              mmdo <- 
                kmmd(matrix(theta_container[, i+1]),
                           matrix(theta_container[, 1]))
              komo_stat[i, ] <- 
                c(mmdo@mmdstats[1], mmdo@Radbound, mmdo@Asymbound)
                
                # prior_mixing_geweke_perm(
                #   theta_container[, i+1],
                #   prior = theta_container[, 1])
            }
            # komo_conf <-
            #   qksone(0.95, dim(theta_container)[1])$root
            
            mixing_tvdist_list[[file_handle]] <- komo_stat
          }
          
          #mixing_tvdist_list[[file_handle]]$uppr <- komo_conf
          
          print(paste0(file_handle, " Done!"))
        }
      }
    }
  }
  
  save(mixing_tvdist_list,
       file = paste0(tar_dir, "mixing_tvdist_geweke_pval.RData")
  )
  
  # ks-stat
  array_list <- list.files(tar_dir)
  array_list <- array_list[grep("snr", array_list)]
  mixing_tvdist_list <- list()
  
  for (family_name in FAMILY){
    for (k in K[2]){
      for (snr in SNR){
        for (sampler in SAMPLR){
          # capture relevant file index
          file_handle <-
            paste0(family_name, "_k", k, "_snr", snr, "_", sampler)
          # load theta_container, row for data, col for iteration
          print(paste0("loading ", file_handle, ".."))
          
          load(paste0(tar_dir, file_handle, ".RData"))
          print("Done! Calculating Test Statistic..")
          
          komo_stat <-
            rep(NaN, length = ncol(theta_container)-1)
          
          pb <-
            txtProgressBar(
              min = 1, max = ncol(theta_container)-1, style = 3)
          for (i in 1:(ncol(theta_container)-1)){
            setTxtProgressBar(pb, i)
            komo_stat[i] <-
              prior_mixing_geweke_stat(
                theta_container[, i+1],
                prior = theta_container[, 1])
          }
          # komo_conf <-
          #   qksone(0.95, dim(theta_container)[1])$root
          
          mixing_tvdist_list[[file_handle]] <- komo_stat
          #mixing_tvdist_list[[file_handle]]$uppr <- komo_conf
          
          print(paste0(file_handle, " Done!"))
        }
      }
    }
  }
  
  save(mixing_tvdist_list,
       file = paste0(tar_dir, "mixing_tvdist_geweke_stat.RData")
  )
}

plot(mixing_tvdist_list[[1]][, 1], type = "l", 
     ylab = "MMD Statistic for Theta[1,1]", xlab = "Iteration")
abline(h = mixing_tvdist_list[[1]][, 2], lty = 2)
lines(mixing_tvdist_list[[2]][, 1], col = 2)

# #### 3. plot ecdf #### 
# type <- c("pval", "stat")[2]
# load(paste0(tar_dir, "mixing_tvdist_geweke_", type, ".RData"))
# 
# quant_crit <- qksone(0.95, 1000) 
# 
# par(mfrow = c(3, 1))
# for (family_name in FAMILY){
#   for (k in K){
#     for (snr in SNR){
#       for (sampler in SAMPLR){
#         # capture relevant file index
#         file_handle <-
#           paste0(family_name, "_k", k, "_snr", snr, "_", sampler)
#         # load theta_container, row for data, col for iteration
#         setting <- 
#           grep(file_handle, names(mixing_tvdist_list), value = TRUE)
#         
#         if (sampler == "gibbs"){
#           plot(mixing_tvdist_list[[setting]], type = "l", 
#                main = file_handle, 
#                xlab = "Iteration", 
#                ylab = "Kolmogorov-Smirnov Statistic", 
#                ylim = c(0, 0.65)
#           )
#           abline(h = quant_crit$root)
#         } else {
#           lines(mixing_tvdist_list[[setting]], col = 2)
#         }
#       }
#     }
#   }
# }
# 
# 
# prior_mixing_jere(theta_container)

# iter_idx <- 
#   c(seq(1, round(iter_max/2), length.out = 24),
#     seq(round(iter_max/2)+100, iter_max, length.out = 16)) %>%
#   round %>% unique
# plot(1:50, mixing_tvdist_list[[2]])
