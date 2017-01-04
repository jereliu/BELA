rm(list = ls()) #WoG

library(magrittr)
library(dplyr)
library(ggplot2)
library(abind)
library(gplots)
library(data.table)
library(psych)
source("./func/util/source_Dir.R")
sourceDir("./func")

raw_dir <- "./result/mixing_stat/"
tar_dir <- "./result/mixing_model/"

#### 0.specify model categories ####
FAMILY <- c("gaussian", "poisson")[1]
K <- c(2, 10, 15)[1] #c(2, 10, 15)
SNR <- 100
SAMPLR <- c("gibbs", "hmc_stan")
#c("gibbs", "hmc_stan", "vi_stan")

#### 1. restructure raw file into array #### 
file_list <- list.files(raw_dir)
cfig_list <- read.csv("cfigList.csv")
data_id <- unique(cfig_list$data_seed)
trial_id <- unique(cfig_list$trial_seed)


for (family in FAMILY){
  for (k in K){
    for (snr in SNR){
      for (sampler in SAMPLR){
        # capture relevant file index
        sett_handle <-
          paste0(family, "_k", k, "_snr", snr, "_", sampler)
        data_list <- 
          grep(sett_handle, file_list, value = TRUE)
        #initiate container for all quantiles 
        quat_container <- NULL
        
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
          data_file <- 
            read.csv(paste0(raw_dir, data_file_names))
          
          if (is.null(quat_container)){
            quat_container <- 
              matrix(NaN, 
                     nrow = length(data_list), 
                     ncol = ncol(data_file)-1) 
          }
          
          #### 1.2. compute quantile then store ----
          # compute quantile estimated from data
          quantile_data <- 
            apply(data_file[, -1], 2, 
                  function(x) mean(x<data_file[, 1]))
          
          # add quantile computed from data to 
          # container for current setting
          quat_container[d_id, ] <- quantile_data
          
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
        save(quat_container, 
             file = paste0(tar_dir, sett_handle, ".RData"))
      }
    }
  }
}

#### 2. analyze array using ecdf #### 
array_list <- list.files(tar_dir)
mixing_tvdist_list <- list()

for (family in FAMILY){
  for (k in K){
    for (snr in SNR){
      for (sampler in SAMPLR){
        # capture relevant file index
        file_handle <-
          paste0(family, "_k", k, "_snr", snr, "_", sampler)
        # load quat_container, row for data, col for iteration
        load(paste0(tar_dir, file_handle, ".RData"))
        
        komo_stat <-
          apply(quat_container, 2, prior_mixing_cook)
        komo_conf <- 
          qksone(0.95, dim(quat_container)[2])$root 
        
        mixing_tvdist_list[[file_handle]]$stat <- komo_stat
        mixing_tvdist_list[[file_handle]]$uppr <- komo_conf
        
        print(paste0(file_handle, " Done!"))
      }
    }
  }
}

save(mixing_tvdist_list, 
     file = paste0(tar_dir, "mixing_tvdist_list.RData")
)



file_handle <- paste0(family, "_k", k, "_snr", snr, "_gibbs")
obs_idx <- seq(1, 2e3, 10)

plot(obs_idx, 
     mixing_tvdist_list[[file_handle]]$stat[obs_idx], 
     type = "l", ylim = c(0, 0.2), 
     xlab = "iteration", ylab = "Kolmogorov–Smirnov Test Statistic", 
     main = paste0(family, ", k=", k, ", snr=", snr))
abline(h = mixing_tvdist_list[[file_handle]]$uppr, 
       lty = 2, col = 2)

file_handle <- paste0(family, "_k", k, "_snr", snr, "_hmc_stan")
lines(obs_idx, 
      mixing_tvdist_list[[file_handle]]$stat[obs_idx], 
      col = 4)

legend("topright", lty = 1, 
       col = c(1, 4), legend = c("Gibbs", "HMC"))
