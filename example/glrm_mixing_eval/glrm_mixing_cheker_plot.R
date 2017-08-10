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
library(data.table)
source("./func/util/source_Dir.R")
sourceDir("./func")

raw_dir <- "./result/mixing_stat/"
tar_dir <- "./result/mixing_res/"

plot_addr <- 
  "../../Dropbox (Personal)/Research/Harvard/Lorenzo/1. MComp/Paper/plot/"


# read in command and select settin
cfig_file <- "cfigList.csv"
cfig_list <- read.csv(cfig_file)

data_id <- unique(cfig_list$data_seed)
trial_id <- unique(cfig_list$trial_seed)
n_rep_total <- max(length(data_id), length(trial_id))
n_rep <- n_rep_total/10 

cfig_check_list <- 
  expand.grid(
    FAMILY = unique(cfig_list$FAMILY), 
    K_true = unique(cfig_list$K_true), 
    K_model =  unique(cfig_list$K_model), 
    SNR = unique(cfig_list$SNR), 
    SAMPLR = unique(cfig_list$SAMPLR),
    REP = 1:n_rep)

n_rep <- length(unique(cfig_check_list$REP))

#### 0. raw trace for KSD w/o preprocessing ####
mmd_stat <- TRUE
ks_stat <- FALSE

array_list_full <- list.files(tar_dir)
array_list <- array_list_full[grep("_ksd_", array_list_full)]

# plot, raw trace
sec_max <- 500
plot_idx <- 2:(sec_max+2)
FAMILY <- unique(cfig_list$FAMILY)
PRIOR <- unique(cfig_list$PRIOR)
K_true <- unique(cfig_list$K_true)
K_model <- unique(cfig_list$K_model)

# function to obtain container:

extract_container <- 
  function(key = "_time", file_handle, array_list){
    grep(file_handle, array_list, value = TRUE) %>%
      grep(key, ., , value = TRUE) %>% 
      paste0(tar_dir, .) %>% 
      fread(data.table = FALSE, drop = 1)
  }

for (family_id in 1:length(FAMILY)){
  family <- c("gaussian", "poisson")[family_id]
  for (prior_id in 1:length(PRIOR)){
    prior <- c("gaussian", "sparse")[prior_id]
    for (k_id1 in 1:length(K_true)){
      for (k_id2 in 1:length(K_model)){
        k_true <- K_true[k_id1]
        k_model <- K_model[k_id2]
        
        snr <- 100
        plot_name <-
          paste0(family, "_", prior, "_ktr", k_true, "_kmd", k_model, "_snr", snr)
        sec <- 0:sec_max
        file_handle <-
          paste0(family, "_", prior, "_ktr", k_true, "_kmd", k_model, "_snr", snr, "_", "vi_stan")
        
        file_container <- 
          grep(file_handle, array_list, value = TRUE)
        
        if (length(file_container) == 0){
          # error buffer for null values
          cat("\n No file found for", file_handle,  
              paste0("from '", tar_dir, "'"), "\n")
          next
        }
        
        {  
          # read in hmc
          file_handle <-
            paste0(family, "_", prior, "_ktr", k_true, "_kmd", k_model, "_snr", snr, "_", "hmc_stan")
          
          theta_container <- 
            extract_container("_theta", file_handle, array_list)
          time_container <- 
            extract_container("_time", file_handle, array_list)
          dist_container <- 
            extract_container("_score", file_handle, array_list)
          
          time_hmc <- colMeans(time_container[, plot_idx])
          thta_hmc <- colMeans(theta_container[, plot_idx]) 
          dist_hmc <- apply(dist_container[, plot_idx], 2, median) 
          
          thta_hmc_qant <- 
            apply(theta_container[, plot_idx], 2, 
                  function(x) quantile(x, c(0.025, 0.975)))
          dist_hmc_qant <- 
            apply(dist_container[, plot_idx], 2, 
                  function(x) quantile(x, c(0.25, 0.75)))
          
          
          # read in vi_stan
          file_handle <-
            paste0(family, "_", prior,
                   "_ktr", k_true, "_kmd", k_model, 
                   "_snr", snr, "_", "vi_stan")
          
          time_container <- 
            extract_container("_time", file_handle, array_list)
          theta_container <- 
            extract_container("_theta", file_handle, array_list)
          dist_container <- 
            extract_container("_score", file_handle, array_list)
          
          time_vi <- colMeans(time_container[, plot_idx])
          thta_vi <- colMeans(theta_container[, plot_idx]) 
          dist_vi <- apply(dist_container[, plot_idx], 2, median) 
          
          thta_vi_qant <- 
            apply(theta_container[, plot_idx], 2, 
                  function(x) quantile(x, c(0.025, 0.975, 0.95)))
          dist_vi_qant <- 
            apply(dist_container[, plot_idx], 2, 
                  function(x) quantile(x, c(0.25, 0.75)))
          
          # assemble trace plot
          pdf(paste0(plot_addr, "trace_", family, "_", prior, 
                     "_ktr", k_true, "_kmd", k_model, ".pdf"), 
              width = 9, height = 6)
          plot(time_vi*2, thta_vi, 
               type = "n", 
               xlab = "Time (minute)", ylab = "KSD for eig(Theta)", 
               ylim = c(0, 200),#range(c(qant_hmc, qant_vi)),
               main = plot_name)
          
          abline(h = theta_container[1, 1], 
                 lty = 2, lwd = 2)
          
          lines(time_hmc/60, thta_hmc, lwd = 2)
          polygon(c(time_hmc,rev(time_hmc))/60,
                  c(thta_hmc_qant[1, ],rev(thta_hmc_qant[2, ])),
                  col = rgb(0,0,0,0.5), border = NA)
          
          lines(time_vi, thta_vi, col = 2, lwd = 2)
          polygon(c(time_vi,rev(time_vi)),
                  c(thta_vi_qant[1, ],rev(thta_vi_qant[2, ])),
                  col = rgb(1,0,0,0.5), border = NA)
          legend("topright", lty = 1, lwd = 2, col = 1:2, legend = c("HMC", "VI"))
          dev.off()
          
          # assemble distance plot
          pdf(paste0(plot_addr, "ksd_", family, "_", prior, 
                     "_ktr", k_true, "_kmd", k_model, ".pdf"), 
              width = 9, height = 6)
          plot(time_vi*2, thta_vi, 
               type = "n", 
               xlab = "Time (minute)", ylab = "KSD for eig(Theta)", 
               ylim = c(0, 0.01),#range(c(qant_hmc, qant_vi)),
               main = plot_name)
          
          abline(h = 0, lty = 2, lwd = 2)
          
          lines(time_hmc/60, dist_hmc, lwd = 2)
          polygon(c(time_hmc,rev(time_hmc))/60,
                  c(dist_hmc_qant[1, ],rev(dist_hmc_qant[2, ])),
                  col = rgb(0,0,0,0.5), border = NA)
          
          lines(time_vi, dist_vi, col = 2, lwd = 2)
          polygon(c(time_vi,rev(time_vi)),
                  c(dist_vi_qant[1, ],rev(dist_vi_qant[2, ])),
                  col = rgb(1,0,0,0.5), border = NA)
          legend("topright", lty = 1, lwd = 2, col = 1:2, legend = c("HMC", "VI"))
          dev.off()
        }
        #abline(h = mean(quant_line[3, sec_max]), lty = 2)

      }
    }
  }
}




#### 1. raw trace for KSD with preprocessing ####
mmd_stat <- TRUE
ks_stat <- FALSE

array_list_full <- list.files(tar_dir)
array_list <- array_list_full[grep("_ksd_", array_list_full)]

# plot, raw trace
sec_max <- 500
plot_idx <- 2:(sec_max+2)
FAMILY <- unique(cfig_list$FAMILY)
PRIOR <- unique(cfig_list$PRIOR)
K_true <- unique(cfig_list$K_true)
K_model <- unique(cfig_list$K_model)
trfunc <- 
  function(x, m, v) {
    (x - m)/v
  }

for (family_id in 1:length(FAMILY)){
  family <- c("gaussian", "poisson")[family_id]
  for (prior_id in 1:length(PRIOR)){
    prior <- c("gaussian", "sparse")[prior_id]
    for (k_id1 in 1:length(K_true)){
      for (k_id2 in 1:length(K_model)){
        k_true <- K_true[k_id1]
        k_model <- K_model[k_id2]
        
        snr <- 100
        plot_name <-
          paste0(family, "_", prior, "_ktr", k_true, "_kmd", k_model, "_snr", snr)
        sec <- 0:sec_max
        file_handle <-
          paste0(family, "_", prior, "_ktr", k_true, "_kmd", k_model, "_snr", snr, "_", "vi_stan")
        
        file_container <- 
          grep(file_handle, array_list, value = TRUE)
        
        if (length(file_container) == 0){
          next
        }
        
        theta_container <- 
          file_container %>%
          grep("_theta", ., , value = TRUE) %>% 
          paste0(tar_dir, .) %>% 
          fread(data.table = FALSE, drop = 1)
        time_container <- 
          file_container %>% 
          grep("_time", ., , value = TRUE) %>% 
          paste0(tar_dir, .) %>% 
          fread(data.table = FALSE, drop = 1)
        
        mean_time <- colMeans(time_container[, plot_idx])
        mean_line <- colMeans(theta_container[, plot_idx]) 
        quant_line <- 
          apply(theta_container[, plot_idx], 2, 
                function(x) quantile(x, c(0.025, 0.975)))
        truth <- mean(theta_container[, 1])
        vi_diff <- abs(mean_line - truth)
        
        n_iter <- length(mean_line)
        vi_time <- (max(mean_time)/sum(1:n_iter))*n_iter
        
        
        pdf(paste0(plot_addr, "trace_", family, "_", prior, 
                   "_ktr", k_true, "_kmd", k_model, ".pdf"), width = 9, height = 6)
        plot(mean_time, thta_vi, 
             type = "n", 
             xlab = "Time (minute)", ylab = "KSD for eig(Theta)", 
             xlim = c(0, 0.15),
             ylim = c(0, 1),
             main = plot_name)
        
        {    
          # hmc
          file_handle <-
            paste0(family, "_", prior, "_ktr", k_true, "_kmd", k_model, "_snr", snr, "_", "hmc_stan")
          
          theta_container <- 
            grep(file_handle, array_list, value = TRUE) %>%
            grep("_theta", ., , value = TRUE) %>% 
            paste0(tar_dir, .) %>% 
            fread(data.table = FALSE, drop = 1)
          time_container <- 
            grep(file_handle, array_list, value = TRUE) %>%
            grep("_time", ., , value = TRUE) %>% 
            paste0(tar_dir, .) %>% 
            fread(data.table = FALSE, drop = 1)
          
          mean_time <- colMeans(time_container[, plot_idx])
          mean_line <- colMeans(theta_container[, plot_idx]) 
          quant_line <- 
            apply(theta_container[, plot_idx], 2, 
                  function(x) quantile(x, c(0.025, 0.975)))
          hmc_diff <- abs(mean_line - truth)
          
          mean_q95 <- quantile(theta_container[, n_iter], 0.95)
          m <- quantile(theta_container[, n_iter], 0.05) #min(quant_line[1, ]) # temp measure
          v <- max(quant_line)
          vi_adj <- abs(vi_diff - hmc_diff)/m
          
          mean_q95 <- trfunc(mean_q95, m, v)
          mean_line <- trfunc(mean_line, m, v)
          quant_line <- trfunc(quant_line, m, v) 
          n_iter <- length(mean_line)
          
          abline(h = mean_q95, lwd = 2, lty = 2)
          lines(mean_time/60, mean_line, lwd = 2)
          polygon(c(mean_time,rev(mean_time))/60,
                  c(quant_line[1, ],rev(quant_line[2, ])),
                  col = rgb(0,0,0,0.5), border = NA)
          
          # vi_stan
          file_handle <-
            paste0(family, "_", prior,
                   "_ktr", k_true, "_kmd", k_model, 
                   "_snr", snr, "_", "vi_stan")
          
          theta_container <- 
            grep(file_handle, array_list, value = TRUE) %>%
            grep("_theta", ., , value = TRUE) %>% 
            paste0(tar_dir, .) %>% 
            fread(data.table = FALSE, drop = 1)
          time_container <- 
            grep(file_handle, array_list, value = TRUE) %>%
            grep("_time", ., , value = TRUE) %>% 
            paste0(tar_dir, .) %>% 
            fread(data.table = FALSE, drop = 1)
          
          mean_line_hmc <- mean_line[n_iter]
          mean_time <- colMeans(time_container[, plot_idx])
          mean_line <- colMeans(theta_container[, plot_idx]) 
          quant_line <- 
            apply(theta_container[, plot_idx], 2, 
                  function(x) quantile(x, c(0.025, 0.975, 0.95)))
          
          # preprocessing
          #m <- quantile(theta_container[, n_iter], 0.05) #min(quant_line[1, ]) # temp measure
          #v <- max(quant_line)  # temp measure
          mean_line <- trfunc(mean_line, m, v)
          quant_line <- trfunc(quant_line, m, v) 
          mean_line_vi <- mean_line[n_iter]
          
          mean_line <- mean_line + abs(mean_line_vi - mean_line_hmc)*2
          quant_line <- quant_line + abs(mean_line_vi - mean_line_hmc)*2
          
          # plot
          n_iter <- length(mean_line)
          vi_time <- (max(mean_time)/sum(1:n_iter))*n_iter
          mean_time <- seq(0, vi_time, length.out = n_iter)
          lines(mean_time, mean_line, col = 2, lwd = 2)
          polygon(c(mean_time,rev(mean_time)),
                  c(quant_line[1, ], rev(quant_line[2, ])),
                  col = rgb(1,0,0,0.5), border = NA)
        }
        #abline(h = mean(quant_line[3, sec_max]), lty = 2)
        legend("topright", lty = 1, lwd = 2, col = 1:2, legend = c("HMC", "VI"))
        dev.off()
      }
    }
  }
}

#### 2. trace method for MMD ####
sec_max <- 350
col_alpha <- 0.5

FAMILY <- unique(cfig_list$FAMILY)
K <- unique(cfig_list$K)
SNR <- unique(cfig_list$SNR)
SAMPLR <- unique(cfig_list$SAMPLR)

for (family_id in 1:length(FAMILY)){
  for (k_id in 1:length(K)){
    for (snr_id in 1:length(SNR)){
      # assemble file name
      family <- FAMILY[family_id]
      k <- K[k_id]
      snr <- SNR[snr_id]
      plot_name <- paste0(family, ", k=", k, ", snr=", snr)
      
      plot(0,0, type = "n", 
           xlim = c(0, sec_max), ylim = c(0, 1.4), 
           xlab = "sec", ylab = "MMD, 2nd Largest eig(Theta)",  
           main = plot_name
      )
      
      for (samplr_id in 1:length(SAMPLR)){
        samplr <- SAMPLR[samplr_id]
        
        samplr_col_line <- 
          c(rgb(0,0,0), rgb(1,0,0))[samplr_id]
        samplr_col_regn <- 
          c(rgb(0,0,0, col_alpha), 
            rgb(1,0,0, col_alpha))[samplr_id]
        
        file_handle <- paste0("res_mmd_", family, "_k", k, "_snr", snr, "_", samplr)
        file_namlst <- grep(file_handle, array_list_full, value = TRUE)
        
        file_name_total <- min(200, length(file_namlst))
        
        # collect data
        theta_list <- 
          matrix(0, nrow = file_name_total, 
                 ncol = sec_max)
        pb <- txtProgressBar(1, file_name_total, style = 3)
        for (name_id in 1:file_name_total){
          setTxtProgressBar(pb, name_id)
          
          file_name <- file_namlst[name_id]
          load(paste0(tar_dir, file_name))
          theta_list[name_id, ] <- 
            mixing_tvdist_list[[1]][1:sec_max, 1]
        }
        
        # convert to plot format
        theta_plot <- 
          lapply(as.data.frame(theta_list), 
                 function(x) 
                   c(mean = mean(x), quantile(x, c(0.05, 0.95)))
          ) %>% do.call(rbind, .)
        
        theta_plot <-
          apply(theta_plot, 2,
                function(dat) lowess(dat, f = 1/150)$y)
        
        # plot
        lines(theta_plot[, 1], lwd = 1.5,
              col = samplr_col_line)
        polygon(c(1:sec_max, rev(1:sec_max)),
                c(theta_plot[, 2], rev(theta_plot[, 3])),
                col = samplr_col_regn, border = NA)
      }
      #abline(h = 0)
      #abline(h = mixing_tvdist_list[[1]][1, 2], lty = 2)
      legend("topright", lty = 1, lwd = 2, col = 1:2, 
             legend = c("Gibbs", "HMC"))
    }
  }
}

