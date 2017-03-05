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
n_rep_total <- length(data_id)
n_rep <- 1000

cfig_check_list <- 
  expand.grid(
    FAMILY = unique(cfig_list$FAMILY), 
    K = unique(cfig_list$K), 
    SNR = unique(cfig_list$SNR), 
    SAMPLR = unique(cfig_list$SAMPLR),
    REP = 1:n_rep)

n_rep <- length(unique(cfig_check_list$REP))

#### 1. raw trace ####
mmd_stat <- TRUE
ks_stat <- FALSE

array_list_full <- list.files(tar_dir)
array_list <- array_list[grep("geweke.RData", array_list_full)]

# plot, raw trace
sec_max <- 350
plot_idx <- 2:(sec_max+2)
FAMILY <- unique(cfig_list$FAMILY)
K <- unique(cfig_list$K)
trfunc <- function(x, m, v) (x - m)/v

for (family_id in 1:length(FAMILY)){
  family <- c("gaussian", "poisson")[family_id]
  for (k_id in 1:length(K)){
    k <- K[k_id]
    pdf(paste0(tar_addr, "trace_", family, "_k", k, ".pdf"), width = 9, height = 6)
    snr <- 100
    plot_name <-
      paste0(family, "_k", k, "_snr", snr)
    sec <- 0:sec_max
    file_handle <-
      paste0(family, "_k", k, "_snr", snr, "_", "hmc_stan")
    
    plot(sec, quant_line[1, ], 
         type = "n", ylab = "2nd Largest eig(Theta)", 
         ylim = c(0, 0.7),
         main = file_handle)
    
    {    
      # hmc
      load(paste0(tar_dir, 
                  grep(file_handle, array_list, value = TRUE)))
      
      plot_name <-
        paste0(family, "_k", k, "_snr", snr)
      
      mean_line <- colMeans(theta_container[, plot_idx]) 
      quant_line <- 
        apply(theta_container[, plot_idx], 2, 
              function(x) quantile(x, c(0.025, 0.975)))
      m <- mean(quant_line[1, ])
      v <- max(quant_line)
      
      mean_line <- trfunc(mean_line, m, v)
      quant_line <- trfunc(quant_line, m, v) 
        
      lines(mean_line, lwd = 2)
      polygon(c(sec,rev(sec)),
              c(quant_line[1, ],rev(quant_line[2, ])),
              col = rgb(0,0,0,0.5), border = NA)
      
      # gibbs
      file_handle <-
        paste0(family, "_k", k, "_snr", snr, "_", "gibbs")
      load(paste0(tar_dir, 
                  grep(file_handle, array_list, value = TRUE)))
      
      # preprocessing
      mean_line <- colMeans(theta_container[, plot_idx])
      quant_line <- 
        apply(theta_container[, plot_idx], 2, 
              function(x) quantile(x, c(0.025, 0.975, 0.95)))
      mean_line <- trfunc(mean_line, m, v)
      quant_line <- trfunc(quant_line, m, v) 
      
      # plot
      lines(mean_line, col = 2, lwd = 2)
      polygon(c(sec,rev(sec)),
              c(quant_line[1, ], rev(quant_line[2, ])),
              col = rgb(1,0,0,0.5), border = NA)
    }
    abline(h = mean(quant_line[3, sec_max]), lty = 2)
    legend("topright", lty = 1, lwd = 2, col = 1:2, legend = c("HMC", "Gibbs"))
    dev.off()
  }
}

#### 2. raw trace ####
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
      abline(h = 0)
      abline(h = mixing_tvdist_list[[1]][1, 2], lty = 2)
      legend("topright", lty = 1, lwd = 2, col = 1:2, 
             legend = c("Gibbs", "HMC"))
    }
  }
}

