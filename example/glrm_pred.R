rm(list = ls()) #WoG

library(magrittr)
library(dplyr)
library(ggplot2)
library(abind)
library(gplots)
source("./func/util/source_Dir.R")
sourceDir("./func")

#### 1. Data Generation ####
n <- 10
p <- 100
k <- 2
family_name <- c("gaussian", "poisson", "poisson_softplus")[2]
samplr_name <- "stein"
snr <- 100
edge_max <- 2
record_freq <- 10

rand_seeds <- list(data = 7200, samplr = 4200)
rec_plot <- FALSE
marg_eig_plot <- TRUE
cond_dens_plot_d1_U <- FALSE
cond_dens_plot_d1_V <- FALSE
marg_dens_plot <- FALSE
marg_dens_plot_slice <- FALSE

# addr_targ <-
#  paste0("../../Dropbox (Personal)/Research/Harvard/Lorenzo/1. BayesOpt/Report/Progress/2017_Nov_Week_4/plot/")
# pdf(paste0(addr_targ, "debug_marginal_varV_poisson.pdf"), height= 4, width = 10)
# par(mfrow = c(1, 2))
for (family_name in c("gaussian", "poisson")[2]){
  #for (snr in c(100, 10, 1, 0.5)){
  for (k in c(1, 2, 5, 10, 15, 20)[2]){
    #for (lambda in c(0.5, 1, 3, 5, 10, 20)[1]){
    lambda = 10
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
    for (samplr_name in c("gibbs", "hmc_stan", "vi_stan", "slice", "stein")[c(1:2)]){
      # choose iter based on 
      if (length(grep("gibbs|slice", samplr_name)) > 0){
        iter_max <- c(1e5, 1e4) # 1e3)
      } else {
        # if sampler name contain "hmc"...
        iter_max <- c(1e5, 1e4) # 1e4)
      }
      
      if (TRUE){
        # if (is.null(rec$init)){
        init_MAP <- FALSE
        init <- NULL
        # init$U <- t(data.sim$U)
        init$V <- t(data.sim$V)
        parm_updt <- c("U", "V")[1]
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
             # slice parameters
             edge_max = edge_max,
             # hmc parameters
             frog_step = 5,
             rotn_freq = iter_max[2],
             mmtm_freq = iter_max[2],
             parm_updt = parm_updt,
             # stein parameters
             n_particle = 1000
        )
      
      rec$data.sim <- data.sim
      
      if (samplr_name == "hmc_stan"){
        rec_hmc <- rec
      } else if (samplr_name == "gibbs"){
        rec_gibbs <- rec
      }
      
      #### 3. Evaluation ####
      # plot
      task_title <- 
        paste0(
          family_name, " ", 
          samplr_name, " k=", k, " snr = ", snr, " ",
          paste(parm_updt, collapse = "_")
          )
      
      if (rec_plot){
        #plot(rec$time, rec$error, type = "l", 
        #     main = paste0("sample error ", task_title))
        time_size <- length(rec$time)
        plot(rec$pred_error, type = "l",
             main = paste0("prediction error ", task_title))
        
        # objective function
        plot(rec$time[(time_size-length(rec$obj)+1):time_size], 
             rec$obj, 
             type = "l", xlab = "Time (min)",
             main = paste0("obj ", task_title))
        T_y <- glrm_family(family_name)$sufficient(Y)
        abline(h = 
                 glrm_family(family_name)$negloglik(T_y, (rec$init$U) %*% t(rec$init$V)) + 
                 (1/phi_sd) * (sum(rec$init$U^2) + sum(rec$init$V^2))
        )
      }
      
      if (cond_dens_plot_d1_U){ # U
        if (k !=1){
          stop("conditional density plot work only for k=1")
        }
        # exact posterior density matching for U[1, 1]
        u_likhd <-
          function(u_range = seq(-1,1, 1e-2), data.sim, lambda, 
                   p, i = 1, family_name){
            # exact likelihood for k=1
            if (family_name == "gaussian"){
              var <- 1/(sum(data.sim$V^2) + lambda)
              mu <- var * data.sim$V %*% data.sim$Y[i, ]
              density_out <- dnorm(u_range, mu, sd = sqrt(var))
            } else {
              family <- glrm_family(family_name)
              T_suff <- family$sufficient(data.sim$Y)
              VT <- data.sim$V %*% data.sim$Y[i, ]
              lik_val <- 
                sapply(u_range, 
                       function(u){
                         u * VT - sum(exp(u * data.sim$V)) -
                           0.5 * lambda * u^2
                       }
                ) %>% exp 
              density_out <- lik_val/max(lik_val)
            }
            density_out
          }
        
        emprden_hmc <- density(rec$U[, 1, 1])
        u_range <-
          seq(min(emprden_hmc$x), max(emprden_hmc$x),
              length.out = 1e3)
        postden_hmc <-
          u_likhd(u_range = u_range,
                  data.sim, lambda = lambda, p, 
                  family_name = family_name)
        
        plot(emprden_hmc$x, emprden_hmc$y/max(emprden_hmc$y),
             type = "l", #xlim = c(-1, 0),
             main = paste0(task_title, " lambda = ", lambda))
        lines(u_range, postden_hmc/max(postden_hmc),
              col = 2)
      }
      
      if (cond_dens_plot_d1_V){ # V
        if (k !=1){
          stop("conditional density plot work only for k=1")
        }
        # exact posterior density matching for V[1, 1]
        v_likhd <-
          function(v_range = seq(-1,1, 1e-2), data.sim, lambda, 
                   p, i = 1, family_name){
            # exact likelihood for k=1
            if (family_name == "gaussian"){
              var <- 1/(sum(data.sim$U^2) + lambda)
              mu <- var * data.sim$U %*% data.sim$Y[, i]
              density_out <- dnorm(v_range, mu, sd = sqrt(var))
            } else {
              family <- glrm_family(family_name)
              T_suff <- family$sufficient(data.sim$Y)
              UT <- data.sim$U %*% T_suff[, i]
              lik_val <- 
                sapply(v_range, 
                       function(v){
                         v * UT - sum(exp(v * data.sim$U)) -
                           0.5 * lambda * v^2
                       }
                ) %>% exp 
              density_out <- lik_val/max(lik_val)
            }
            density_out
          }
        
        
        emprden_hmc <- density(rec$V[, 1, 1])
        v_range <-
          seq(min(emprden_hmc$x), max(emprden_hmc$x),
              length.out = 1e3)
        postden_hmc <-
          v_likhd(v_range = v_range,
                  data.sim, lambda = lambda, p, 
                  family_name = family_name)
        
        plot(emprden_hmc$x, emprden_hmc$y/max(emprden_hmc$y),
             type = "l", #xlim = c(-3, 3),
             main = paste0(task_title, " lambda = ", lambda)
        )
        lines(v_range, postden_hmc/max(postden_hmc),
              col = 2)
      }
      
      if (marg_dens_plot){
        # exact posterior density matching for U[1, 1]
        # u_likhd <- 
        #   function(u_range = seq(-3,3, 1e-2), data.sim, lambda, p, i = 1){
        #     # exact likelihood for k=1 normal data
        #     y2 <- sum(data.sim$Y^2)
        #     
        #     exp(- 0.5 * (y2 + lambda * u_range^2 + lambda^2) * u_range^2 / (u_range^2 + lambda)) * 
        #       (u_range^2 + lambda)^(-p/2)
        #   }
        
        emprden_hmc <- density(rec$U[, 1, 1])
        # u_range <- 
        #   seq(min(emprden_hmc$x, -3), max(emprden_hmc$x, 3), 
        #       length.out = 1e3)
        # postden_hmc <- 
        #   u_likhd(u_range = u_range, 
        #           data.sim, lambda = lambda, p)
        
        plot(emprden_hmc$x, emprden_hmc$y/max(emprden_hmc$y),
             xlim = c(-3, 3), #range(u_range), 
             type= "l", main = task_title)
        # lines(u_range, postden_hmc/max(postden_hmc), 
        #       col = 2)
      }
      
      if (marg_eig_plot){
        if (samplr_name == "gibbs") {
          burn_idx <- 1:round(iter_max[2]/(2*record_freq))
          rec_gibbs <- rec
          plot(density(rec_gibbs$eig_list[-burn_idx]), 
               xlim = c(0, 4), ylim = c(0, 1),
               main = samplr_name)
          abline(v = svd(data.sim$theta)$d[2])
        } else if (samplr_name == "hmc_stan") {
          rec_hmc <- rec
          lines(density(rec_hmc$eig_list), col = 2)
        }
      }
      
      if (marg_dens_plot_slice){
        if (samplr_name == "slice") {
          rec_slice <- rec
        } else if (samplr_name == "hmc_stan") {
          rec_hmc <- rec
          par(mfrow = c(4, 4))
          idx_list <- 
            # cbind(sample(1:n, 4, replace = TRUE), 
            #       sample(1:p, 4, replace = TRUE))
            expand.grid(1:n, 1:k)
          
          for (i in 1:nrow(idx_list)){
            idx1 <- idx_list[i, 1]
            idx2 <- idx_list[i, 2]
            densest_slc <- 
              density(rec_slice$U[ , idx1, idx2])
            densest_hmc <- 
              density(rec_hmc$U[ , idx1, idx2])
            
            plot(densest_slc$x, densest_slc$y, type = "l",
                 xlab = "Theta", ylab = "density",
                 main = paste0("Theta[", idx1, ", ", idx2, "]"))
            lines(densest_hmc$x, densest_hmc$y, col = 2)
            
            # plot(rec_hmc$Theta[ , idx1, idx2], 
            #      type = "l", col = 2, 
            #      xlab = "Iteration",
            #      ylab = paste0("Theta[", idx1, ", ", idx2, "]"))
            # lines(rec_slice$Theta[ , idx1, idx2])
          }
        }
      }
      # save
      var_name <- paste0(family_name,
                         "_k", k, "_snr", snr,
                         "_", samplr_name)
      paste0("rec_", var_name, " <- rec") %>%
        parse(text = .) %>% eval()
      save(rec, file = paste0("./result/rec_", var_name, ".RData"))
      
      # if (samplr_name == "hmc_stan") {
      #   theta_true <- 
      #     apply(rec$Theta[round(iter_max[2]/2):(iter_max[2]), , ], c(2, 3), mean)
      #   U_true <- 
      #     apply(rec$U[round(iter_max[2]/2):(iter_max[2]), , ], c(2, 3), mean)
      #   V_true <-
      #     apply(rec$V[round(iter_max[2]/2):(iter_max[2]), , ], c(2, 3), mean)
      #   true_par$theta <- theta_true
      #   init_hmc <- list(theta = rec$theta[length(rec$time), , ], 
      #                    U = rec$U[length(rec$time), , ], 
      #                    V = rec$V[length(rec$time), , ])
      #   save(init_hmc, file = paste0("./result/theta_", var_name, ".RData"))
      # }
      #}
    }
  }
}
# dev.off()

# source("./example/glrm_pred.R")
