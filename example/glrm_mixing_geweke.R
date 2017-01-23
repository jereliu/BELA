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
k <- 1
family_name <- c("gaussian", "poisson", "poisson_softplus")[2]
samplr_name <- "hmc_stan"
snr <- 100

rand_seeds <- list(data = 4200, samplr = 1300)
pred_plot <- TRUE
prior_dens_plot <- FALSE

for (family_name in c("gaussian", "poisson")){
  #for (snr in c(100, 10, 1, 0.5)){
  for (k in c(1, 2, 5, 10, 15, 20)[c(1)]){
    lambda = 1
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
    for (samplr_name in c("hmc_stan_debug", "gibbs_debug")[c(1:2)]){
      # choose iter based on 
      if (length(grep("gibbs", samplr_name)) > 0) {
        iter_max <- c(1e5, 1e3) # 1e3)
      } else {
        # if sampler name contain "hmc"...
        iter_max <- c(1e5, 1e3) # 1e4)
      }
      
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
        glrm(Y, lambda = lambda, k = k, 
             true_par = true_par,
             init = init, init_MAP = init_MAP,
             samplr_name = samplr_name,
             family_name = family_name,
             iter_max = iter_max, 
             record_freq = 1,
             time_max = 60, 
             step_size = c(step_optim, step_sampl), 
             frog_step = 5,
             rotn_freq = iter_max[2],
             mmtm_freq = iter_max[2])
      
      #rec$mix_score <- prior_mixing_score(rec, 1/sqrt(lambda))
      
      #### 3. Evaluation ####
      # plot
      task_title <- 
        paste0(
          family_name, " ", 
          samplr_name, " k=", k, " snr = ", snr)
      
      if (pred_plot){
        #plot(rec$time, rec$error, type = "l", 
        #     main = paste0("sample error ", task_title))
        time_size <- length(rec$time)
        plot(rec$pred_error, type = "l", 
             main = paste0("prediction error ", task_title))
        
        # objective function
        plot(rec$time[(time_size-length(rec$obj)):time_size], 
             rec$obj, 
             type = "l", xlab = "Time (min)",
             main = paste0("obj ", task_title))
        T_y <- glmr_family(family_name)$sufficient(Y)
        abline(h = 
                 glmr_family(family_name)$negloglik(T_y, (rec$init$U) %*% t(rec$init$V)) + 
                 (1/phi_sd) * (sum(rec$init$U^2) + sum(rec$init$V^2))
        )
      }
      
      if (prior_dens_plot){
        # exact posterior density matching for U[1, 1]
        emprden <- density(rec$U[, 1, 1])
        u_range <- range(c(emprden$x, -3, 3))
        u_range <- seq(min(u_range), max(u_range), length.out = 1e3)
        
        plot(emprden, xlim = range(u_range), type = "l", main = task_title)
        lines(u_range, dnorm(u_range), col = 2)
      }
      # save
      var_name <- paste0(family_name,
                         "_k", k, "_snr", snr,
                         "_", samplr_name)
      paste0("rec_", var_name, " <- rec") %>%
        parse(text = .) %>% eval()
      save(rec, file = paste0("./result/mixing/rec_", var_name, ".RData"))
      
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
    }
  }
}



post_diag <- FALSE

if (post_diag){
  
  # norm 
  norm <- sapply(1:iter_max[2], function(i) sum(rec$U[i, , ]^2) + sum(rec$V[i, , ]^2))
  plot(norm, type = "l")
  
  # likelihood
  likhd <-
    sapply(1:iter_max[2], function(i) 
      glmr_family(family_name)$negloglik(T_y, (rec$U[i, , ]) %*% t(rec$V[i, , ]))
    )
  plot(likhd, type = "l")
  abline(
    glmr_family(family_name)$negloglik(T_y, 
                                       (rec$init$U) %*% t(rec$init$V))
  )
  
  # objective
  # likelihood
  l_2 <- 
    sapply(1:iter_max[2], function(i) 
      sqrt(mean((rec$true_theta - (rec$Theta[i, , ]))^2)/
             mean(rec$true_theta^2))
    )
  plot(l_2, type = "l")
}
