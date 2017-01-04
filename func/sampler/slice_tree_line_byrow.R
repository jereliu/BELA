# Slice Sampler, sample by row using doubling procedure on line sampler
glrm_sampler_slice <- 
  function(Y, lambda, family_name, 
           init, config, rec, info, 
           # whether visualize tree when dim = 2
           freq_slicing = FALSE,
           dir_hess = TRUE,
           hit_n_run = FALSE,
           hit_n_run_exact = FALSE,
           visual = FALSE){
    # unpack family properties
    n <- info$n
    p <- info$p
    k <- info$k
    true_theta <- info$true_par$theta
    
    family <- glrm_family(family_name)
    T_suff <- family$sufficient(Y)
    d1 <- family$partition$d
    d2 <- family$partition$d2
    negloglik <- family$negloglik
    
    target_lik <- function(U, V, T_suff, lambda){
      negloglik(T_suff, U %*% t(V)) + 
        (lambda/2) * (sum(V^2) + sum(U^2))
    }
    
    # unpack mcmc parameters
    iter_max <- config$sampler$iter_max
    record_freq <- config$record_freq
    time_max <- config$time_max
    rotn_freq <- config$sampler$rotn_freq
    edge_max <- config$sampler$edge_max 
    
    
    U_cur <- init$U
    V_cur <- init$V
    Theta_cur <- U_cur %*% t(V_cur)
    time0 <- proc.time()[3]
    
    # initiate sampler
    pb <- txtProgressBar(min = 1, max = iter_max, style = 3)
    
    for (iter in 1:iter_max) {
      setTxtProgressBar(pb, iter)
      # if (iter %% rotn_freq == 0) {
      #   R <- rortho(k)
      #   U_cur <- U_cur %*% R
      #   V_cur <- V_cur %*% R
      # }
      
      U_old <- U_cur; V_old <- V_cur
      ####  y  ################
      slice_lik <- 
        target_lik(
          U_cur, V_cur, T_suff, lambda) + rexp(1)
      
      
      ####  U  ################
      acc_U <- rep(NaN, n)
      Theta_cur <- U_cur %*% t(V_cur)
      
      # Loops to Sample U
      for (i in 1:n){
        if (freq_slicing){
          slice_lik <- 
            target_lik(
              U_cur, V_cur, T_suff, lambda) + rexp(1)
        }
        
        # grow tree, stop after first step out
        tree <- matrix(U_cur[i, ], nrow = 1)
        acc <- TRUE
        
        # optionally, visualize
        if (visual){
          if (k == 1){
            grid_power <- 10^(0.2)
            grid <- seq(-1, 1, 1e-2)/grid_power + U_cur
            grid_dens <- 
              sapply(grid,
                     function(x){
                       U_try <- as.matrix(x)
                       target_lik(U_try, V_cur, T_suff, lambda)
                     })
          } else if (k == 2){
            grid_power <- 10^(0.5)
            grid <- list(seq(-1, 1, 1e-2)/grid_power + U_cur[1], 
                         seq(-1, 1, 1e-2)/grid_power + U_cur[2])
            grid_dens <- 
              outer(grid[[1]], grid[[2]],
                    Vectorize(function(x1, x2){
                      U_try <- matrix(c(x1, x2), ncol = 2)
                      target_lik(U_try, V_cur, T_suff[i, ], lambda)
                    }))
          }
        }
        
        # 1. select direction
        #(TODO: optionally, 
        #       sample principal direction using 
        #       eigenvectors of hessian matrix)
        if (dir_hess) {
          Hess_U <-
            t(V_cur) %*% diag(d2(Theta_cur[i, ])) %*% V_cur +
            lambda * diag(k)
          eig_Hess <- eigen(Hess_U)
          dir_id <- which(rmultinom(1, 1, 1/eig_Hess$values)==1)
          U_dir <- eig_Hess$vectors[, dir_id]
        } else {
          U_dir <- rmvnorm(1, rep(0, k))
          U_dir <- U_dir/sqrt(sum(U_dir^2))
        }
        
        if (visual){
          library(shape)
          if (k==1){
            plot(grid, grid_dens, type = "l")
            abline(h = slice_lik)
            points(U_cur, slice_lik, pch = 19)
          } else if (k==2){
            contour(grid[[1]], grid[[2]], grid_dens)
            contour(grid[[1]], grid[[2]], grid_dens < slice_lik, add = TRUE)
            points(U_cur[1], U_cur[2], pch = 19)
            
            coord1 <- as.matrix(tree[1, 1:2] - 10 * U_dir)
            coord2 <- as.matrix(tree[1, 1:2] + 10 * U_dir)
            
            segments(coord1[1], coord1[2], 
                     coord2[1], coord2[2], 
                     col = rgb(0,0,0,0.4))
            
            # center <- c(-0.4, -0.3)
            # Arrows(center[1], center[2], 
            #        center[1] + 0.2 * eig_Hess$vectors[1, 2], 
            #        center[2] + 0.2 * eig_Hess$vectors[2, 2])
            # Arrows(center[1], center[2], 
            #        center[1] + 0.1 * eig_Hess$vectors[1, 1], 
            #        center[2] + 0.1 * eig_Hess$vectors[2, 1])
            
          }
        }
        
        # 2. select step size, use Hessian (diagonal)
        # obtain directional gradient
        
        # grad_U <-
        #   grad_hmc_U(T_suff, U_cur, V_cur, lambda, d1)
        step_dir <- 0.05 * t(U_dir)
        # / as.numeric(grad_U %*% U_dir)
        
        # Hess_U <-
        #   t(V_cur) %*% diag(d2(Theta_cur[i, ])) %*% V_cur +
        #   lambda * diag(k)
        # step_dir <- solve(Hess_U, U_dir)
        # 
        # U_dir_step <- 
        #   t(eigen(Hess_U)$vectors[, 2])     
        # 
        # if (visual_2d){
        #   Arrows(
        #     x0 = U_cur[1], y0 = U_cur[2], 
        #     x1 = U_cur[1] + U_dir_step[1], 
        #     y1 = U_cur[2] + U_dir_step[2], pch = 19)
        # }
        # 
        
        # Hess_U_diag <-
        #   apply(V_cur, 2,
        #         function(row)
        #           sum(row^2 * d2(Theta_cur[i, ]))
        #   ) + lambda
        # step_mat <- 1/Hess_U_diag
        # step_dir <- step_mat * U_dir
        if (hit_n_run) {
          
          U_cur[i, ] <- 
            sampleLevelSet(
              var = "U", idx = i, 
              U_cur, V_cur, T_suff, lambda,
              U_dir, V_dir = NULL, 
              slice_lik, target_lik
            )
          
        } else if (hit_n_run_exact) {
          # 3.1 hit and run!
          # find the two solutions where the line "hits" level set
          solMat <- 
            solveLevelSet_bisect(
              var = "U", idx = i, 
              U_cur, V_cur, T_suff, lambda,
              U_dir, V_dir = NULL, 
              slice_lik, target_lik
            )
          
          if (visual){
            if (k == 1){
              points(solMat[1, ], slice_lik, pch = 19, col = 2)
              points(solMat[2, ], slice_lik, pch = 19, col = 2)
            } else if (k == 2){
              points(solMat[1, 1], solMat[1, 2], pch = 19, col = 2)
              points(solMat[2, 1], solMat[2, 2], pch = 19, col = 2)
              segments(solMat[1, 1], solMat[1, 2], 
                       solMat[2, 1], solMat[2, 2], col = 2)
            }
          }
          
          # uniform sampling on direction
          runif_line <- runif(1)
          U_cur[i, ] <-
            runif_line * solMat[1, ] +
            (1 - runif_line)* solMat[2, ]
        } else {
          # 3.2 or, build doubling tree till step out
          new_node_count <- 1
          
          while(acc){
            # propose node
            #new_node_count <- new_node_count * 2
            tree_dir <- sample(c(-1, 1), 1)
            
            new_node <-
              prop_node_line(new_node_count, tree,
                             tree_dir, step_dir)
            
            # examine node
            U_prop <- U_cur
            U_prop[i, ] <- as.matrix(new_node[new_node_count, ])
            acc <-
              (target_lik(U_prop, V_cur, T_suff, lambda) < slice_lik)
            
            # accept/reject node
            if (acc){
              # if accept, grow tree according to direction
              if (tree_dir == 1){
                tree <- rbind(tree, new_node)
              } else {
                tree <-
                  rbind(new_node[new_node_count:1, ], tree)
              }
            }
            
            if (visual){
              if (k == 1){
                # points(parent_node[, 1:k], cex = 2)
                for (new_node_i in new_node){
                  points(new_node_i, slice_lik, pch = 19,
                         col = ifelse(acc, "red", "grey"))
                }
              } else if (k == 2){
                # points(parent_node[, 1:k], cex = 2)
                points(new_node, pch = 19, cex = 1,
                       col = ifelse(acc, "red", "grey"))
                Arrows(U_cur[1], U_cur[2], 
                       new_node[1], new_node[2], 
                       arr.length = 0.15,
                       arr.width = 0.1)
              }
              Sys.sleep(1)
            }
            
            # user_input <- readline("Continue? (y/n)")
            # if (user_input == "y"){
            #   next
            # } else {
            #   break
            # }
          }
          
          # select tree node
          node_id <- sample(x = 1:nrow(tree), size = 1)
          U_cur[i, ] <- as.matrix(tree[node_id, ])
        }
        
        # after tree grown, select proposal
        # # select tree node
        # node_id <- sample(x = 1:nrow(tree), size = 1)
        # U_cur[i, ] <- as.matrix(tree[node_id, ])
        
        if (visual){
          if (k == 1){
            points(U_cur[i, ], slice_lik, cex = 2)
          } else if (k == 2){          
            points(U_cur[i, 1], U_cur[i, 2], cex = 2)
            text(U_cur[i, 1] + 0.05, U_cur[i, 2] + 0.05, labels = node_id)
          }
        }
      }
      
      ####  V  #################
      acc_V <- rep(NaN, p)
      
      for (j in 1:p){
        
        # grow tree, stop after first step out
        tree <- matrix(U_cur[i, ], nrow = 1)
        acc <- TRUE
        
        # 1. select direction
        #(TODO: optionally, 
        #       sample principal direction using 
        #       eigenvectors of hessian matrix)
        if (dir_hess) {
          Hess_V <-
            t(U_cur) %*% diag(d2(Theta_cur[, j])) %*% U_cur +
            lambda * diag(k)
          eig_Hess <- eigen(Hess_V)
          dir_id <- which(rmultinom(1, 1, 1/eig_Hess$values)==1)
          V_dir <- eig_Hess$vectors[, dir_id]
        } else {
          V_dir <- rmvnorm(1, rep(0, k))
          V_dir <- V_dir/sqrt(sum(V_dir^2))
        }
        
        # 2. select step size, use Hessian (diagonal)
        # obtain directional gradient
        
        # grad_U <-
        #   grad_hmc_U(T_suff, U_cur, V_cur, lambda, d1)
        step_dir <- 0.05 * t(V_dir)
        # / as.numeric(grad_U %*% U_dir)

        V_cur[j, ] <- 
          sampleLevelSet(
            var = "V", idx = j, 
            U_cur, V_cur, T_suff, lambda,
            U_dir = NULL, V_dir, 
            slice_lik, target_lik
          )
      }
      
      
      # record
      if (iter %% record_freq == 0){
        rec$U[iter/record_freq, , ] <- U_cur
        rec$V[iter/record_freq, , ] <- V_cur
        rec$Theta[iter/record_freq, , ] <- U_cur %*% t(V_cur)
        
        rec$acc[iter/record_freq, ] <- 
          c(mean(acc_U), mean(acc_V))
        rec$obj[iter/record_freq] <- 
          negloglik(T_suff, U_cur %*% t(V_cur)) + 
          (lambda/2) * (sum(V_cur^2) + sum(U_cur^2))
        
        rec$time[iter/record_freq] <- 
          (proc.time()[3] - time0)/60
      }
      
      # stop if time used up
      # if ((proc.time()[3] - time0)/60 > time_max){
      #   break
      # }
    }
    
    # return
    rec
  }