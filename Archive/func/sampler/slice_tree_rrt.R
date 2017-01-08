# Slice Sampler, sample by row, use rapid-exploring random tree

glrm_sampler_slice_rrt <- 
  function(Y, lambda, family_name, 
           init, config, rec, info, 
           # whether visualize tree when dim = 2
           visual_2d = FALSE){
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
      # Loops to Sample U
      for (i in 1:n){
        # grow tree, stop after first step out
        tree <- init_tree(U_cur[i, ])
        acc <- TRUE
        
        # optionally, visualize
        if (visual_2d){
          grid_power <- 10^(0.2)
          grid <- list(seq(-1, 1, 1e-2)/grid_power + U_cur[1], 
                       seq(-1, 1, 1e-2)/grid_power + U_cur[2])
          grid_dens <- 
            outer(grid[[1]], grid[[2]],
                  Vectorize(function(x1, x2){
                    U_try <- matrix(c(x1, x2), ncol = 2)
                    target_lik(U_try, V_cur, T_suff, lambda)
                  }))
          
          contour(grid[[1]], grid[[2]], grid_dens)
          contour(grid[[1]], grid[[2]], grid_dens < slice_lik, add = TRUE)
          points(U_cur[1], U_cur[2], pch = 19)
        }
        
        while(acc | (nrow(tree) < 2)){
          # propose node
          U_dir <- prop_dir_rrt(tree[nrow(tree), 1:k])
          new_node <- 
            prop_node_rrt(U_dir, tree, 
                      edge_max = edge_max, 
                      random = TRUE)
          # examine node
          U_prop <- U_cur
          U_prop[i, ] <- as.matrix(new_node[, 1:k])
          acc <- 
            (target_lik(U_prop, V_cur, T_suff, lambda) < slice_lik)
          # accept/reject node
          if (acc){
            tree <- rbind(tree, new_node)
          } 
          
          if (visual_2d){
            parent_node <- tree[new_node$parent, ]
            
            # points(parent_node[, 1:k], cex = 2)            
            points(new_node[, 1:k], pch = 19, 
                   col = ifelse(acc, "red", "grey"))
            # segments(x0 = parent_node[, 1], y0 =  parent_node[, 2],
            #          x1 = U_dir[, 1], y1 = U_dir[, 2],
            #          lwd = 0.5, lty = 2)
            segments(x0 = parent_node[, 1], y0 =  parent_node[, 2],
                     x1 = new_node[, 1], y1 = new_node[, 2],
                     lwd = 0.5, col =  ifelse(acc, "black", "grey"))
            Sys.sleep(0.1)
          }
          
          reflect <- FALSE
          if (reflect & (!acc)){
            # if step out, then try reflection once, produce new_node2
            # if still not work then give up
            parent_node <- tree[new_node$parent, ]
            prop_grad <-
              grad_hmc_U(T_suff, U_prop, V_cur, lambda, d1)[i, ]
            new_node2 <-
              flex_node(new_node, parent_node, prop_grad)
            U_prop[i, ] <- as.matrix(new_node2[, 1:k])
            acc2 <-
              (target_lik(U_prop, V_cur, T_suff, lambda) < slice_lik)
            
            if (acc2){
              tree <- rbind(tree, new_node2)
            }
            if (visual_2d){
              points(new_node2[, 1:k], pch = 19,
                     col = ifelse(acc, "red", "grey"))
              segments(x0 = parent_node[, 1],
                       y0 =  parent_node[, 2],
                       x1 = U_dir[, 1], y1 = U_dir[, 2],
                       lwd = 0.5)
              # segments(x0 = new_node[, 1], y0 =  new_node[, 2],
              #          x1 = new_node2[, 1], y1 = new_node2[, 2],
              #          lwd = 0.5, col =  ifelse(acc, "black", "grey"))
              #Sys.sleep(0.1)
            }
            
          }
          
          # user_input <- readline("Continue? (y/n)")
          # if (user_input == "y"){
          #   next
          # } else {
          #   break
          # }
        }
        
        # after tree grown, select proposal
        node_id <- 
          ifelse(
            nrow(tree) == 2, 2, 
            sample(x = 2:nrow(tree), size = 1)
          )
        U_cur[i, ] <- as.matrix(tree[node_id, 1:k])
        
        if (visual_2d){
          points(U_cur[i, 1], U_cur[i, 2], cex = 2)
          text(U_cur[i, 1] + 0.05, U_cur[i, 2] + 0.05, labels = node_id)
        }
      }
      
      ####  V  #################
      acc_V <- rep(NaN, p)
      # warning("only V updated")
      # Loops to Sample V
      # for (j in 1:p){
      #   V_old <- V_cur
      #   
      #   # generate A', A'', B
      #   Theta_cur <- U_cur %*% t(V_cur)
      #   A_d1 <- d1(Theta_cur)
      #   A_d2 <- d2(Theta_cur)
      #   B <- T_suff - A_d1 + A_d2 * Theta_cur
      #   lnr_coef_v <- t(B) %*% U_cur # p x k, each row lnr coef for v_j
      #   
      #   # sample for i^th row of U
      #   V_prop <- V_cur
      #   sigma <- solve(
      #     t(U_cur) %*% diag(A_d2[, j]) %*% U_cur +
      #       lambda * diag(k))
      #   mu <- sigma %*% lnr_coef_v[j, ]
      #   V_prop[j, ] <-
      #     rmvnorm(1,
      #             mean = mu,
      #             sigma = sigma)
      #   
      #   # V_prop[j, ] <-
      #   #   rmvnorm(1,
      #   #           mean = V_cur[j, ],
      #   #           sigma = sigma * diag(k))
      #   
      #   # metroplis step for V
      #   acc_prob <-
      #     acc_prob_V(U_cur, U_cur, V_prop, V_old, j,
      #                lambda, family, T_suff)
      #   # warning("no rejection for V")
      #   acc_V[j] <- (runif(1) < acc_prob)
      #   
      #   if (acc_V[j])
      #     V_cur <- V_prop
      # }
      
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