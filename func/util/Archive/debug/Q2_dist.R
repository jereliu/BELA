Q2_dist <- 
  function(info, N){
    dist <- 
      apply(N, 1, function(x) x/sum(x)) - 
      apply(info$mean$Q2, 1, function(x) x/sum(x))
    mean(abs(dist))
  }

plot_diag <- FALSE
iter_idx <- 150:211

if (plot_diag){
  res_name <- "Q2"
  res_type <- "par"
  # graphically examing change in mean_Q2 over iteration
  data_plot_raw <- iter[[res_type]][[res_name]]
  
  if (length(dim(data_plot_raw)) == 3){
    data_plot <- 
      data_plot_raw[iter_idx, , ]
  } else if (length(dim(data_plot_raw)) == 2) {
    data_plot <- 
      data_plot_raw[iter_idx, ]
  }
  
  iter_cur <- dim(data_plot)[1]
  plot(0, type = "n", 
       xlim = c(1, iter_cur), 
       ylim = range(data_plot), 
       main = paste(res_type, res_name)
  )
  
  if (length(dim(data_plot)) == 3){
    for (j in 1:J){
      for (i in 1:I){
        lines(1:iter_cur, 
              data_plot[, j, i], 
              col = 
                c(rgb(0,0,0,0.3), 
                  rgb(1,0,0,0.3)
                )[2 - (N[j, i]>0)]
        )
      }
    }
  } else if (length(dim(data_plot)) == 2) {
    for (j in 1:ncol(data_plot)){
      lines(1:iter_cur, data_plot[, j], 
            col = rgb(0,0,0,0.5))
    }
  }
}

if (plot_diag){
  mode_est <- info$plot$Q1$mode_est
  normalizer <- info$plot$Q1$normalizer
  range_max <- info$plot$Q1$range_max
  s_j <- prior$Q$Sg_cond
  sd <- sqrt(s_j)
  n_ij <- info$stat$n_ij
  
  lambda_Q1 <- lambda$Q1
  lambda_Q2 <- lambda$Q2
  
  file_dir <- "./temp_plot/"
  # graphically examing mean for Q2 est 
  for (j in 1:J){
    for (i in 1:I){
      pdf(paste0(file_dir, "|Q|_+^2 ", j, "-", i, ", N = ", N[j, i], ".pdf"))
      Q <- seq(0, mode_est[j, i]+ 50*sd[j], 1e-3)
      kern_val <- 
        Q_kernel(
          Q,
          lambda_Q1 = lambda_Q1[j, i],
          lambda_Q2 = lambda_Q2[j, i],
          n_ij = n_ij[j, i],
          s_j = s_j[j],
          log_const_adj = normalizer[j, i]
        )
      idx <- which(kern_val > 1e-10)
      plot(pmax(0, Q)[idx]^2,
           kern_val[idx],
           type = "l", ylab = "", 
           main = paste0("|Q|_+^2 ", j, "-", i, ", N = ", N[j, i]))
      abline(v = res[j, i], col = 2)
      dev.off()
    }
  }
}

