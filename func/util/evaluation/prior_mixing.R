prior_mixing_score <- 
  function(rec){
    iter_max <- dim(rec$U)[1]
    score_out <- 
      sapply(3:iter_max, function(index) 
        pmscore_index(index, rec, prior_sd)
      )
    score_out
  }


pmscore_index <- 
  function(index_max, rec, prior_sd = 1){
    # quant <- pnorm(rec$U[1:index_max, 1, 1], sd = prior_sd)
    # ecdf_quant <- ecdf(quant)
    # # test for normality for Phi^{-1}(quantile)
    # # log(1 - ks.test(quantile, "pnorm", 0, 1)$p.value)
    # # shapiro.test(qnorm(quantile))$p.value
    # 
    # # check deviation from uniform cdf
    # prob <- seq(0, 1, 1e-3)
    # max(abs(ecdf_quant(prob) - prob))
    dens_est <- density(rec$U[1:index_max, 1, 1])
    dens_real <- dnorm(dens_est$x, sd = prior_sd)
    max(abs(dens_est$y - dens_real))
  }