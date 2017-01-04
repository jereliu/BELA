library(KSD)

#### 1. geweke, converge to prior ####
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

#### 2. Cook, quantile method ####
prior_mixing_cook <- 
  function(sample, conf.level = 0.95){
    # taken from 
    # http://www.stat.umn.edu/geyer/old/5601/examp/kolmogorov.html#one-ci
    n <- length(sample)
    stat <- 
      (ecdf(sample)(seq(0, n)/n) - seq(0, n)/n) %>% 
      abs %>% max
    stat
  }

pksone <- function(x) {
  k <- seq(1:20)
  1 + 2 * sum( (-1)^k * exp(- 2 * k^2 * x^2))
}

qksone <- function(p, n) {
  foo <- function(x) pkolm(x, n) - p
  uniroot(foo, c(1e-10, 1))
}


#### 3. Geweke, KS-test ####
prior_mixing_geweke_perm <- 
  function(sample, prior, n_perm = 100){
    # taken from 
    # http://www.stat.umn.edu/geyer/old/5601/examp/kolmogorov.html#one-ci
    n <- length(sample)
    
    stat_perm <- 
      sapply(
        1:n_perm, 
        function(i) {
          perm_prob <- rbinom(n, 1, 0.5)
          
          sample_1 <- (perm_prob==1)*sample + (perm_prob==0)*prior
          sample_2 <- (perm_prob==0)*sample + (perm_prob==1)*prior
          
          ks.test(sample_1, sample_2)$statistic
        }
      )
    
    stat_obs <- ks.test(sample, prior)$statistic
    
    mean(stat_obs > stat_perm)
  }

prior_mixing_geweke_exact <- 
  function(sample, prior, n_perm = 100){
    ks.test(sample, "pnorm")$"statistic"
  }

prior_mixing_geweke_stat <- 
  function(sample, prior){
    # taken from 
    # http://www.stat.umn.edu/geyer/old/5601/examp/kolmogorov.html#one-ci
    ks.test(sample, prior)$statistic
  }

#### 4. Jeremiah ####
prior_mixing_jere <- 
  function(theta_container, length = 10){
    n_iter <- ncol(theta_container) - 1
    prior <- theta_container[, 1]
    sample_iter <- 
      round(seq(2, n_iter, length.out = length))
    
    quantile_byIter <- 
      sapply(sample_iter, 
             function(i)
               rowMeans((theta_container[, 2:(i+1)] > prior))
      ) %>% set_colnames(sample_iter)
    
    stat_byIter <- 
      apply(quantile_byIter, 2, 
            function(sample){
              ks.test(sample, "punif")$statistic
            }
      )
    # thres_byIter <- 
    #   sapply(sample_iter,
    #         function(n) qksone(0.95, n)$root 
    #   )
    
    plot(log(stat_byIter), type = "l")
    lines(thres_byIter, col = 2)
    
    stat_byIter
  }


#### 5. Geweke, KL-distance ####

mixing_KL_dist <- 
  function(par_list, loglikfunc){
    # X is aggregated sample at iter n, a list contain U and V
    # loglikfunc is the likelihood
    self_ep <- entropy(discretize(par_list))
    cross_ep <- 
      apply(par_list, 1, function(y_u_v){
        fac_dim <- (length(y_u_v)-1)/2
        y <- y_u_v[1]
        u <- y_u_v[2:(fac_dim+1)]
        v <- y_u_v[(fac_dim+2):(2*fac_dim+1)]
        loglikfunc(u, v, y)
      }) %>% mean
    self_ep + cross_ep
  }

#### 6. Stein ####
mixing_stein <- 
  function(negrad_func, theta_list, theta_row = NULL, 
           # hyperparameter for Gaussian kernel
           width = -1, nboot = 1000, ...
  ){
    # gaussian kernel with bandwidth equal to 1/alpha
    neggrad_list <- 
      apply(theta_list, 1, 
            function(theta){
              negrad_func(matrix(theta, nrow = theta_row), ...)
            }
      ) %>% t 
    
    KSD(theta_list, -neggrad_list, 
        width = width, nboot = nboot)
  }

result <-
  mixing_stein(grad_neglik_S,
               S_cur, theta_row = sum(dim(T_suff)),
               nboot = 1000,
               T_suff = T_suff,
               lambda = lambda,
               dist_family = family)
