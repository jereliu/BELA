library(VineCopula)
library(clusterGeneration)

biom_sample <- 
  function(n_sample = 50, n_OTU = 20, 
           n_obs = 100,
           pois_eps = 1e-3, pois_sample = c(1,1),
           cop_type = "Gauss", cop_par, 
           inv_cdf = qnorm, 
           n_block = 5, 
           theta = 0.5
  ){
    # INPUT: 
    # > data parameters
    # >> n_sample: num of total samples, 
    #              each sample contain OTU counts
    # >> n_OTU: num of possible OTU categories
    # >> n_obs: ttl num of OTU count per sample
    
    #
    # > model parameters
    # >> pois_eps:    posterior pois process intensity pois_eps = alpha/m
    # >> pois_sample: effective posterior adj to poisson process 
    #                 (observed (n_1, n_2) success/failures)  
    # >> cov_type:    Copula type
    # >> cov_par:     Copula parameters
    
    # DOCUMENTATION:
    # Generate samples from correlated random processes 
    # (index OTU with i, and sample with j)
    # using below steps:
    # 1. Generate sigma [I x 1 vec] from Pois Process (pois_alpha)
    # 2. Generate Q [I x J matrix]:
    #   2.1 Generate correlation structure from Copula (cop_*)
    #   2.2 Generate sample using inverse CDF
    # 3. Assemble, 
    #   P_{ji} = sigma_i * Q_{ji}
    # 4. Generate n_{ji} ~ P_{ji}
    
    
    #### 1. Generate sigma ####
    # generate poisson process with par (a=pois_rate, b=0.5+pois_rate)
    # (don't use prior since its degenerate)
    # using beta approximation
    sigma <- 
      rbeta(n_OTU, 
            pois_eps + pois_sample[1], 
            0.5 - pois_eps + pois_sample[2]) %>% sort %>%
      matrix(., ncol = 1)
    
    #### 2. Generate Q_{ji} ####
    Q_out <- rmvtnorm_sparse(
      n_OTU, n_sample, 
      n_block = n_block, 
      theta = theta)
    Q <- Q_out$Q %>% t
    Sigma <- Q_out$Sigma
    
    #### 3. Assemble ####
    P <- 
      # assemble
      apply(Q, 1, function(q) pmax(q, 0)^2*sigma) %>%
      # normalize
      apply(2, function(p_j) p_j/sum(p_j)) %>% t
    
    #### 4. Generate ####
    sample <- 
      sapply(
        1:n_sample, 
        function(j) 
          rmultinom(1, n_obs, P[j, ]))
    
    # re-structure then return 
    sample <- 
      sample %>% t %>% as.data.frame %>% 
      set_names(paste0("OTU_", 1:n_OTU))
    
    list(N = sample, Q = Q, Sigma = Sigma)
  }