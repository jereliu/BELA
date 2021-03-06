boyu_sample <- 
  function( lscounts, n, p, m, K, 
            a.er, b.er, sigma.value, sigma.weight, 
            a=2, 
            type = c("mult", "pois")[1], 
            link = "exp"){
    # simulate from the prior
    link_func <- 
      paste0("link_", link) %>% 
      parse(text = .) %>% eval
    
    # INPUT #
    # lscounts: list of total counts, will generate one dataset for each total counts in the list
    # n: number of biological samples
    # p: number of species
    # m: number of factors
    # K: number of blocks in the underlying covariance matrix
    # allow for misspecification test
    # biological samples are independent between blocks
    # a.er, b.er: hyperparameter in the prior of e_r
    # sigma.value, sigma.weight: discretized prior of small sigma
    # a: random weights of species are generated by sigma*Q^a*(Q>0)
    
    # OUTPUT #
    # List with entries: sigma, Q, X.tru, Y.tru, er, data
    # sigma: vector of simulated small sigma's, length = p
    # Q: matrix of normal latent variables, dim = p*n
    # X.tru, Y.tru: latent factor and loadings, dim(X.tru) = m*p, dim(Y.tru) = m*n
    # er: error term
    # data: a list of matrix, length = length( lscounts )
    # each matrix is p*n.
    
    sigma = sample( sigma.value, p, replace = T, sigma.weight )
    X = matrix( rnorm( p*m ), nrow = m )
    #split the m factors into K blocks
    factor.indx.group = split( 1:m, cut( 1:m, breaks = K ) )
    pop.indx.group = split( 1:n, cut( 1:n, breaks = K ) )
    Y = matrix( 0, nrow = m, ncol = n )
    
    for( x in 1:K){
      raw = matrix( rnorm( length(factor.indx.group[[x]])*length(pop.indx.group[[x]]) ), 
                    ncol = length(pop.indx.group[[x]]) )
      Y[factor.indx.group[[x]],pop.indx.group[[x]]] = raw
    }
    
    er = 1/rgamma( 1, a.er, b.er)
    Q = t( apply( t(Y)%*%X, 2, function(x) rnorm( length(x), mean = x, sd = sqrt(er) ) ) )
    final.weights = sigma*link_func(Q)
    
    if (type == "mult"){
    data = 
      lapply(lscounts, 
             function(counts) 
               apply(final.weights, 2, 
                     function(x) rmultinom( 1, counts, prob = x ) ) )[[1]]
    } else if (type == "pois"){
      # final.weights = 
      #   (sigma %*% t(rep(lscounts, n))) *
      #   link_func(Q)
      data = round(final.weights)
    }
    
    return( list( sigma = sigma, Q = Q, data = data,
                  X.tru=X, Y.tru = Y, er = er, 
                  final.weights = final.weights ) )
  }
