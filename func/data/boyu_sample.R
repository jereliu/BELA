boyu_sample <- 
  function(lscounts, n, p, m, K, 
           a.er, b.er, 
           sigma.value, sigma.prob, 
           strength)
  {
    #we consider K block simulation stucture
    #browser()
    
    #### 2.1. Q ----
    sigma <-  
      sample(sigma.value, p, replace = T, 
             sigma.prob )
    
    #### 2.2. X and Y ----
    X <-  matrix( rnorm( p*m ), nrow = m )
    
    #split the m factors into K blocks
    pop.indx.group <-  split( 1:n, cut( 1:n, breaks = K ) )
    Y <-  matrix( 0, nrow = m, ncol = n )
    
    for( x in 1:K){
      Y[,pop.indx.group[[x]]] <- 
        (2*x - (K+1))*strength + 
        matrix( rnorm( m*length(pop.indx.group[[x]]) ), 
                ncol = length(pop.indx.group[[x]]) )
    }
    
    #### 2.3. Q ----
    er <-  1/rgamma(1, a.er, b.er) # sample Q error
    Q <-  
      apply(t(Y)%*%X, 2, 
            function(x) 
              rnorm(length(x), mean = x, sd = sqrt(er) ) 
      ) %>% t
    
    #### 3. Final Assemble and Sampling ====
    final.weights <-  sigma*(Q*(Q>=0))
    data <-  
      lapply(lscounts, 
             function(counts) 
               apply(final.weights, 2, 
                     function(x) 
                       rmultinom( 1, counts, prob = x ) 
               )
      )
    
    #### 4. Return ====
    return( 
      list(sigma = sigma, Q = Q, data = data,
           X.tru=X, Y.tru = Y, er = er) 
    )
  }
