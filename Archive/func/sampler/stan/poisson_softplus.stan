data {
  int<lower=1> N;               // number of population
  int<lower=1> P;               // number of species
  int<lower=1> K;               // dimension of latent factors 
  int<lower=0> Y[N,P];          // data matrix of order [N,P]    
  
  real<lower=0> lambda_u;       // penalty for prior
  real<lower=0> lambda_v;  
}

transformed data {
  // helper parameters for constructing prior
  vector[K] mu_K;    
  cov_matrix[K] Sig_K;
  cov_matrix[P] Sig_P;

  mu_K = rep_vector(0.0, K);
  Sig_K = diag_matrix(rep_vector(lambda_u, K));
  Sig_P = diag_matrix(rep_vector(lambda_v, P));
}

parameters {
  matrix[N, K] U;   // factor loading for population
  matrix[P, K] V;   // latent factor for species     
  /*  
  matrix[N, P] e;     // factor noise
  */  
}

transformed parameters{
  matrix[N, P] Theta;             // parameter matrix
  //Theta = U * V' + e;   
  Theta = log(1 + exp(U * V')) ;  // softmax link function
}

model {
  // the priors 
  /*
  for(i in 1:N){
    e[i] ~ multi_normal(mu_K, Sig_K); 
    U[i] ~ multi_normal(mu_K, Sig_K);  
  }
  */
  
  for(i in 1:P){
    V[i] ~ multi_normal(mu_K, Sig_K);   
  }

  //The likelihood
  for (j in 1:P){  
    for(i in 1:N){
      Y[i, j] ~  poisson(Theta[i, j]); 
   }
  }
}
