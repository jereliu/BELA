data {
  int<lower=1> N;                // number of 
  int<lower=1> P;                // number of 
  int<lower=1> K;              // number of latent dimensions 
  real<lower=0> lambda_u;      // penalty for prior
  real<lower=0> lambda_v;  
  matrix[N,P] Y;                 // data matrix of order [N,P]  

  row_vector[K] V[P];   // latent factor for species
}

parameters {    
  row_vector[K] U[N];   // factor loading for population
  matrix[N, P] e;   // factor noise
  vector[K] mu_K;  
  cov_matrix[K] Sig_K;
  cov_matrix[N] Sig_N;

  mu_K <- rep_vector(0.0, K);  
  Sig_K <- diag_matrix(rep_vector(1.0, K));
  Sig_N <- diag_matrix(rep_vector(1.0, N));  
}

transformed parameters{
  matrix[N, P] Theta;             //parameter matrix
  Theta <- U * V' + e; 
}

model {
  // the priors 
  for(i in 1:N){
    e[i] ~ multi_normal(mu_K, Sig_K); 
    //U[i] ~ multi_normal(mu_K, Sig_K);  
  }
  for(i in 1:P){
    V[i] ~ multi_normal(mu_K, Sig_K);   
  }

  //The likelihood
  for(j in 1:N){
    Y[j] ~ multi_normal(Theta, Sig_N); 
  }
}