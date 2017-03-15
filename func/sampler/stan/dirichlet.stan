data {
  int<lower=1> N;               // number of population
  int<lower=1> P;               // number of species
  int<lower=1> K;               // dimension of latent factors 
  matrix[N,P] Y;                // data matrix of order [N,P]    
  
  real<lower=0> lambda_v;             // observation noise precision
  real<lower=0> lambda_u;             // penalty (inverse variance) for source emission uncertainty

  matrix[K,P] p_V;                  // source profile prior
}

transformed data {
  // helper parameters for constructing prior
  vector[P] mu_P;    
  vector[K] mu_K;      
  cov_matrix[K] Sig_K;  
  cov_matrix[P] Sig_P;  

  mu_K = rep_vector(0.0, K);
  mu_P = rep_vector(0.0, P);  
  Sig_K = diag_matrix(rep_vector(1/lambda_u, K));   
  Sig_P = diag_matrix(rep_vector(1/lambda_v, P));
}

parameters {
  simplex[P] V_spx[K];                // latent factor for species/source profile (matrix)        
  matrix<lower=0>[N, K] U;            // factor loading for population/source emission  
  //vector<lower=0>[K] U_unc;           // uncertainty for source emission
}

transformed parameters{
  matrix<lower=0,upper=1>[P, K] V;   // latent factor for species/source profile (matrix)      
  matrix[N, P] Theta;                // parameter matrix

  //Sig_K = diag_matrix(U_unc);  

  for (k in 1:K){
    V[, k] = V_spx[k];   
  }

  Theta = U * V';
}

model {
  // the priors 
  for (k in 1:K){
    //U_unc[k] ~ exponential(lambda_u);
    V_spx[k] ~ dirichlet(to_vector(p_V[k]));   
  }

  for (i in 1:N){
    log(U[i]) ~ multi_normal(mu_K, Sig_K);
  }
  target += -sum(log(U));

  //The likelihood
  for(j in 1:N){
    Y[j] ~ multi_normal(log(Theta[j]), Sig_P); 
  }
}
