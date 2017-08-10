functions {
  vector cum_prod(vector vec_input){

    return exp(cumulative_sum(log(vec_input)));

  }
}

data {
  int<lower=1> N;               // number of population
  int<lower=1> P;               // number of species
  int<lower=1> Q;               // number of covariates    
  int<lower=1> K;               // dimension of latent factors 
  matrix[N,P] Y;                // data matrix of order [N,P]  
  matrix[N,Q] X;                // covar matrix of order [N,Q]    
  
  real<lower=0> lambda_v;             // observation noise precision
  real<lower=0> lambda_u;             // penalty (inverse variance) for source emission uncertainty

  matrix[K,P] p_V;          // source profile prior
  real<lower=0> v_p;        // hyperparameter for gamma_V   
  real<lower=0> v_k;        // hyperparameter for gamma_rank    
}

transformed data {
  // helper parameters for constructing prior
  vector[Q] mu_Q;    
  vector[K] mu_K;      
  cov_matrix[K] Sig_K;  
  cov_matrix[P] Sig_P;  
  cov_matrix[Q] Sig_Q;  

  mu_Q = rep_vector(0.0, Q);
  mu_K = rep_vector(0.0, K);
  Sig_K = diag_matrix(rep_vector(1/lambda_u, K));   
  Sig_P = diag_matrix(rep_vector(1/lambda_v, P));
  Sig_Q = diag_matrix(rep_vector(1, Q)); 
}

parameters {
  vector[Q] B;        	   		  // covariate coeffcient   
  real<lower=0> sigma_b;         // regularization for covariate coeffcient   

  simplex[P] V_spx[K];                // latent factor for source profile (matrix)        
  matrix<lower=0>[N, K] U;            // factor loading for population/source emission  
  vector<lower=0>[K] gamma_rank;      // penalty factor control matrix rank  
  //vector<lower=0>[K] U_unc;         // uncertainty for source emission
}

transformed parameters{
  matrix<lower=0,upper=1>[P, K] V;   // latent factor for species/source profile (matrix)      
  vector[K] gamma_rank_cumprod;      // cumulative product of gamma_rank parameter  
  cov_matrix[K] D;                   // shrinkage matrix
  matrix[N, P] Theta;                // parameter matrix

  gamma_rank_cumprod = cum_prod(gamma_rank);  // cumulative product of gamma rank
  D = diag_matrix(gamma_rank_cumprod);  	  // cumulative product of gamma rank

  for (k in 1:K){
    V[, k] = V_spx[k];   
  }

  Theta = U * D * V';
  // Theta = X * B + U * D * V';
}

model {
  // the covariate coefficients 
  B ~ multi_normal(mu_Q, sigma_b * Sig_Q);
  sigma_b ~ inv_gamma(1, 1);
  

  // the priors 
  for (k in 1:K){
    //U_unc[k] ~ exponential(lambda_u);
    gamma_rank[k] ~ inv_gamma(v_k, 1);        // D matrix
    V_spx[k] ~ dirichlet(to_vector(p_V[k]));  // V matrix   
  }

  for (i in 1:N){
    log(U[i]) ~ multi_normal(mu_K, Sig_K);    // U matrix
  }
  target += -sum(log(U));

  //The likelihood
  for(j in 1:N){
    Y[j] ~ multi_normal(X[j] * B + log(Theta[j]), Sig_P);   
  }
}

