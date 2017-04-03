functions {
  vector cum_prod(vector vec_input){

    return exp(cumulative_sum(log(vec_input)));

  }
}


data {
  int<lower=1> N;               // number of population
  int<lower=1> P;               // number of species
  int<lower=1> K;               // dimension of latent factors 
  int<lower=0> Y[N,P];          // data matrix of order [N,P]    
  
  real<lower=0> lambda_u;       // penalty for prior
  real<lower=0> lambda_v;  

  real<lower=0> v_k;  			// hyperparameter for: gamma_rank ~ Gamma(v_k, 1)
  real<lower=0> v_p;			  // hyperparameter for: gamma_V    ~ Gamma(v_p/2, v_p/2)
}

transformed data {
  // helper parameters for constructing prior
  vector[K] mu_K;    
  cov_matrix[K] Sig_K;
  // cov_matrix[P] Sig_P;

  mu_K = rep_vector(0.0, K);
  Sig_K = diag_matrix(rep_vector(1/lambda_u, K));
  // Sig_P = diag_matrix(rep_vector(1, P));
}

parameters {
  matrix[N, K] U;   // factor loading for population
  matrix[P, K] V;   // latent factor for species       
  matrix<lower=0>[P, K] gamma_V; // penalty factor control magnitude of V
  vector<lower=0>[K] gamma_rank; // penalty factor control matrix rank

  //matrix[N, P] e;     // factor noise  
}

transformed parameters{
  vector[K] gamma_rank_cumprod;   // cumulative product of gamma_rank parameter
  matrix[N, P] Theta;             	// parameter matrix

  gamma_rank_cumprod = cum_prod(gamma_rank);  // cumulative product of gamma rank
  Theta = U * V'; 
  //Theta = U * V' + e; 
}

model {
  // the priors 
  for(i in 1:N){
    //e[i] ~ multi_normal(mu_K, Sig_K); 
    U[i] ~ multi_normal(mu_K, Sig_K);  
  }
  
  for (j in 1:K){
  	gamma_rank[j] ~ inv_gamma(v_k, 1);
  	for(i in 1:P){
  	 	gamma_V[i, j] ~ inv_gamma(v_p/2, v_p/2); // remove element-wise variance for now
    	V[i,j] ~ normal(0, gamma_V[i, j] * gamma_rank_cumprod[j]); 
    }
  }

  //The likelihood
  for(i in 1:N){
    for (j in 1:P){
      Y[i, j] ~  poisson_log(Theta[i, j]); 
   }
  }
}
