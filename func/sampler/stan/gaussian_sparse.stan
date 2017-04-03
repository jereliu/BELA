functions {
  vector cum_prod(vector vec_input){

    return exp(cumulative_sum(log(vec_input)));

  }
}


data {
  int<lower=1> N;               // number of population
  int<lower=1> P;               // number of species
  int<lower=1> K;               // dimension of latent factors 
  matrix[N,P] Y;                // data matrix of order [N,P]    
  
  real<lower=0> lambda_u;       // penalty for prior
  real<lower=0> lambda_v;  

  real<lower=0> v_p;        // hyperparameter for gamma_V   
  real<lower=0> v_k;  			// hyperparameter for gamma_rank 	
}

transformed data {
  // helper parameters for constructing prior
  vector[K] mu_K;    
  cov_matrix[K] Sig_U;
  cov_matrix[P] Sig_Y;

  mu_K = rep_vector(0.0, K);
  Sig_U = diag_matrix(rep_vector(1/lambda_u, K));
  Sig_Y = diag_matrix(rep_vector(1, P));
}

parameters {
  matrix[N, K] U;   // factor loading for population
  matrix[P, K] V;   // latent factor for species       
  matrix<lower=0>[P, K] gamma_V; // penalty factor control magnitude of V
  vector<lower=0>[K] gamma_rank; // penalty factor control matrix rank

  //matrix[N, P] e;     // factor noise  
}

transformed parameters{
  vector[K] gamma_rank_cumprod;		// cumulative product of gamma_rank parameter
  matrix[N, P] Theta;               // parameter matrix  

  gamma_rank_cumprod = cum_prod(gamma_rank);  // cumulative product of gamma rank
  Theta = U * V'; 
  //Theta = U * V' + e; 
}

model {
  // the priors 
  for(i in 1:N){
    //e[i] ~ multi_normal(mu_K, Sig_U); 
    U[i] ~ multi_normal(mu_K, Sig_U);  
  }
  
  for (j in 1:K){
  	gamma_rank[j] ~ inv_gamma(v_k, 1);
  	for(i in 1:P){
  	 	gamma_V[i, j] ~ inv_gamma(v_p/2, v_p/2);
    	V[i,j] ~ normal(0, gamma_V[i, j] * gamma_rank_cumprod[j]);   
    }
  }

  //The likelihood
  for(j in 1:N){
    Y[j] ~ multi_normal(Theta[j], Sig_Y); 
  }
}
