data {
  int<lower=1> N;         // number of caterpillars
  int<lower=1> J;         // number of strains
  int<lower=1,upper=J> K; // number of capsids
  
  int<lower=1,upper=J> sid[N]; // strain indices
  int<lower=1,upper=K> cid[J]; // capsid indices
  
  int<lower=0> x_obs[N];          // dose in ob/uL  
  int<lower=0,upper=1> x_tree[N]; // tree_sp
  int<lower=0,upper=1> y[N];      // dependent variables 
}

parameters{
  vector[J] raw_alpha;
  vector[J] raw_beta_obs;
  vector[J] raw_beta_tree;
  
  vector[J] mu_alpha;
  vector[J] mu_obs;
  vector[J] mu_tree;
  
  vector<lower=0>[K] sigma_alpha;
  vector<lower=0>[K] sigma_obs;
  vector<lower=0>[K] sigma_tree;
}

transformed parameters {
  vector[J] alpha;
  vector[J] beta_obs;
  vector[J] beta_tree;
  
  vector[N] y_hat; //predicted y
  
  for(j in 1:J) {
    alpha[j] = mu_alpha[j] + sigma_alpha[cid[j]]*raw_alpha[j];
    beta_obs[j] = mu_obs[j] + sigma_obs[cid[j]]*raw_beta_obs[j];
    beta_tree[j] = mu_tree[j] + sigma_tree[cid[j]]*raw_beta_tree[j];
  }
  
  for(n in 1:N) {
    y_hat[n] = alpha[sid[n]] + beta_obs[sid[n]]*x_obs[n] + beta_tree[sid[n]]*x_tree[n];
  }
}

model {
  //priors
  mu_alpha ~ normal(0,5);
  mu_obs ~ normal(0,5);
  mu_tree ~ normal(0,5);
  
  sigma_alpha ~ cauchy(0,2.5);
  sigma_obs ~ cauchy(0,2.5);
  sigma_tree ~ cauchy(0,2.5);

  raw_alpha ~ normal(0,1);
  raw_beta_obs ~ normal(0,1);
  raw_beta_tree ~ normal(0,1);
  
  //likelihood
  y ~ bernoulli_logit(y_hat);
}

