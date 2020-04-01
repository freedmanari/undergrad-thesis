data {
  int<lower=1> N; // number of treatments
  int<lower=1> I; // number of tree species
  int<lower=1> J; // number of strains
  int<lower=1> K; // number of capsids
  
  int<lower=1,upper=I> tid[N]; // tree indices
  int<lower=1,upper=J> sid[N]; // strain indices
  int<lower=1,upper=K> cid[J]; // capsid indices
  
  int<lower=0> x[N];     // dose in ob/uL
  int<lower=0> y[N];     // dependent variable, number of virus-killed
  int<lower=0> total[N]; // total caterpillars in treatment
  
  matrix[J,I] prior_mu_alpha;
  matrix[J,I] prior_mu_beta;
  
  vector<lower=0>[K] prior_sigma_alpha;
  vector<lower=0>[K] prior_sigma_beta;
}

parameters{
  matrix[J,I] raw_alpha;
  matrix[J,I] raw_beta;
  
  vector<lower=0>[K] sigma_alpha;
  vector<lower=0>[K] sigma_beta;
}

transformed parameters {
  matrix[J,I] alpha;
  matrix[J,I] beta;
  
  vector[N] theta; //predicted proportion virus-killed on logit scale
  vector[N] inv_logit_theta;
  
  for(j in 1:J) {
      alpha[j] = prior_mu_alpha[j] + raw_alpha[j] * sigma_alpha[cid[j]];
      beta[j] = prior_mu_beta[j] + raw_beta[j] * sigma_beta[cid[j]];
  }
  
  for(n in 1:N) {
    theta[n] = alpha[sid[n],tid[n]] + beta[sid[n],tid[n]] * x[n];
  }
  
  inv_logit_theta = inv_logit(theta);
}

model {
  //priors
  sigma_alpha ~ normal(prior_sigma_alpha, 3*prior_sigma_alpha);
  sigma_beta ~ normal(prior_sigma_beta, 3*prior_sigma_beta);
  
  for (j in 1:J) {
    raw_alpha[j] ~ normal(0,1);
    raw_beta[j] ~ normal(0,1);
  }
  

  //likelihood
  y ~ binomial_logit(total, theta);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_logit_lpmf(y[n] | total[n], theta[n]);
  }
}
