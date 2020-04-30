data {
  int<lower=1> N; // number of treatments
  int<lower=1> J; // number of tree species
  
  int<lower=1,upper=J> tid[N]; // tree indices
  
  int<lower=0> x[N];     // dose in ob
  int<lower=0> y[N];     // dependent variable, number of virus-killed
  int<lower=0> total[N]; // total caterpillars in treatment
  
  vector[J] prior_mu_alpha;
  vector[J] prior_mu_beta;
  
  real<lower=0> prior_sigma_alpha;
  real<lower=0> prior_sigma_beta;
}

parameters{
  vector[J] raw_alpha;
  vector[J] raw_beta;
  
  real<lower=0> sigma_alpha;
  real<lower=0> sigma_beta;
}

transformed parameters {
  vector[J] alpha;
  vector[J] beta;
  
  vector[N] theta; //predicted proportion virus-killed on logit scale
  vector[N] inv_logit_theta;
  
  alpha = prior_mu_alpha + raw_alpha * sigma_alpha;
  beta = prior_mu_beta + raw_beta * sigma_beta;

  for(n in 1:N) {
    theta[n] = alpha[tid[n]] + beta[tid[n]] * x[n];
  }
  
  inv_logit_theta = inv_logit(theta);
}

model {
  //priors
  sigma_alpha ~ normal(prior_sigma_alpha, prior_sigma_alpha*3);
  sigma_beta ~ normal(prior_sigma_beta, prior_sigma_beta*3);
  
  raw_alpha ~ normal(0,1);
  raw_beta ~ normal(0,1);
  

  //likelihood
  y ~ binomial_logit(total, theta);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = binomial_logit_lpmf(y[n] | total[n], theta[n]);
  }
}
