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
  
  vector[N] theta; //predicted proportion virus-killed
  
  for(j in 1:J) {
    alpha[j] = raw_alpha[j] * sigma_alpha[cid[j]];
    beta[j] = raw_beta[j] * sigma_beta[cid[j]];
  }
  
  for(n in 1:N) {
    theta[n] = alpha[sid[n],tid[n]] + beta[sid[n],tid[n]]*x[n];
  }
}

model {
  //priors
  sigma_alpha ~ cauchy(0,2.5);
  sigma_beta ~ cauchy(0,2.5);
  
  for (j in 1:J) {
    raw_alpha[j] ~ normal(0,1);
    raw_beta[j] ~ normal(0,1);
  }
  
  
  //likelihood
  y ~ binomial_logit(total, theta);
}

