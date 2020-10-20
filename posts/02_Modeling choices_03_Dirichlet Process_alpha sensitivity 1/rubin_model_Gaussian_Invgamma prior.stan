data {
  int<lower = 0> K;         // number of sites 
  real tau_hat_k[K];        // estimated treatment effects
  real<lower=0> se_k[K];    // s.e. of effect estimates 
}
parameters {
  real tau; 
  real<lower=0> sigma_tau2;
  real tau_k[K];
}
transformed parameters {
  real<lower=0> sigma_tau;
  sigma_tau = sqrt(sigma_tau2);
}
model {
  sigma_tau2 ~ inv_gamma(1, 1);     // vague prior
  tau ~ normal(0, 100);             // vague prior
  tau_k ~ normal(tau, sigma_tau);   // second-level normal
  tau_hat_k ~ normal(tau_k, se_k);  // third-level normal
}
generated quantities{
  real predicted_tau_k;
  predicted_tau_k = normal_rng(tau, sigma_tau);  //mixed predictive 
}
