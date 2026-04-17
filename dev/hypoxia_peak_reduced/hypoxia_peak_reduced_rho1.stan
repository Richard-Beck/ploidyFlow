data {
  int<lower=1> N;
  int<lower=1> N_date;
  array[N] int<lower=1, upper=N_date> date_id;
  vector[N] x_ratio;
  vector<lower=0>[N] y_cen;
  vector<lower=0>[N] y_g1;
}

parameters {
  real alpha;
  real<lower=0> beta;
  real log_M_cen;
  real<lower=log(1.25), upper=log(30)> log_R_2n;
  real<lower=log(1.1), upper=log(4)> log_R_4n_over_2n;
  vector[N_date] log_phi_raw;
  real<lower=0.01, upper=1> sigma_cen;
  real<lower=0.01, upper=1> sigma_g1;
  real<lower=0.01, upper=0.99> p_4n;
}

transformed parameters {
  real M_cen = exp(log_M_cen);
  real R_2n = exp(log_R_2n);
  real R_4n = R_2n * exp(log_R_4n_over_2n);
  real rho = 1;
  vector[N_date] log_phi = log_phi_raw - mean(log_phi_raw);
  vector[N_date] phi_date = exp(log_phi);
  vector[N] u;
  vector[N] mu_cen;
  vector[N] mu_g1_2n;
  vector[N] mu_g1_4n;

  for (n in 1:N) {
    real ell = alpha + beta * x_ratio[n];
    real distortion;
    u[n] = exp(-ell);
    mu_cen[n] = M_cen * u[n] / (1 + u[n]);
    distortion = (1 + u[n]) / (1 + rho * u[n]);
    mu_g1_2n[n] = mu_cen[n] * phi_date[date_id[n]] * R_2n * distortion;
    mu_g1_4n[n] = mu_cen[n] * phi_date[date_id[n]] * R_4n * distortion;
  }
}

model {
  alpha ~ normal(0, 1.5);
  beta ~ normal(0.6, 0.5);
  log_M_cen ~ normal(log(5000), 0.75);
  log_R_2n ~ normal(log(8), 0.35);
  log_R_4n_over_2n ~ normal(log(2), 0.2);
  log_phi_raw ~ normal(0, 0.35);
  sigma_cen ~ normal(0, 0.15);
  sigma_g1 ~ normal(0, 0.15);
  p_4n ~ beta(1.5, 1.5);

  y_cen ~ lognormal(log(mu_cen), sigma_cen);
  for (n in 1:N) {
    target += log_mix(
      p_4n,
      lognormal_lpdf(y_g1[n] | log(mu_g1_4n[n]), sigma_g1),
      lognormal_lpdf(y_g1[n] | log(mu_g1_2n[n]), sigma_g1)
    );
  }
}

generated quantities {
  vector[N] prob_4n;
  vector[N] mu_g1_mixture;

  for (n in 1:N) {
    real log_w4 = log(p_4n) + lognormal_lpdf(y_g1[n] | log(mu_g1_4n[n]), sigma_g1);
    real log_w2 = log1m(p_4n) + lognormal_lpdf(y_g1[n] | log(mu_g1_2n[n]), sigma_g1);
    prob_4n[n] = exp(log_w4 - log_sum_exp(log_w4, log_w2));
    mu_g1_mixture[n] = prob_4n[n] * mu_g1_4n[n] + (1 - prob_4n[n]) * mu_g1_2n[n];
  }
}
