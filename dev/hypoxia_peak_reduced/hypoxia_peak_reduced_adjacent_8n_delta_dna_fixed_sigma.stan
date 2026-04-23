data {
  int<lower=1> N;
  int<lower=1> N_date;
  array[N] int<lower=1, upper=N_date> date_id;
  vector[N] x_ratio;
  vector<lower=0>[N] y_cen;
  vector<lower=0>[N] y_peak_lower;
  array[N] int<lower=0, upper=1> has_upper_peak;
  vector<lower=0>[N] y_peak_upper;
  real<lower=0.001, upper=0.5> sigma_delta_dna_fixed;
}

parameters {
  real alpha;
  real<lower=0> beta;
  real log_M_cen;
  real<lower=log(1.25), upper=log(30)> log_R_2n;
  real<lower=log(1.1), upper=log(4)> log_R_4n_over_2n;
  real<lower=log(1.1), upper=log(4)> log_R_8n_over_4n;
  real log_rho;
  vector[N_date] log_phi_raw;
  vector[N] delta_dna;
  real<lower=0.01, upper=1> sigma_cen;
  real<lower=0.01, upper=1> sigma_g1;
  real<lower=0.01, upper=0.99> p_4n;
}

transformed parameters {
  real M_cen = exp(log_M_cen);
  real R_2n = exp(log_R_2n);
  real R_4n = R_2n * exp(log_R_4n_over_2n);
  real R_8n = R_4n * exp(log_R_8n_over_4n);
  real rho = exp(log_rho);
  vector[N_date] log_phi = log_phi_raw - mean(log_phi_raw);
  vector[N_date] phi_date = exp(log_phi);
  vector[N] u;
  vector[N] delta_dna_multiplier;
  vector[N] mu_cen;
  vector[N] mu_peak_2n;
  vector[N] mu_peak_4n;
  vector[N] mu_peak_8n;

  for (n in 1:N) {
    real ell = alpha + beta * x_ratio[n];
    real distortion;
    u[n] = exp(-ell);
    delta_dna_multiplier[n] = exp(delta_dna[n]);
    mu_cen[n] = M_cen * u[n] / (1 + u[n]);
    distortion = (1 + u[n]) / (1 + rho * u[n]);
    mu_peak_2n[n] = mu_cen[n] * phi_date[date_id[n]] * R_2n * distortion * delta_dna_multiplier[n];
    mu_peak_4n[n] = mu_cen[n] * phi_date[date_id[n]] * R_4n * distortion * delta_dna_multiplier[n];
    mu_peak_8n[n] = mu_cen[n] * phi_date[date_id[n]] * R_8n * distortion * delta_dna_multiplier[n];
  }
}

model {
  alpha ~ normal(0, 1.5);
  beta ~ normal(0.6, 0.5);
  log_M_cen ~ normal(log(5000), 0.75);
  log_R_2n ~ normal(log(8), 0.35);
  log_R_4n_over_2n ~ normal(log(2), 0.2);
  log_R_8n_over_4n ~ normal(log(2), 0.2);
  log_rho ~ normal(0, 0.5);
  log_phi_raw ~ normal(0, 0.35);
  delta_dna ~ normal(0, sigma_delta_dna_fixed);
  sigma_cen ~ normal(0, 0.15);
  sigma_g1 ~ normal(0, 0.15);
  p_4n ~ beta(1.5, 1.5);

  y_cen ~ lognormal(log(mu_cen), sigma_cen);
  for (n in 1:N) {
    real lower_if_2n = log1m(p_4n) + lognormal_lpdf(y_peak_lower[n] | log(mu_peak_2n[n]), sigma_g1);
    real lower_if_4n = log(p_4n) + lognormal_lpdf(y_peak_lower[n] | log(mu_peak_4n[n]), sigma_g1);

    if (has_upper_peak[n] == 1) {
      lower_if_2n += lognormal_lpdf(y_peak_upper[n] | log(mu_peak_4n[n]), sigma_g1);
      lower_if_4n += lognormal_lpdf(y_peak_upper[n] | log(mu_peak_8n[n]), sigma_g1);
    }

    target += log_sum_exp(lower_if_2n, lower_if_4n);
  }
}

generated quantities {
  vector[N] prob_4n;
  vector[N] mu_lower_mixture;
  vector[N] mu_upper_mixture;

  for (n in 1:N) {
    real log_w2 = log1m(p_4n) + lognormal_lpdf(y_peak_lower[n] | log(mu_peak_2n[n]), sigma_g1);
    real log_w4 = log(p_4n) + lognormal_lpdf(y_peak_lower[n] | log(mu_peak_4n[n]), sigma_g1);
    if (has_upper_peak[n] == 1) {
      log_w2 += lognormal_lpdf(y_peak_upper[n] | log(mu_peak_4n[n]), sigma_g1);
      log_w4 += lognormal_lpdf(y_peak_upper[n] | log(mu_peak_8n[n]), sigma_g1);
    }

    prob_4n[n] = exp(log_w4 - log_sum_exp(log_w4, log_w2));
    mu_lower_mixture[n] = prob_4n[n] * mu_peak_4n[n] + (1 - prob_4n[n]) * mu_peak_2n[n];
    mu_upper_mixture[n] = prob_4n[n] * mu_peak_8n[n] + (1 - prob_4n[n]) * mu_peak_4n[n];
  }
}
