data {
  // N is the number of histogram bins after stacking all samples together.
  // Each row i represents one bin center and its observed event/cell count.
  int<lower=1> N;

  // N_exp is the number of experiments/batches. The centering/control peak is
  // estimated at this level because the DNA-area scale can shift by experiment.
  int<lower=1> N_exp;

  // N_sample is the number of biological/flow samples. Mixture weights and
  // sample-specific ploidy peak positions are estimated at this level.
  int<lower=1> N_sample;

  // Number of bins per sample histogram. This is passed for bookkeeping by the
  // R data-preparation code; the likelihood below uses the stacked length N.
  int<lower=1> N_bin;

  // Retained for compatibility with existing prepared Stan data. The current
  // simplified S-phase likelihood does not use subpeaks.
  int<lower=1> N_s_phase;

  // Center of each DNA-area histogram bin and the observed count in that bin.
  vector[N] areaDNA_bin_center;
  array[N] int<lower=0> count;

  // sampleExpID maps sample s -> experiment e. expID and sampleID map stacked
  // histogram row i -> experiment e and sample s.
  array[N_sample] int<lower=1, upper=N_exp> sampleExpID;
  array[N] int<lower=1, upper=N_exp> expID;
  array[N] int<lower=1, upper=N_sample> sampleID;
}

transformed data {
  // The observed DNA-area range defines the parameter bounds for the
  // centering peak.
  real min_x = min(areaDNA_bin_center);
  real max_x = max(areaDNA_bin_center);

  // Small separation added to the tumor/control ratio. This keeps the 2n peak
  // width positive even if exp(.) underflows during an optimizer line search.
  real min_tumor_cen_ratio_gap = 1e-3;

  // Minimum extra spacing above 1 for ploidy-to-ploidy ratios. A 0.5 gap means
  // 4n/2n and 8n/4n cannot collapse below 1.5, while the prior below still
  // centers them near 2.
  real min_ploidy_ratio_gap = 0.5;
}

parameters {
  // Experiment-level location of the low-DNA centering/control peak. The upper
  // bound keeps this peak in the lower third of the observed DNA-area range so
  // it cannot swap labels with the tumor ploidy peaks.
  vector<lower=min_x, upper=max_x / 3>[N_exp] mu_cen;

  // Global mean log ratio from the centering/control peak to the 2n tumor peak.
  real log_ratio_tumor_cen_mu;

  // Global log increments for the 4n/2n and 8n/4n peak ratios. These ratios
  // are shared across samples to keep the first-pass model stable.
  real log_R_4n_mu;
  real log_R_8n_mu;

  // Shared coefficient of variation for all Gaussian peak widths. A peak with
  // mean mu has standard deviation mu * cv.
  real<lower=0.01, upper=0.5> cv;

  // Sample-specific mixture weights for three latent populations: centering
  // chicken/control, diploid tumor, and WGD tumor.
  array[N_sample] simplex[3] pop_theta;

  // Global allocation of the centering/control population into singlet and
  // doublet noise peaks.
  simplex[2] cen_phase;

  // Global cell-cycle allocations. Entries are G0/G1, S phase, and G2/M.
  // Diploid G2/M and WGD G0/G1 both contribute to the observed 4n peak.
  simplex[3] diploid_phase;
  simplex[3] wgd_phase;
}

transformed parameters {
  // Sample-specific derived peak ratios, locations, and standard deviations.
  // Ratios are global and repeated here for plotting/output convenience.
  vector[N_sample] ratio_tumor_cen_sample;
  vector[N_sample] R_4n_sample;
  vector[N_sample] R_8n_sample;
  vector[N_sample] mu_2n_sample;
  vector[N_sample] mu_4n_sample;
  vector[N_sample] mu_8n_sample;
  vector[N_sample] sig_cen_sample;
  vector[N_sample] sig_2n_sample;
  vector[N_sample] sig_4n_sample;
  vector[N_sample] sig_8n_sample;
  array[N_sample] vector[7] theta;

  for (s in 1:N_sample) {
    int e = sampleExpID[s];

    // Convert global log-ratios into positive multiplicative ratios.
    ratio_tumor_cen_sample[s] = 1 + min_tumor_cen_ratio_gap + exp(log_ratio_tumor_cen_mu);
    R_4n_sample[s] = 1 + min_ploidy_ratio_gap + exp(log_R_4n_mu);
    R_8n_sample[s] = 1 + min_ploidy_ratio_gap + exp(log_R_8n_mu);

    // Build the ordered ploidy peak locations. The 2n peak is a multiplicative
    // offset from the experiment's centering peak; 4n and 8n are multiplicative
    // offsets from the previous tumor peak.
    sig_cen_sample[s] = mu_cen[e] * cv;
    mu_2n_sample[s] = mu_cen[e] * ratio_tumor_cen_sample[s];
    mu_4n_sample[s] = mu_2n_sample[s] * R_4n_sample[s];
    mu_8n_sample[s] = mu_4n_sample[s] * R_8n_sample[s];

    // Peak widths scale with peak location through the common CV.
    sig_2n_sample[s] = mu_2n_sample[s] * cv;
    sig_4n_sample[s] = mu_4n_sample[s] * cv;
    sig_8n_sample[s] = mu_8n_sample[s] * cv;

    // Deterministic observed-component masses from the latent population model:
    // 1 CEN singlet, 2 diploid G0/G1 at 2n, 3 combined 4n peak, 4 WGD G2/M at
    // 8n, 5 diploid S phase, 6 WGD S phase, 7 CEN doublet.
    theta[s][1] = pop_theta[s][1] * cen_phase[1];
    theta[s][2] = pop_theta[s][2] * diploid_phase[1];
    theta[s][3] = pop_theta[s][2] * diploid_phase[3] + pop_theta[s][3] * wgd_phase[1];
    theta[s][4] = pop_theta[s][3] * wgd_phase[3];
    theta[s][5] = pop_theta[s][2] * diploid_phase[2];
    theta[s][6] = pop_theta[s][3] * wgd_phase[2];
    theta[s][7] = pop_theta[s][1] * cen_phase[2];
  }
}

model {
  // Priors for the global tumor/control ratio and the population
  // mean ploidy-to-ploidy ratios. log(9) says the 2n tumor peak is expected to
  // sit about 9x above the centering peak. log_R priors near log(0.5) imply
  // 1 + min_ploidy_ratio_gap + exp(log(0.5)) = 2, so 4n is expected near 2x
  // 2n and 8n near 2x 4n.
  log_ratio_tumor_cen_mu ~ normal(log(8), 0.35);
  log_R_4n_mu ~ normal(log(0.5), 0.05);
  log_R_8n_mu ~ normal(log(0.5), 0.05);

  for (s in 1:N_sample) {
    pop_theta[s] ~ dirichlet([2, 4, 2]');
  }
  cen_phase ~ dirichlet([8, 0.5]');
  diploid_phase ~ dirichlet([4, 0.75, 1]');
  wgd_phase ~ dirichlet([4, 0.75, 1]');
  cv ~ normal(0.05, 0.05);

  // Likelihood: each observed histogram bin is modeled as repeated draws from
  // a sample-specific seven-group mixture density evaluated at the bin center.
  // The two S-phase groups are simple broad Gaussians centered between their
  // adjacent ploidy peaks. Multiplying by count[i] is equivalent to adding that
  // log density once for each event/cell that fell into the bin.
  for (i in 1:N) {
    if (count[i] > 0) {
      int s = sampleID[i];
      real x = areaDNA_bin_center[i];
      vector[7] comp;

      // Each comp[k] is log mixture weight + log component density at x.
      comp[1] = log(theta[s][1]) + normal_lpdf(x | mu_cen[sampleExpID[s]], sig_cen_sample[s]);
      comp[2] = log(theta[s][2]) + normal_lpdf(x | mu_2n_sample[s], sig_2n_sample[s]);
      comp[3] = log(theta[s][3]) + normal_lpdf(x | mu_4n_sample[s], sig_4n_sample[s]);
      comp[4] = log(theta[s][4]) + normal_lpdf(x | mu_8n_sample[s], sig_8n_sample[s]);
      comp[5] = log(theta[s][5]) + normal_lpdf(
        x |
        0.5 * (mu_2n_sample[s] + mu_4n_sample[s]),
        0.5 * (mu_4n_sample[s] - mu_2n_sample[s])
      );
      comp[6] = log(theta[s][6]) + normal_lpdf(
        x |
        0.5 * (mu_4n_sample[s] + mu_8n_sample[s]),
        0.5 * (mu_8n_sample[s] - mu_4n_sample[s])
      );
      comp[7] = log(theta[s][7]) + normal_lpdf(
        x |
        2 * mu_cen[sampleExpID[s]],
        sqrt(2) * sig_cen_sample[s]
      );

      // Marginalize over the unknown component label for this bin center.
      target += count[i] * log_sum_exp(comp);
    }
  }
}
