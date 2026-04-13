functions {
  real solve_reduced_mechanistic_u(real D_tilde, real A_C_tilde, real A_T_tilde, real rho, real W_s) {
    real u = D_tilde / (1 + A_C_tilde + A_T_tilde * W_s * rho);

    for (iter in 1:50) {
      real f = u +
        A_C_tilde * u / (1 + u) +
        A_T_tilde * W_s * (rho * u) / (1 + rho * u) -
        D_tilde;
      real fp = 1 +
        A_C_tilde / square(1 + u) +
        A_T_tilde * W_s * rho / square(1 + rho * u);
      real u_next = u - f / fp;

      if (!(u_next > 0)) {
        u_next = 0.5 * u;
      }
      if (abs(u_next - u) < 1e-10 * fmax(1, u)) {
        return u_next;
      }
      u = u_next;
    }

    return u;
  }
}

data {
  // N is the number of histogram bins after stacking all samples together.
  // Each row i represents one bin center and its observed event/cell count.
  int<lower=1> N;

  // N_exp is the number of experiments/batches. The baseline centering/control
  // peak is estimated at this level because the DNA-area scale can shift by
  // experiment.
  int<lower=1> N_exp;

  // N_sample is the number of biological/flow samples. Mixture weights and
  // sample-specific peak locations are estimated at this level.
  int<lower=1> N_sample;

  // N_group is the number of sample groups derived from the FACS filenames
  // (for example AOA, AN, N1). sampleGroupID maps sample s -> group g.
  int<lower=1> N_group;

  // Number of bins per sample histogram. This is passed for bookkeeping by the
  // R data-preparation code; the likelihood below uses the stacked length N.
  int<lower=1> N_bin;

  // Retained for compatibility with existing prepared Stan data. The current
  // simplified S-phase likelihood does not use subpeaks.
  int<lower=1> N_s_phase;

  // Fixed DNA-area threshold used to estimate each sample's observed CEN
  // fraction from the histogram itself. A value near 20000 is a good default.
  real<lower=0> cen_threshold;

  // Center of each DNA-area histogram bin and the observed count in that bin.
  vector[N] areaDNA_bin_center;
  array[N] int<lower=0> count;

  // sampleExpID maps sample s -> experiment e. expID and sampleID map stacked
  // histogram row i -> experiment e and sample s.
  array[N_sample] int<lower=1, upper=N_exp> sampleExpID;
  array[N_sample] int<lower=1, upper=N_group> sampleGroupID;
  array[N] int<lower=1, upper=N_exp> expID;
  array[N] int<lower=1, upper=N_sample> sampleID;
}

transformed data {
  // The observed DNA-area range defines the parameter bounds for the baseline
  // centering peak.
  real min_x = min(areaDNA_bin_center);
  real max_x = max(areaDNA_bin_center);

  // Small separation added to the tumor/control ratio. This keeps the 2n peak
  // width positive even if exp(.) underflows during an optimizer line search.
  real min_tumor_cen_ratio_gap = 1e-3;

  // Minimum extra spacing above 1 for ploidy-to-ploidy ratios. A 0.5 gap means
  // 8n/4n cannot collapse below 1.5, while the prior below still centers it
  // near 2.
  real min_ploidy_ratio_gap = 0.5;

  // Hard bounds for the 4n/2n peak ratio.
  real R_4n_lower = 1.9;
  real R_4n_upper = 2.05;
  real log_R_4n_lower = log(R_4n_lower - 1 - min_ploidy_ratio_gap);
  real log_R_4n_upper = log(R_4n_upper - 1 - min_ploidy_ratio_gap);

  // Histogram-derived estimate of the observed CEN fraction in each sample.
  array[N_sample] int<lower=0> sample_total_count;
  array[N_sample] int<lower=0> sample_cen_count;
  vector[N_sample] observed_cen_fraction;

  for (s in 1:N_sample) {
    sample_total_count[s] = 0;
    sample_cen_count[s] = 0;
  }

  for (i in 1:N) {
    int s = sampleID[i];
    sample_total_count[s] += count[i];
    if (areaDNA_bin_center[i] < cen_threshold) {
      sample_cen_count[s] += count[i];
    }
  }

  for (s in 1:N_sample) {
    if (sample_total_count[s] > 0) {
      observed_cen_fraction[s] = sample_cen_count[s] / (1.0 * sample_total_count[s]);
    } else {
      observed_cen_fraction[s] = 0;
    }
  }
}

parameters {
  // Experiment-level asymptotic high-dye location of the low-DNA
  // centering/control peak. The sample-specific free-dye state rescales this
  // down to the observed CEN location through u_s / (1 + u_s).
  vector<lower=min_x, upper=max_x / 3>[N_exp] mu_cen;

  // Global asymptotic high-dye log ratio from the centering/control peak to the
  // 2n tumor peak.
  real<lower=log(1), upper=log(20)> log_ratio_tumor_cen_mu;

  // Global log increments for the 4n/2n and 8n/4n peak ratios. The 4n/2n
  // increment is bounded so the derived ratio stays within [1.95, 2.05].
  real<lower=log_R_4n_lower, upper=log_R_4n_upper> log_R_4n_mu;
  real<lower=log(0.1), upper=log(2)> log_R_8n_mu;

  // Reduced mechanistic dye-depletion parameters on the log scale.
  real<lower=-8, upper=8> log_D_tilde;
  real<lower=-8, upper=8> log_A_C_tilde;
  real<lower=-12, upper=4> log_A_T_tilde;
  real<lower=-4, upper=4> log_rho;
  real<lower=1e-6> sigma_group_extra_load;
  vector[N_group] group_extra_load;

  // Shared coefficient of variation for all Gaussian peak widths. A peak with
  // mean mu has standard deviation mu * cv.
  real<lower=0.01, upper=0.5> cv;

  // Sample-specific mixture weights for three latent populations: centering
  // chicken/control, diploid tumor, and WGD tumor.
  array[N_sample] simplex[3] pop_theta;

  // Global allocation of the centering/control population into singlet and
  // doublet noise peaks.
  simplex[2] cen_phase;

  // Global cell-cycle baselines in additive log-ratio coordinates. Entries are
  // interpreted as offsets relative to G0/G1, then mapped through softmax to
  // yield group-specific simplexes over G0/G1, S phase, and G2/M.
  vector[2] diploid_phase_global_logit;
  vector[2] wgd_phase_global_logit;

  // Group-level random effects on the tumor cell-cycle composition. Tight
  // priors on these scales keep each group close to the global ratios.
  vector<lower=1e-6>[2] sigma_group_diploid_phase;
  vector<lower=1e-6>[2] sigma_group_wgd_phase;
  matrix[2, N_group] z_group_diploid_phase;
  matrix[2, N_group] z_group_wgd_phase;
}

transformed parameters {
  // Sample-specific derived burden/load proxies, reduced mechanistic free-dye
  // states, peak ratios, locations, and standard deviations.
  vector[N_sample] load_proxy;
  vector[N_sample] observed_load_proxy;
  vector[N_sample] W_sample;
  vector[N_sample] u_sample;
  vector[N_sample] sample_extra_load;
  vector[N_sample] dye_scale_sample;
  vector[N_sample] ratio_tumor_cen_sample;
  vector[N_sample] R_4n_sample;
  vector[N_sample] R_8n_sample;
  vector[N_sample] mu_cen_sample;
  vector[N_sample] mu_2n_sample;
  vector[N_sample] mu_4n_sample;
  vector[N_sample] mu_8n_sample;
  vector[N_sample] sig_cen_sample;
  vector[N_sample] sig_2n_sample;
  vector[N_sample] sig_4n_sample;
  vector[N_sample] sig_8n_sample;
  vector[N_sample] tumor_s_phase_fraction;
  array[N_group] simplex[3] diploid_phase_group;
  array[N_group] simplex[3] wgd_phase_group;
  real D_tilde = exp(log_D_tilde);
  real A_C_tilde = exp(log_A_C_tilde);
  real A_T_tilde = exp(log_A_T_tilde);
  real rho = exp(log_rho);
  array[N_sample] vector[7] theta;

  for (g in 1:N_group) {
    vector[3] diploid_eta;
    vector[3] wgd_eta;

    diploid_eta[1] = 0;
    diploid_eta[2] = diploid_phase_global_logit[1] +
      sigma_group_diploid_phase[1] * z_group_diploid_phase[1, g];
    diploid_eta[3] = diploid_phase_global_logit[2] +
      sigma_group_diploid_phase[2] * z_group_diploid_phase[2, g];
    diploid_phase_group[g] = softmax(diploid_eta);

    wgd_eta[1] = 0;
    wgd_eta[2] = wgd_phase_global_logit[1] +
      sigma_group_wgd_phase[1] * z_group_wgd_phase[1, g];
    wgd_eta[3] = wgd_phase_global_logit[2] +
      sigma_group_wgd_phase[2] * z_group_wgd_phase[2, g];
    wgd_phase_group[g] = softmax(wgd_eta);
  }

  for (s in 1:N_sample) {
    int e = sampleExpID[s];
    int g = sampleGroupID[s];
    real tumor_dna_load;
    real cen_fraction_safe;
    real distortion_factor;

    // Deterministic observed-component masses from the latent population model:
    // 1 CEN singlet, 2 diploid G0/G1 at 2n, 3 combined 4n peak, 4 WGD G2/M at
    // 8n, 5 diploid S phase, 6 WGD S phase, 7 CEN doublet.
    theta[s][1] = pop_theta[s][1] * cen_phase[1];
    theta[s][2] = pop_theta[s][2] * diploid_phase_group[g][1];
    theta[s][3] = pop_theta[s][2] * diploid_phase_group[g][3] + pop_theta[s][3] * wgd_phase_group[g][1];
    theta[s][4] = pop_theta[s][3] * wgd_phase_group[g][3];
    theta[s][5] = pop_theta[s][2] * diploid_phase_group[g][2];
    theta[s][6] = pop_theta[s][3] * wgd_phase_group[g][2];
    theta[s][7] = pop_theta[s][1] * cen_phase[2];
    tumor_s_phase_fraction[s] = (theta[s][5] + theta[s][6]) / (pop_theta[s][2] + pop_theta[s][3]);

    // Load proxy now uses all tumor-associated observed components, weighted by
    // their effective DNA content, and divides by an external observed CEN
    // fraction derived from the sample histogram below a fixed threshold.
    tumor_dna_load =
      2 * theta[s][2] +
      4 * theta[s][3] +
      8 * theta[s][4] +
      3 * theta[s][5] +
      6 * theta[s][6];
    cen_fraction_safe = fmax(observed_cen_fraction[s], 1e-9);
    observed_load_proxy[s] = log(tumor_dna_load / cen_fraction_safe);
    sample_extra_load[s] = group_extra_load[sampleGroupID[s]];
    load_proxy[s] = observed_load_proxy[s] + sample_extra_load[s];
    W_sample[s] = exp(load_proxy[s]);

    // Solve the reduced mechanistic mass-balance equation for the
    // non-dimensional free-dye state u_s, then map it into the CEN scale and
    // the shared tumor-vs-CEN distortion factor.
    u_sample[s] = solve_reduced_mechanistic_u(D_tilde, A_C_tilde, A_T_tilde, rho, W_sample[s]);
    dye_scale_sample[s] = u_sample[s] / (1 + u_sample[s]);
    distortion_factor = (1 + u_sample[s]) / (1 + rho * u_sample[s]);
    ratio_tumor_cen_sample[s] =
      (1 + min_tumor_cen_ratio_gap + exp(log_ratio_tumor_cen_mu)) * distortion_factor;
    R_4n_sample[s] = 1 + min_ploidy_ratio_gap + exp(log_R_4n_mu);
    R_8n_sample[s] = 1 + min_ploidy_ratio_gap + exp(log_R_8n_mu);

    // Build the ordered ploidy peak locations. The sample-specific CEN peak is
    // the asymptotic experiment baseline scaled by the free-dye state, while
    // all tumor peaks share the same sample-specific distortion relative to
    // CEN and therefore preserve constant tumor-to-tumor spacing.
    mu_cen_sample[s] = mu_cen[e] * dye_scale_sample[s];
    sig_cen_sample[s] = mu_cen_sample[s] * cv;
    mu_2n_sample[s] = mu_cen_sample[s] * ratio_tumor_cen_sample[s];
    mu_4n_sample[s] = mu_2n_sample[s] * R_4n_sample[s];
    mu_8n_sample[s] = mu_4n_sample[s] * R_8n_sample[s];

    // Peak widths scale with peak location through the common CV.
    sig_2n_sample[s] = mu_2n_sample[s] * cv;
    sig_4n_sample[s] = mu_4n_sample[s] * cv;
    sig_8n_sample[s] = mu_8n_sample[s] * cv;
  }
}

model {
  // Priors for the global asymptotic tumor/control ratio, the reduced
  // mechanistic mass-balance parameters, and the population mean
  // ploidy-to-ploidy ratios. log(9) says the 2n tumor peak is expected to sit
  // about 9x above the centering peak in the high-dye limit. log_R priors near
  // log(0.5) imply
  // 1 + min_ploidy_ratio_gap + exp(log(0.5)) = 2, so 4n is expected near 2x
  // 2n and 8n near 2x 4n.
  log_ratio_tumor_cen_mu ~ normal(log(8), 0.35);
  log_R_4n_mu ~ normal(log(0.5), 0.05);
  log_R_8n_mu ~ normal(log(0.5), 0.05);
  log_D_tilde ~ normal(log(5), 1);
  log_A_C_tilde ~ normal(log(1), 1);
  log_A_T_tilde ~ normal(log(0.05), 1);
  log_rho ~ normal(0, 0.35);
  sigma_group_extra_load ~ normal(0, 0.35);
  group_extra_load ~ normal(0, sigma_group_extra_load);

  for (s in 1:N_sample) {
    pop_theta[s] ~ dirichlet([2, 4, 2]');

    // S phase should remain a small share of the latent tumor fraction. This
    // soft upper penalty is zero at 20% and increasingly discourages larger
    // values without imposing a hard boundary on the simplex parameters.
    target += -0.5 * square(fmax(tumor_s_phase_fraction[s] - 0.2, 0) / 0.02);
  }
  cen_phase ~ dirichlet([8, 0.5]');
  diploid_phase_global_logit[1] ~ normal(log(0.75 / 4), 0.35);
  diploid_phase_global_logit[2] ~ normal(log(1.0 / 4), 0.35);
  wgd_phase_global_logit[1] ~ normal(log(0.75 / 4), 0.35);
  wgd_phase_global_logit[2] ~ normal(log(1.0 / 4), 0.35);
  sigma_group_diploid_phase ~ normal(0, 0.05);
  sigma_group_wgd_phase ~ normal(0, 0.05);
  to_vector(z_group_diploid_phase) ~ normal(0, 1);
  to_vector(z_group_wgd_phase) ~ normal(0, 1);
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
      comp[1] = log(theta[s][1]) + normal_lpdf(x | mu_cen_sample[s], sig_cen_sample[s]);
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
        2 * mu_cen_sample[s],
        sqrt(2) * sig_cen_sample[s]
      );

      // Marginalize over the unknown component label for this bin center.
      target += count[i] * log_sum_exp(comp);
    }
  }
}
