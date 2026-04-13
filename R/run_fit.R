library(cmdstanr)
library(ggplot2)

source(file.path("R", "fit_specs.R"))

make_simplex_vector <- function(weights, expected_length, label) {
  if (length(weights) != expected_length) {
    stop(sprintf("%s initialization requires exactly %d weights.", label, expected_length))
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop(sprintf("%s initialization weights must be finite and non-negative.", label))
  }
  as.numeric(weights / sum(weights))
}

simplex_to_baseline_logits <- function(simplex_weights, label) {
  simplex_weights <- make_simplex_vector(simplex_weights, 3L, label)
  log(simplex_weights[2:3] / simplex_weights[1])
}

make_pop_theta_init <- function(n_sample, weights = c(1, 2, 1)) {
  simplex_weights <- make_simplex_vector(weights, 3L, "pop_theta")
  matrix(rep(simplex_weights, times = n_sample), nrow = n_sample, byrow = TRUE)
}

make_fixed_init <- function(stan_data, spec) {
  cen_phase <- make_simplex_vector(c(8, 0.5), 2L, "cen_phase")
  diploid_phase <- make_simplex_vector(c(4, 0.75, 1), 3L, "diploid_phase")
  wgd_phase <- make_simplex_vector(c(4, 0.75, 1), 3L, "wgd_phase")

  if (identical(spec$model_variant, "first_pass")) {
    return(list(
      mu_cen = rep(5000, stan_data$N_exp),
      log_ratio_tumor_cen_mu = log(8),
      log_R_4n_mu = log(0.5),
      log_R_8n_mu = log(0.5),
      cv = 0.3,
      pop_theta = make_pop_theta_init(stan_data$N_sample),
      cen_phase = cen_phase,
      diploid_phase = diploid_phase,
      wgd_phase = wgd_phase
    ))
  }

  init <- list(
    mu_cen = rep(5000, stan_data$N_exp),
    log_ratio_tumor_cen_mu = log(8),
    log_R_4n_mu = log(0.5),
    log_R_8n_mu = log(0.5),
    beta = 1e-1,
    cv = 0.3,
    pop_theta = make_pop_theta_init(stan_data$N_sample),
    cen_phase = cen_phase
  )

  if (identical(spec$model_variant, "original")) {
    init$diploid_phase <- diploid_phase
    init$wgd_phase <- wgd_phase
    return(init)
  }

  init$log_D_tilde <- log(5)
  init$log_A_C_tilde <- log(1)
  init$log_A_T_tilde <- log(0.05)
  init$log_rho <- 0
  init$sigma_group_extra_load <- 0.05
  init$group_extra_load <- rep(0, stan_data$N_group)
  init$diploid_phase_global_logit <- simplex_to_baseline_logits(diploid_phase, "diploid_phase")
  init$wgd_phase_global_logit <- simplex_to_baseline_logits(wgd_phase, "wgd_phase")
  init$sigma_group_diploid_phase <- c(0.02, 0.02)
  init$sigma_group_wgd_phase <- c(0.02, 0.02)
  init$z_group_diploid_phase <- matrix(0, nrow = 2, ncol = stan_data$N_group)
  init$z_group_wgd_phase <- matrix(0, nrow = 2, ncol = stan_data$N_group)
  init
}

extract_indexed_vector <- function(row_df, prefix, n) {
  vapply(seq_len(n), function(i) row_df[[sprintf("%s.%d", prefix, i)]], numeric(1))
}

extract_theta_matrix <- function(row_df, n_sample) {
  theta <- matrix(NA_real_, nrow = n_sample, ncol = 7)
  for (s in seq_len(n_sample)) {
    for (k in seq_len(7)) {
      theta[s, k] <- row_df[[sprintf("theta.%d.%d", s, k)]]
    }
  }
  theta
}

build_group_phase_matrix <- function(global_logit, sigma_group, z_group) {
  n_group <- ncol(z_group)
  group_phase <- matrix(NA_real_, nrow = n_group, ncol = 3)
  for (g in seq_len(n_group)) {
    eta <- c(0, global_logit + sigma_group * z_group[, g])
    eta <- eta - max(eta)
    group_phase[g, ] <- exp(eta) / sum(exp(eta))
  }
  group_phase
}

build_theta_from_population <- function(pop_theta, cen_phase, diploid_phase, wgd_phase) {
  if (is.null(dim(diploid_phase))) {
    diploid_phase <- matrix(rep(diploid_phase, each = nrow(pop_theta)), nrow = nrow(pop_theta))
  }
  if (is.null(dim(wgd_phase))) {
    wgd_phase <- matrix(rep(wgd_phase, each = nrow(pop_theta)), nrow = nrow(pop_theta))
  }

  theta <- matrix(NA_real_, nrow = nrow(pop_theta), ncol = 7)
  theta[, 1] <- pop_theta[, 1] * cen_phase[1]
  theta[, 2] <- pop_theta[, 2] * diploid_phase[, 1]
  theta[, 3] <- pop_theta[, 2] * diploid_phase[, 3] + pop_theta[, 3] * wgd_phase[, 1]
  theta[, 4] <- pop_theta[, 3] * wgd_phase[, 3]
  theta[, 5] <- pop_theta[, 2] * diploid_phase[, 2]
  theta[, 6] <- pop_theta[, 3] * wgd_phase[, 2]
  theta[, 7] <- pop_theta[, 1] * cen_phase[2]
  theta
}

compute_load_proxy <- function(pop_theta) {
  log((2 * pop_theta[, 2] + 4 * pop_theta[, 3]) / pop_theta[, 1])
}

compute_observed_cen_fraction <- function(stan_data) {
  sample_totals <- rowsum(stan_data$count, group = stan_data$sampleID, reorder = FALSE)
  cen_totals <- rowsum(stan_data$count * (stan_data$areaDNA_bin_center < stan_data$cen_threshold), group = stan_data$sampleID, reorder = FALSE)
  as.numeric(cen_totals[, 1] / pmax(sample_totals[, 1], 1))
}

compute_load_proxy_observed_cen_scaled <- function(theta, observed_cen_fraction) {
  tumor_dna_load <- 2 * theta[, 2] + 4 * theta[, 3] + 8 * theta[, 4] + 3 * theta[, 5] + 6 * theta[, 6]
  log(tumor_dna_load / pmax(observed_cen_fraction, 1e-9))
}

solve_reduced_mechanistic_u <- function(D_tilde, A_C_tilde, A_T_tilde, rho, W_s, max_iter = 100L, tol = 1e-10) {
  u <- D_tilde / (1 + A_C_tilde + A_T_tilde * W_s * rho)
  for (iter in seq_len(max_iter)) {
    f <- u + A_C_tilde * u / (1 + u) + A_T_tilde * W_s * (rho * u) / (1 + rho * u) - D_tilde
    fp <- 1 + A_C_tilde / (1 + u)^2 + A_T_tilde * W_s * rho / (1 + rho * u)^2
    u_next <- u - f / fp
    if (!is.finite(u_next) || u_next <= 0) u_next <- 0.5 * u
    if (abs(u_next - u) < tol * max(1, u)) return(u_next)
    u <- u_next
  }
  u
}

s_phase_density <- function(x, lower_peak, upper_peak, n_subpeaks) {
  width <- upper_peak - lower_peak
  if (!is.finite(width) || width <= 0) return(0)
  dnorm(x, mean = 0.5 * (lower_peak + upper_peak), sd = 0.5 * width)
}

build_overlay_df <- function(stan_data, lut, mu_cen, ratio_tumor_cen_sample, R_4n_sample, R_8n_sample, cv, theta, dye_scale_sample = NULL) {
  if (is.null(dye_scale_sample)) dye_scale_sample <- rep(1, stan_data$N_sample)
  x <- stan_data$areaDNA_bin_center
  bin_width <- diff(sort(unique(x)))[1]

  plot_df <- data.frame(areaDNA_bin_center = stan_data$areaDNA_bin_center, count = stan_data$count, expID = stan_data$expID, sampleID = stan_data$sampleID)
  sample_totals <- aggregate(count ~ sampleID, data = plot_df, FUN = sum)
  names(sample_totals)[2] <- "sample_total"
  plot_df <- merge(plot_df, sample_totals, by = "sampleID", sort = FALSE)

  plot_df$mixture_count <- mapply(function(x_i, sample_id, exp_id, sample_total) {
    mu_cen_sample <- mu_cen[exp_id] * dye_scale_sample[sample_id]
    sig_cen <- mu_cen_sample * cv
    mu_2n <- mu_cen_sample * ratio_tumor_cen_sample[sample_id]
    mu_4n <- mu_2n * R_4n_sample[sample_id]
    mu_8n <- mu_4n * R_8n_sample[sample_id]
    sig_2n <- mu_2n * cv
    sig_4n <- mu_4n * cv
    sig_8n <- mu_8n * cv
    density_val <- theta[sample_id, 1] * dnorm(x_i, mean = mu_cen_sample, sd = sig_cen) +
      theta[sample_id, 2] * dnorm(x_i, mean = mu_2n, sd = sig_2n) +
      theta[sample_id, 3] * dnorm(x_i, mean = mu_4n, sd = sig_4n) +
      theta[sample_id, 4] * dnorm(x_i, mean = mu_8n, sd = sig_8n) +
      theta[sample_id, 5] * s_phase_density(x_i, mu_2n, mu_4n, stan_data$N_s_phase) +
      theta[sample_id, 6] * s_phase_density(x_i, mu_4n, mu_8n, stan_data$N_s_phase) +
      theta[sample_id, 7] * dnorm(x_i, mean = 2 * mu_cen_sample, sd = sqrt(2) * sig_cen)
    sample_total * bin_width * density_val
  }, plot_df$areaDNA_bin_center, plot_df$sampleID, plot_df$expID, plot_df$sample_total)

  plot_df <- merge(plot_df, lut$sample[, c("sample_id", "sample_name", "experiment_name")], by.x = "sampleID", by.y = "sample_id", all.x = TRUE, sort = FALSE)
  plot_df$sample_label <- sprintf("%s | %s", plot_df$experiment_name, plot_df$sample_name)
  plot_df
}

save_overlay_plot <- function(plot_df, plot_path, title_text) {
  bin_width <- diff(sort(unique(plot_df$areaDNA_bin_center)))[1]
  p <- ggplot(plot_df, aes(x = areaDNA_bin_center)) +
    geom_col(aes(y = count), width = bin_width, fill = "grey82", color = "grey55", linewidth = 0.15) +
    geom_line(aes(y = mixture_count), color = "#d94801", linewidth = 0.5) +
    facet_wrap(~ sample_label, scales = "free_y") +
    labs(title = title_text, x = "Raw DNA area bin center", y = "Bin count") +
    theme_bw(base_size = 9) +
    theme(strip.text = element_text(size = 7), panel.grid.minor = element_blank())
  ggsave(plot_path, p, width = 16, height = 12, dpi = 200)
}

run_fit <- function(spec_name = "first_pass", project_root = ".") {
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  spec <- get_fit_spec(spec_name, project_root = project_root)
  prepared <- prepare_fit_inputs(spec, project_root = project_root)
  stan_data <- prepared$stan_data
  lut <- prepared$lut
  dir.create(spec$fit_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(spec$plot_dir, recursive = TRUE, showWarnings = FALSE)

  fixed_init <- make_fixed_init(stan_data, spec)
  min_tumor_cen_ratio_gap <- 1e-3
  min_ploidy_ratio_gap <- 0.5

  if (identical(spec$model_variant, "first_pass")) {
    theta_init <- build_theta_from_population(fixed_init$pop_theta, fixed_init$cen_phase, fixed_init$diploid_phase, fixed_init$wgd_phase)
    ratio_init <- rep(1 + min_tumor_cen_ratio_gap + exp(fixed_init$log_ratio_tumor_cen_mu), stan_data$N_sample)
    dye_scale_init <- rep(1, stan_data$N_sample)
  } else if (identical(spec$model_variant, "original")) {
    theta_init <- build_theta_from_population(fixed_init$pop_theta, fixed_init$cen_phase, fixed_init$diploid_phase, fixed_init$wgd_phase)
    ratio_init <- 1 + min_tumor_cen_ratio_gap + exp(fixed_init$log_ratio_tumor_cen_mu - fixed_init$beta * compute_load_proxy(fixed_init$pop_theta))
    dye_scale_init <- rep(1, stan_data$N_sample)
  } else {
    diploid_phase_group_init <- build_group_phase_matrix(fixed_init$diploid_phase_global_logit, fixed_init$sigma_group_diploid_phase, fixed_init$z_group_diploid_phase)
    wgd_phase_group_init <- build_group_phase_matrix(fixed_init$wgd_phase_global_logit, fixed_init$sigma_group_wgd_phase, fixed_init$z_group_wgd_phase)
    theta_init <- build_theta_from_population(fixed_init$pop_theta, fixed_init$cen_phase, diploid_phase_group_init[stan_data$sampleGroupID, , drop = FALSE], wgd_phase_group_init[stan_data$sampleGroupID, , drop = FALSE])
    u_sample_init <- vapply(
      exp(compute_load_proxy_observed_cen_scaled(theta_init, compute_observed_cen_fraction(stan_data))),
      function(W_s) solve_reduced_mechanistic_u(exp(fixed_init$log_D_tilde), exp(fixed_init$log_A_C_tilde), exp(fixed_init$log_A_T_tilde), exp(fixed_init$log_rho), W_s),
      numeric(1)
    )
    dye_scale_init <- u_sample_init / (1 + u_sample_init)
    ratio_init <- rep(1 + min_tumor_cen_ratio_gap + exp(fixed_init$log_ratio_tumor_cen_mu), stan_data$N_sample) * (1 + u_sample_init) / (1 + exp(fixed_init$log_rho) * u_sample_init)
  }

  R_4n_init <- rep(1 + min_ploidy_ratio_gap + exp(fixed_init$log_R_4n_mu), stan_data$N_sample)
  R_8n_init <- rep(1 + min_ploidy_ratio_gap + exp(fixed_init$log_R_8n_mu), stan_data$N_sample)
  init_plot_df <- build_overlay_df(stan_data, lut, fixed_init$mu_cen, ratio_init, R_4n_init, R_8n_init, fixed_init$cv, theta_init, dye_scale_init)
  init_plot_path <- file.path(spec$plot_dir, sprintf("%s_overlay.png", spec$plot_label_prefix))
  opt_plot_path <- file.path(spec$plot_dir, sprintf("%s_optimized_overlay.png", spec$plot_label_prefix))
  save_overlay_plot(init_plot_df, init_plot_path, sprintf("Fixed initialization overlay (%s)", spec$experiment_filter))

  mod <- cmdstan_model(spec$stan_file, compile = TRUE, exe_file = spec$exe_file, force_recompile = FALSE)
  fit <- mod$optimize(
    data = stan_data,
    algorithm = "lbfgs",
    iter = 5000,
    seed = 123,
    init = list(fixed_init),
    output_dir = spec$fit_root,
    output_basename = spec$output_basename
  )

  opt_row <- read.csv(fit$output_files(), comment.char = "#", check.names = FALSE)[1, , drop = FALSE]
  dye_scale_opt <- if (identical(spec$model_variant, "observed_cen_scaled")) extract_indexed_vector(opt_row, "dye_scale_sample", stan_data$N_sample) else rep(1, stan_data$N_sample)
  opt_plot_df <- build_overlay_df(
    stan_data,
    lut,
    extract_indexed_vector(opt_row, "mu_cen", stan_data$N_exp),
    extract_indexed_vector(opt_row, "ratio_tumor_cen_sample", stan_data$N_sample),
    extract_indexed_vector(opt_row, "R_4n_sample", stan_data$N_sample),
    extract_indexed_vector(opt_row, "R_8n_sample", stan_data$N_sample),
    opt_row$cv[1],
    extract_theta_matrix(opt_row, stan_data$N_sample),
    dye_scale_opt
  )
  save_overlay_plot(opt_plot_df, opt_plot_path, sprintf("Optimization result from fixed initialization (%s)", spec$experiment_filter))

  manifest <- list(
    generated_at = Sys.time(),
    spec_name = spec$spec_name,
    model_variant = spec$model_variant,
    experiment_filter = spec$experiment_filter,
    exclude_sample_patterns = spec$exclude_sample_patterns,
    excluded_samples = prepared$excluded_samples,
    cen_threshold = if (!is.null(spec$cen_threshold)) stan_data$cen_threshold else NULL,
    experiment_lut = lut$experiment,
    sample_lut = lut$sample,
    stan_file = normalizePath(spec$stan_file, winslash = "/", mustWork = TRUE),
    stan_data_path = normalizePath(file.path(project_root, "processed_data", "stan_data.Rds"), winslash = "/", mustWork = TRUE),
    output_file = normalizePath(fit$output_files(), winslash = "/", mustWork = TRUE),
    init = fixed_init,
    init_plot_path = normalizePath(init_plot_path, winslash = "/", mustWork = TRUE),
    opt_plot_path = normalizePath(opt_plot_path, winslash = "/", mustWork = TRUE)
  )
  saveRDS(manifest, spec$manifest_path)
  cat(sprintf("Saved fit manifest to %s\n", spec$manifest_path))
}

if (sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  run_fit(
    spec_name = if (length(args) >= 1) args[[1]] else "first_pass",
    project_root = if (length(args) >= 2) args[[2]] else "."
  )
}
