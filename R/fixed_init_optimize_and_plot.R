library(cmdstanr)
library(ggplot2)

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
stan_file <- file.path(project_root, "stan", "ploidy_histogram_first_pass.stan")
stan_data_path <- file.path(project_root, "processed_data", "stan_data.Rds")
experiment_filter <- "anoxia"
fit_root <- file.path(project_root, "processed_data", sprintf("stan_optimize_fixed_init_%s", experiment_filter))
plot_dir <- file.path(project_root, "figure", sprintf("stan_optimize_fixed_init_%s", experiment_filter))
exe_file <- file.path(fit_root, "ploidy_histogram_first_pass_fixed_init.exe")
init_plot_path <- file.path(plot_dir, "fixed_init_overlay.png")
opt_plot_path <- file.path(plot_dir, "fixed_init_optimized_overlay.png")
output_basename <- "fixed_init_optimize"
manifest_path <- file.path(fit_root, "fixed_init_optimize_manifest.rds")

dir.create(fit_root, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

stan_bundle <- readRDS(stan_data_path)
stan_data <- stan_bundle$stan_data
if (is.null(stan_data$N_s_phase)) {
  stan_data$N_s_phase <- 3L
}
lut <- stan_bundle$lut
min_tumor_cen_ratio_gap <- 1e-3
min_ploidy_ratio_gap <- 0.5

subset_stan_bundle_by_experiment <- function(stan_data, lut, experiment_pattern) {
  keep_experiments <- grepl(experiment_pattern, lut$experiment$experiment_name, ignore.case = TRUE)
  if (!any(keep_experiments)) {
    stop(sprintf("No experiments matched pattern '%s'.", experiment_pattern))
  }

  old_exp_ids <- lut$experiment$experiment_id[keep_experiments]
  keep_samples <- lut$sample$experiment_id %in% old_exp_ids
  old_sample_ids <- lut$sample$sample_id[keep_samples]
  keep_bins <- stan_data$sampleID %in% old_sample_ids

  exp_map <- setNames(seq_along(old_exp_ids), old_exp_ids)
  sample_map <- setNames(seq_along(old_sample_ids), old_sample_ids)

  experiment_lut <- lut$experiment[keep_experiments, , drop = FALSE]
  experiment_lut$experiment_id <- seq_len(nrow(experiment_lut))

  sample_lut <- lut$sample[keep_samples, , drop = FALSE]
  sample_lut$sample_id <- seq_len(nrow(sample_lut))
  sample_lut$experiment_id <- unname(exp_map[as.character(sample_lut$experiment_id)])

  stan_data$areaDNA_bin_center <- stan_data$areaDNA_bin_center[keep_bins]
  stan_data$count <- stan_data$count[keep_bins]
  stan_data$expID <- unname(exp_map[as.character(stan_data$expID[keep_bins])])
  stan_data$sampleID <- unname(sample_map[as.character(stan_data$sampleID[keep_bins])])
  stan_data$sampleExpID <- unname(exp_map[as.character(lut$sample$experiment_id[keep_samples])])
  stan_data$N <- length(stan_data$count)
  stan_data$N_exp <- nrow(experiment_lut)
  stan_data$N_sample <- nrow(sample_lut)

  list(
    stan_data = stan_data,
    lut = list(
      experiment = experiment_lut,
      sample = sample_lut
    )
  )
}

filtered_bundle <- subset_stan_bundle_by_experiment(stan_data, lut, experiment_filter)
stan_data <- filtered_bundle$stan_data
lut <- filtered_bundle$lut

default_pop_weights <- c(
  centering = 1,
  diploid = 2,
  wgd = 1
)

default_cen_phase <- c(
  singlet = 8,
  doublet = 0.5
)

default_diploid_phase <- c(
  g0g1 = 4,
  s_phase = 0.75,
  g2m = 1
)

default_wgd_phase <- c(
  g0g1 = 4,
  s_phase = 0.75,
  g2m = 1
)

make_simplex_vector <- function(weights, expected_length, label) {
  if (length(weights) != expected_length) {
    stop(sprintf("%s initialization requires exactly %d weights.", label, expected_length))
  }

  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop(sprintf("%s initialization weights must be finite and non-negative.", label))
  }

  weight_sum <- sum(weights)
  if (weight_sum <= 0) {
    stop(sprintf("At least one %s initialization weight must be positive.", label))
  }

  as.numeric(weights / weight_sum)
}

make_pop_theta_init <- function(n_sample, weights = default_pop_weights) {
  simplex_weights <- make_simplex_vector(weights, 3L, "pop_theta")
  pop_theta <- matrix(rep(simplex_weights, times = n_sample), nrow = n_sample, byrow = TRUE)

  if (any(abs(rowSums(pop_theta) - 1) > 1e-12)) {
    stop("Internal error: pop_theta initialization rows do not sum to 1.")
  }

  pop_theta
}

make_fixed_init <- function(stan_data) {
  list(
    mu_cen = rep(5000, stan_data$N_exp),
    log_ratio_tumor_cen_mu = log(8),
    log_R_4n_mu = log(0.5),
    log_R_8n_mu = log(0.5),
    cv = 0.3,
    pop_theta = make_pop_theta_init(stan_data$N_sample),
    cen_phase = make_simplex_vector(default_cen_phase, 2L, "cen_phase"),
    diploid_phase = make_simplex_vector(default_diploid_phase, 3L, "diploid_phase"),
    wgd_phase = make_simplex_vector(default_wgd_phase, 3L, "wgd_phase")
  )
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

build_theta_from_population <- function(pop_theta, cen_phase, diploid_phase, wgd_phase) {
  theta <- matrix(NA_real_, nrow = nrow(pop_theta), ncol = 7)
  theta[, 1] <- pop_theta[, 1] * cen_phase[[1]]
  theta[, 2] <- pop_theta[, 2] * diploid_phase[[1]]
  theta[, 3] <- pop_theta[, 2] * diploid_phase[[3]] + pop_theta[, 3] * wgd_phase[[1]]
  theta[, 4] <- pop_theta[, 3] * wgd_phase[[3]]
  theta[, 5] <- pop_theta[, 2] * diploid_phase[[2]]
  theta[, 6] <- pop_theta[, 3] * wgd_phase[[2]]
  theta[, 7] <- pop_theta[, 1] * cen_phase[[2]]
  theta
}

s_phase_density <- function(x, lower_peak, upper_peak, n_subpeaks) {
  width <- upper_peak - lower_peak
  if (!is.finite(width) || width <= 0) {
    return(0)
  }
  dnorm(x, mean = 0.5 * (lower_peak + upper_peak), sd = 0.5 * width)
}

build_overlay_df <- function(stan_data, lut, mu_cen, ratio_tumor_cen_sample, R_4n_sample, R_8n_sample, cv, theta) {
  x <- stan_data$areaDNA_bin_center
  bin_width <- diff(sort(unique(x)))[1]

  plot_df <- data.frame(
    areaDNA_bin_center = stan_data$areaDNA_bin_center,
    count = stan_data$count,
    expID = stan_data$expID,
    sampleID = stan_data$sampleID,
    stringsAsFactors = FALSE
  )

  sample_totals <- aggregate(count ~ sampleID, data = plot_df, FUN = sum)
  names(sample_totals)[2] <- "sample_total"
  plot_df <- merge(plot_df, sample_totals, by = "sampleID", sort = FALSE)

  plot_df$mixture_count <- mapply(
    FUN = function(x_i, sample_id, exp_id, sample_total) {
      sig_cen <- mu_cen[exp_id] * cv
      mu_2n <- mu_cen[exp_id] * ratio_tumor_cen_sample[sample_id]
      mu_4n <- mu_2n * R_4n_sample[sample_id]
      mu_8n <- mu_4n * R_8n_sample[sample_id]
      sig_2n <- mu_2n * cv
      sig_4n <- mu_4n * cv
      sig_8n <- mu_8n * cv

      density_val <- theta[sample_id, 1] * dnorm(x_i, mean = mu_cen[exp_id], sd = sig_cen) +
        theta[sample_id, 2] * dnorm(x_i, mean = mu_2n, sd = sig_2n) +
        theta[sample_id, 3] * dnorm(x_i, mean = mu_4n, sd = sig_4n) +
        theta[sample_id, 4] * dnorm(x_i, mean = mu_8n, sd = sig_8n) +
        theta[sample_id, 5] * s_phase_density(
          x_i,
          mu_2n,
          mu_4n,
          stan_data$N_s_phase
        ) +
        theta[sample_id, 6] * s_phase_density(
          x_i,
          mu_4n,
          mu_8n,
          stan_data$N_s_phase
        ) +
        theta[sample_id, 7] * dnorm(x_i, mean = 2 * mu_cen[exp_id], sd = sqrt(2) * sig_cen)

      sample_total * bin_width * density_val
    },
    x_i = plot_df$areaDNA_bin_center,
    sample_id = plot_df$sampleID,
    exp_id = plot_df$expID,
    sample_total = plot_df$sample_total
  )

  plot_df <- merge(
    plot_df,
    lut$sample[, c("sample_id", "sample_name", "experiment_name")],
    by.x = "sampleID",
    by.y = "sample_id",
    all.x = TRUE,
    sort = FALSE
  )

  plot_df$sample_label <- sprintf("%s | %s", plot_df$experiment_name, plot_df$sample_name)
  plot_df
}

save_overlay_plot <- function(plot_df, plot_path, title_text) {
  bin_width <- diff(sort(unique(plot_df$areaDNA_bin_center)))[1]
  p <- ggplot(plot_df, aes(x = areaDNA_bin_center)) +
    geom_col(aes(y = count), width = bin_width, fill = "grey82", color = "grey55", linewidth = 0.15) +
    geom_line(aes(y = mixture_count), color = "#d94801", linewidth = 0.5) +
    facet_wrap(~ sample_label, scales = "free_y") +
    labs(
      title = title_text,
      x = "Raw DNA area bin center",
      y = "Bin count"
    ) +
    theme_bw(base_size = 9) +
    theme(
      strip.text = element_text(size = 7),
      panel.grid.minor = element_blank()
    )

  ggsave(plot_path, p, width = 16, height = 12, dpi = 200)
}

fixed_init <- make_fixed_init(stan_data)
ratio_tumor_cen_sample_init <- rep(
  1 + min_tumor_cen_ratio_gap + exp(fixed_init$log_ratio_tumor_cen_mu),
  stan_data$N_sample
)
R_4n_sample_init <- rep(1 + min_ploidy_ratio_gap + exp(fixed_init$log_R_4n_mu), stan_data$N_sample)
R_8n_sample_init <- rep(1 + min_ploidy_ratio_gap + exp(fixed_init$log_R_8n_mu), stan_data$N_sample)
theta_init <- build_theta_from_population(
  fixed_init$pop_theta,
  fixed_init$cen_phase,
  fixed_init$diploid_phase,
  fixed_init$wgd_phase
)

init_plot_df <- build_overlay_df(
  stan_data = stan_data,
  lut = lut,
  mu_cen = fixed_init$mu_cen,
  ratio_tumor_cen_sample = ratio_tumor_cen_sample_init,
  R_4n_sample = R_4n_sample_init,
  R_8n_sample = R_8n_sample_init,
  cv = fixed_init$cv,
  theta = theta_init
)
save_overlay_plot(init_plot_df, init_plot_path, sprintf("Fixed initialization overlay (%s)", experiment_filter))

mod <- cmdstan_model(stan_file, compile = TRUE, exe_file = exe_file, force_recompile = TRUE)

fit <- mod$optimize(
  data = stan_data,
  algorithm ="lbfgs",
  seed = 123,
  init = list(fixed_init),
  output_dir = fit_root,
  output_basename = output_basename
)

opt_row <- read.csv(fit$output_files(), comment.char = "#", check.names = FALSE)[1, , drop = FALSE]
opt_plot_df <- build_overlay_df(
  stan_data = stan_data,
  lut = lut,
  mu_cen = extract_indexed_vector(opt_row, "mu_cen", stan_data$N_exp),
  ratio_tumor_cen_sample = extract_indexed_vector(opt_row, "ratio_tumor_cen_sample", stan_data$N_sample),
  R_4n_sample = extract_indexed_vector(opt_row, "R_4n_sample", stan_data$N_sample),
  R_8n_sample = extract_indexed_vector(opt_row, "R_8n_sample", stan_data$N_sample),
  cv = opt_row$cv[1],
  theta = extract_theta_matrix(opt_row, stan_data$N_sample)
)
save_overlay_plot(opt_plot_df, opt_plot_path, sprintf("Optimization result from fixed initialization (%s)", experiment_filter))

manifest <- list(
  generated_at = Sys.time(),
  experiment_filter = experiment_filter,
  experiment_lut = lut$experiment,
  sample_lut = lut$sample,
  stan_file = normalizePath(stan_file, winslash = "/", mustWork = TRUE),
  stan_data_path = normalizePath(stan_data_path, winslash = "/", mustWork = TRUE),
  output_file = normalizePath(fit$output_files(), winslash = "/", mustWork = TRUE),
  init = fixed_init,
  init_plot_path = normalizePath(init_plot_path, winslash = "/", mustWork = TRUE),
  opt_plot_path = normalizePath(opt_plot_path, winslash = "/", mustWork = TRUE)
)

saveRDS(manifest, manifest_path)

cat(sprintf("Saved fixed-init manifest to %s\n", manifest_path))
