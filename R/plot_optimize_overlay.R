library(ggplot2)

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
stan_data_path <- file.path(project_root, "processed_data", "stan_data.Rds")
fit_root <- file.path(project_root, "processed_data", "stan_optimize_first_pass")
plot_dir <- file.path(project_root, "figure", "stan_optimize_first_pass")

manifest_files <- Sys.glob(file.path(fit_root, "optimize_manifest_seed_*.rds"))
if (length(manifest_files) == 0L) {
  stop("No optimize manifests found under processed_data/stan_optimize_first_pass.")
}

manifest_lp <- vapply(manifest_files, function(path) {
  manifest <- readRDS(path)
  draw_row <- read.csv(manifest$output_file, comment.char = "#", check.names = FALSE)[1, , drop = FALSE]
  draw_row$lp__[1]
}, numeric(1))
manifest_path <- manifest_files[[which.max(manifest_lp)]]
plot_path <- file.path(plot_dir, sprintf("optimize_overlay_%s.png", tools::file_path_sans_ext(basename(manifest_path))))

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

stan_bundle <- readRDS(stan_data_path)
stan_data <- stan_bundle$stan_data
if (is.null(stan_data$N_s_phase)) {
  stan_data$N_s_phase <- 5L
}
lut <- stan_bundle$lut
manifest <- readRDS(manifest_path)
draw_row <- read.csv(manifest$output_file, comment.char = "#", check.names = FALSE)[1, , drop = FALSE]

extract_indexed_vector <- function(draw_row, prefix, n) {
  vapply(seq_len(n), function(i) draw_row[[sprintf("%s.%d", prefix, i)]], numeric(1))
}

extract_theta_matrix <- function(draw_row, n_sample) {
  theta <- matrix(NA_real_, nrow = n_sample, ncol = 7)
  for (s in seq_len(n_sample)) {
    for (k in seq_len(7)) {
      theta[s, k] <- draw_row[[sprintf("theta.%d.%d", s, k)]]
    }
  }
  theta
}

extract_sample_vector <- function(draw_row, prefix, n_sample) {
  vapply(seq_len(n_sample), function(i) draw_row[[sprintf("%s.%d", prefix, i)]], numeric(1))
}

s_phase_density <- function(x, lower_peak, upper_peak, n_subpeaks) {
  width <- upper_peak - lower_peak
  if (!is.finite(width) || width <= 0) {
    return(0)
  }
  dnorm(x, mean = 0.5 * (lower_peak + upper_peak), sd = 0.5 * width)
}

mu_cen <- extract_indexed_vector(draw_row, "mu_cen", stan_data$N_exp)
ratio_tumor_cen_sample <- extract_sample_vector(draw_row, "ratio_tumor_cen_sample", stan_data$N_sample)
R_4n_sample <- extract_sample_vector(draw_row, "R_4n_sample", stan_data$N_sample)
R_8n_sample <- extract_sample_vector(draw_row, "R_8n_sample", stan_data$N_sample)
cv <- draw_row[["cv"]]
theta <- extract_theta_matrix(draw_row, stan_data$N_sample)

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
      theta[sample_id, 5] * s_phase_density(x_i, mu_2n, mu_4n, stan_data$N_s_phase) +
      theta[sample_id, 6] * s_phase_density(x_i, mu_4n, mu_8n, stan_data$N_s_phase) +
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

p <- ggplot(plot_df, aes(x = areaDNA_bin_center)) +
  geom_col(aes(y = count), width = bin_width, fill = "grey82", color = "grey55", linewidth = 0.15) +
  geom_line(aes(y = mixture_count), color = "#d94801", linewidth = 0.5) +
  facet_wrap(~ sample_label, scales = "free_y") +
  labs(
    title = sprintf("Optimization overlay (%s)", basename(manifest_path)),
    x = "Raw DNA area bin center",
    y = "Bin count"
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.text = element_text(size = 7),
    panel.grid.minor = element_blank()
  )

ggsave(plot_path, p, width = 16, height = 12, dpi = 200)

cat(sprintf("Saved optimize overlay plot to %s\n", plot_path))
