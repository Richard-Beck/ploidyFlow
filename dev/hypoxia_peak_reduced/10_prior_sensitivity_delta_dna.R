library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(grid)

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1) args[[1]] else ".", winslash = "/", mustWork = TRUE)
spec_name <- if (length(args) >= 2) args[[2]] else "adjacent_8n_delta_dna"

default_scales <- c(0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3)
scale_grid <- if (length(args) >= 3) {
  as.numeric(strsplit(args[[3]], ",", fixed = TRUE)[[1]])
} else {
  default_scales
}

if (any(!is.finite(scale_grid)) || any(scale_grid <= 0)) {
  stop("All prior scales must be positive finite numbers.")
}

scale_grid <- sort(unique(scale_grid))

out_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output")
nuts_dir <- file.path(out_dir, "extensions", spec_name, "nuts")
draws_path <- file.path(nuts_dir, "draws_matrix.rds")
prepared_path <- file.path(out_dir, "prepared_input.rds")

if (!file.exists(draws_path)) {
  stop("Could not find saved draws matrix: ", draws_path)
}
if (!file.exists(prepared_path)) {
  stop("Could not find prepared input: ", prepared_path)
}

draws_matrix <- readRDS(draws_path)
prepared <- readRDS(prepared_path)
input_tbl <- prepared$input_tbl
n_draw <- nrow(draws_matrix)
n_sample <- nrow(input_tbl)

sigma_name <- "sigma_delta_dna"
delta_names <- sprintf("delta_dna[%d]", seq_len(n_sample))
prob_names <- sprintf("prob_4n[%d]", seq_len(n_sample))
lower_names <- sprintf("mu_lower_mixture[%d]", seq_len(n_sample))
upper_names <- sprintf("mu_upper_mixture[%d]", seq_len(n_sample))

required_names <- c(sigma_name, delta_names, prob_names, lower_names, upper_names)
missing_names <- setdiff(required_names, colnames(draws_matrix))
if (length(missing_names)) {
  stop("Saved draws are missing required variables: ", paste(missing_names, collapse = ", "))
}

sigma_draw <- as.numeric(draws_matrix[, sigma_name])
delta_draws <- as.matrix(draws_matrix[, delta_names, drop = FALSE])
prob_draws <- as.matrix(draws_matrix[, prob_names, drop = FALSE])
lower_draws <- as.matrix(draws_matrix[, lower_names, drop = FALSE])
upper_draws <- as.matrix(draws_matrix[, upper_names, drop = FALSE])

sigma_lower <- 0.001
sigma_upper <- 0.5
reference_scale <- 0.1

truncated_half_normal_logpdf <- function(x, scale, lower, upper) {
  log_density <- dnorm(x, mean = 0, sd = scale, log = TRUE)
  log_norm <- pnorm(upper, mean = 0, sd = scale, log.p = TRUE) -
    pnorm(lower, mean = 0, sd = scale, log.p = TRUE)
  log_density - log_norm
}

weighted_quantile <- function(x, w, probs) {
  ord <- order(x)
  x_ord <- x[ord]
  w_ord <- w[ord]
  cum_w <- cumsum(w_ord)
  total_w <- sum(w_ord)
  if (!is.finite(total_w) || total_w <= 0) {
    return(rep(NA_real_, length(probs)))
  }
  approx(
    x = cum_w / total_w,
    y = x_ord,
    xout = probs,
    method = "linear",
    ties = "ordered",
    rule = 2
  )$y
}

summarise_weighted_vector <- function(x, w) {
  mean_x <- sum(w * x)
  q <- weighted_quantile(x, w, c(0.05, 0.5, 0.95))
  tibble(
    mean = mean_x,
    median = q[[2]],
    q5 = q[[1]],
    q95 = q[[3]]
  )
}

stable_softmax_weights <- function(log_w) {
  centered <- log_w - max(log_w)
  w <- exp(centered)
  w / sum(w)
}

compute_scale_summary <- function(scale_value) {
  log_w <- truncated_half_normal_logpdf(
    x = sigma_draw,
    scale = scale_value,
    lower = sigma_lower,
    upper = sigma_upper
  ) - truncated_half_normal_logpdf(
    x = sigma_draw,
    scale = reference_scale,
    lower = sigma_lower,
    upper = sigma_upper
  )

  w <- stable_softmax_weights(log_w)
  ess <- 1 / sum(w^2)

  sigma_summary <- summarise_weighted_vector(sigma_draw, w)
  mean_abs_delta <- sum(w * rowMeans(abs(delta_draws)))
  expected_num_4n <- sum(colSums(prob_draws * w))
  assigned_4n_count <- sum(colSums(prob_draws * w) >= 0.5)

  tibble(
    prior_scale = scale_value,
    weight_ess = ess,
    weight_ess_fraction = ess / length(w),
    sigma_delta_dna_mean = sigma_summary$mean,
    sigma_delta_dna_median = sigma_summary$median,
    sigma_delta_dna_q5 = sigma_summary$q5,
    sigma_delta_dna_q95 = sigma_summary$q95,
    mean_abs_delta_dna = mean_abs_delta,
    expected_num_4n = expected_num_4n,
    assigned_4n_count = assigned_4n_count
  )
}

compute_sample_summary <- function(scale_value) {
  log_w <- truncated_half_normal_logpdf(
    x = sigma_draw,
    scale = scale_value,
    lower = sigma_lower,
    upper = sigma_upper
  ) - truncated_half_normal_logpdf(
    x = sigma_draw,
    scale = reference_scale,
    lower = sigma_lower,
    upper = sigma_upper
  )

  w <- stable_softmax_weights(log_w)

  prob_mean <- colSums(prob_draws * w)
  lower_mean <- colSums(lower_draws * w)
  upper_mean <- colSums(upper_draws * w)
  delta_summary <- bind_rows(lapply(seq_len(n_sample), function(idx) {
    summarise_weighted_vector(delta_draws[, idx], w)
  }))

  tibble(
    prior_scale = scale_value,
    sample_index = seq_len(n_sample),
    sample_name = input_tbl$sample_name,
    condition = input_tbl$condition,
    latest_match_date = input_tbl$latest_match_date,
    observed_lower_peak = input_tbl$modeled_peak_1,
    observed_upper_peak = input_tbl$modeled_peak_2,
    prob_4n_mean = prob_mean,
    assigned_state = ifelse(prob_mean >= 0.5, "4N", "2N"),
    mu_lower_mixture_mean = lower_mean,
    mu_upper_mixture_mean = upper_mean,
    delta_dna_mean = delta_summary$mean,
    delta_dna_median = delta_summary$median,
    delta_dna_q5 = delta_summary$q5,
    delta_dna_q95 = delta_summary$q95
  )
}

overall_tbl <- bind_rows(lapply(scale_grid, compute_scale_summary)) %>%
  arrange(prior_scale)

sample_tbl <- bind_rows(lapply(scale_grid, compute_sample_summary)) %>%
  arrange(prior_scale, sample_index)

reference_assignments <- sample_tbl %>%
  filter(abs(prior_scale - reference_scale) < 1e-12) %>%
  select(sample_name, reference_prob_4n = prob_4n_mean, reference_state = assigned_state)

sample_tbl <- sample_tbl %>%
  left_join(reference_assignments, by = "sample_name") %>%
  mutate(
    state_changed_vs_reference = assigned_state != reference_state,
    prob_4n_shift = prob_4n_mean - reference_prob_4n
  )

overall_tbl <- overall_tbl %>%
  left_join(
    sample_tbl %>%
      group_by(prior_scale) %>%
      summarise(
        n_state_changes_vs_reference = sum(state_changed_vs_reference),
        max_abs_prob_4n_shift = max(abs(prob_4n_shift)),
        .groups = "drop"
      ),
    by = "prior_scale"
  )

output_dir <- file.path(nuts_dir, "prior_sensitivity")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

overall_csv <- file.path(output_dir, "prior_sensitivity_overall.csv")
sample_csv <- file.path(output_dir, "prior_sensitivity_samples.csv")
plot_path <- file.path(output_dir, "prior_sensitivity_summary.png")

write.csv(overall_tbl, overall_csv, row.names = FALSE)
write.csv(sample_tbl, sample_csv, row.names = FALSE)

plot_tbl <- sample_tbl %>%
  mutate(
    sample_name = factor(sample_name, levels = rev(unique(input_tbl$sample_name))),
    prior_scale_label = sprintf("%.3f", prior_scale)
  )

overall_plot_tbl <- overall_tbl %>%
  mutate(prior_scale_label = sprintf("%.3f", prior_scale))

p1 <- ggplot(overall_plot_tbl, aes(x = prior_scale, y = sigma_delta_dna_median)) +
  geom_ribbon(aes(ymin = sigma_delta_dna_q5, ymax = sigma_delta_dna_q95), fill = "#b7d4ea") +
  geom_line(linewidth = 0.8, color = "#0b4f6c") +
  geom_point(size = 2, color = "#0b4f6c") +
  labs(
    title = "Reweighted posterior for sigma_delta_dna",
    x = "Half-normal prior scale on sigma_delta_dna",
    y = "Posterior median and 90% interval"
  ) +
  theme_minimal(base_size = 11)

p2 <- ggplot(overall_plot_tbl, aes(x = prior_scale, y = weight_ess_fraction)) +
  geom_line(linewidth = 0.8, color = "#8a3b12") +
  geom_point(size = 2, color = "#8a3b12") +
  geom_hline(yintercept = 0.1, linetype = "dashed", color = "grey50") +
  labs(
    title = "Importance-weight ESS fraction",
    x = "Half-normal prior scale on sigma_delta_dna",
    y = "ESS / draws"
  ) +
  theme_minimal(base_size = 11)

p3 <- ggplot(plot_tbl, aes(x = prior_scale, y = sample_name, fill = prob_4n_mean)) +
  geom_tile(color = "white", linewidth = 0.15) +
  scale_fill_gradient2(
    low = "#2166ac",
    mid = "white",
    high = "#b2182b",
    midpoint = 0.5,
    limits = c(0, 1)
  ) +
  labs(
    title = "Reweighted posterior P(4N) by sample",
    x = "Half-normal prior scale on sigma_delta_dna",
    y = NULL,
    fill = "P(4N)"
  ) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

png(filename = plot_path, width = 11, height = 13, units = "in", res = 180)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 1, heights = unit(c(3, 3, 7), "null"))))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
dev.off()

cat("Saved prior-sensitivity outputs under:", output_dir, "\n")
cat("Overall summary:", overall_csv, "\n")
cat("Sample summary:", sample_csv, "\n")
cat("Plot:", plot_path, "\n")
