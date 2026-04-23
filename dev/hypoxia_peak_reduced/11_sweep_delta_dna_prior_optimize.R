library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(grid)
source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1) args[[1]] else ".", winslash = "/", mustWork = TRUE)
spec_name <- if (length(args) >= 2) args[[2]] else "adjacent_8n_delta_dna"
iter <- if (length(args) >= 3) as.integer(args[[3]]) else 4000L
scale_grid <- if (length(args) >= 4) {
  as.numeric(strsplit(args[[4]], ",", fixed = TRUE)[[1]])
} else {
  c(0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4)
}

if (!is.finite(iter) || is.na(iter) || iter < 1) {
  stop("`iter` must be a positive integer.")
}
if (any(!is.finite(scale_grid)) || any(scale_grid <= 0)) {
  stop("Prior scales must be positive finite numbers.")
}

scale_grid <- sort(unique(scale_grid))
prepared <- readRDS(file.path(project_root, "dev", "hypoxia_peak_reduced", "output", "prepared_input.rds"))
spec_registry <- get_extension_registry()

if (!spec_name %in% names(spec_registry)) {
  stop(sprintf("Unknown extension `%s`.", spec_name))
}

spec <- spec_registry[[spec_name]]
if (!isTRUE(spec$sample_delta_dna)) {
  stop("This sweep only makes sense for models with `sample_delta_dna = TRUE`.")
}

format_scale_label <- function(x) {
  gsub("\\.", "p", formatC(x, format = "fg", digits = 3, flag = "#"))
}

sweep_root <- file.path(
  project_root,
  "dev",
  "hypoxia_peak_reduced",
  "output",
  "prior_scale_sweep",
  spec_name
)
dir.create(sweep_root, recursive = TRUE, showWarnings = FALSE)

run_one_scale <- function(scale_value) {
  scale_label <- format_scale_label(scale_value)
  spec_run <- spec
  spec_run$delta_dna_prior_scale <- scale_value
  spec_run$name <- paste0(spec_name, "_prior_", scale_label)

  run_parts <- run_ablation_fit(
    project_root = project_root,
    prepared = prepared,
    spec = spec_run,
    iter = iter,
    output_category = file.path("prior_scale_sweep", spec_name)
  )

  fit_tbl <- run_parts$fit_parts$fit_tbl %>%
    mutate(prior_scale = scale_value)
  param_tbl <- run_parts$fit_parts$param_tbl %>%
    mutate(
      prior_scale = scale_value,
      out_dir = run_parts$out_dir,
      objective_lp = run_parts$opt_row$lp__[[1]]
    )

  list(fit_tbl = fit_tbl, param_tbl = param_tbl)
}

results <- lapply(scale_grid, run_one_scale)
sample_tbl <- bind_rows(lapply(results, `[[`, "fit_tbl"))
param_tbl <- bind_rows(lapply(results, `[[`, "param_tbl"))

late_4n_pattern <- "^Sample_SUM159_4N_O[12]_A(12|19|22)$"

overall_tbl <- sample_tbl %>%
  group_by(prior_scale) %>%
  summarise(
    mean_phi1_counterfactual_peak = mean(phi1_counterfactual_peak, na.rm = TRUE),
    mean_fitted_g1 = mean(fitted_g1, na.rm = TRUE),
    mean_abs_delta_dna = mean(abs(delta_dna), na.rm = TRUE),
    n_assigned_4n = sum(assigned_state == "4N", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    param_tbl %>% select(prior_scale, sigma_delta_dna, p_4n, objective_lp),
    by = "prior_scale"
  )

target_tbl <- sample_tbl %>%
  filter(grepl(late_4n_pattern, sample_name))

overall_csv <- file.path(sweep_root, "prior_scale_sweep_overall.csv")
sample_csv <- file.path(sweep_root, "prior_scale_sweep_samples.csv")
target_csv <- file.path(sweep_root, "prior_scale_sweep_late_4N_O_samples.csv")
plot_path <- file.path(sweep_root, "prior_scale_sweep_summary.png")

write.csv(overall_tbl, overall_csv, row.names = FALSE)
write.csv(sample_tbl, sample_csv, row.names = FALSE)
write.csv(target_tbl, target_csv, row.names = FALSE)

p1 <- ggplot(overall_tbl, aes(x = prior_scale, y = sigma_delta_dna)) +
  geom_line(linewidth = 0.8, color = "#0b4f6c") +
  geom_point(size = 2, color = "#0b4f6c") +
  labs(
    title = "Optimized sigma_delta_dna versus prior scale",
    x = "Half-normal prior scale on sigma_delta_dna",
    y = "Optimized sigma_delta_dna"
  ) +
  theme_minimal(base_size = 11)

p2 <- ggplot(target_tbl, aes(x = prior_scale, y = phi1_counterfactual_peak, color = sample_name)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  labs(
    title = "Late 4N_O* phi1 counterfactual peak",
    x = "Half-normal prior scale on sigma_delta_dna",
    y = "phi1_counterfactual_peak",
    color = "Sample"
  ) +
  theme_minimal(base_size = 11)

p3 <- ggplot(target_tbl, aes(x = prior_scale, y = fitted_g1, color = sample_name)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  labs(
    title = "Late 4N_O* fitted lower peak",
    x = "Half-normal prior scale on sigma_delta_dna",
    y = "fitted_g1",
    color = "Sample"
  ) +
  theme_minimal(base_size = 11)

png(filename = plot_path, width = 10, height = 13, units = "in", res = 180)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 1, heights = unit(c(3, 5, 5), "null"))))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
dev.off()

cat("Saved prior scale sweep outputs under:", sweep_root, "\n")
cat("Overall summary:", overall_csv, "\n")
cat("Sample summary:", sample_csv, "\n")
cat("Late 4N_O summary:", target_csv, "\n")
cat("Plot:", plot_path, "\n")
