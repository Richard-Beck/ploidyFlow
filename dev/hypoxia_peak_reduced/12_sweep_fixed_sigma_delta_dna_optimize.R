library(cmdstanr)
library(ggplot2)
library(dplyr)
library(tibble)
library(grid)
source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1) args[[1]] else ".", winslash = "/", mustWork = TRUE)
iter <- if (length(args) >= 2) as.integer(args[[2]]) else 4000L
sigma_grid <- if (length(args) >= 3) {
  as.numeric(strsplit(args[[3]], ",", fixed = TRUE)[[1]])
} else {
  c(0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2)
}

if (!is.finite(iter) || is.na(iter) || iter < 1) {
  stop("`iter` must be a positive integer.")
}
if (any(!is.finite(sigma_grid)) || any(sigma_grid < 0.001) || any(sigma_grid > 0.5)) {
  stop("All fixed sigma values must lie in [0.001, 0.5].")
}

sigma_grid <- sort(unique(sigma_grid))
prepared <- readRDS(file.path(project_root, "dev", "hypoxia_peak_reduced", "output", "prepared_input.rds"))
input_tbl <- prepared$input_tbl
date_levels <- prepared$date_levels
dev_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced")
sweep_root <- file.path(dev_dir, "output", "fixed_sigma_sweep", "adjacent_8n_delta_dna")
dir.create(sweep_root, recursive = TRUE, showWarnings = FALSE)

build_fixed_sigma_data <- function(prepared, sigma_value) {
  input_tbl <- prepared$input_tbl
  peak_upper <- input_tbl$modeled_peak_2
  list(
    N = nrow(input_tbl),
    N_date = length(prepared$date_levels),
    date_id = input_tbl$date_id,
    x_ratio = input_tbl$ratio_scaled,
    y_cen = input_tbl$cen_peak,
    y_peak_lower = input_tbl$modeled_peak_1,
    has_upper_peak = as.integer(is.finite(peak_upper)),
    y_peak_upper = ifelse(is.finite(peak_upper), peak_upper, input_tbl$modeled_peak_1),
    sigma_delta_dna_fixed = sigma_value
  )
}

build_fixed_sigma_init <- function(input_tbl, date_levels, previous_opt_row = NULL) {
  spec <- get_extension_registry()[["adjacent_8n_delta_dna"]]
  init <- build_ablation_init(input_tbl, date_levels, spec)
  init$sigma_delta_dna <- NULL

  if (is.null(previous_opt_row)) {
    return(init)
  }

  shared_names <- intersect(names(previous_opt_row), names(init))
  for (nm in shared_names) {
    init[[nm]] <- previous_opt_row[[nm]][[1]]
  }

  vector_prefixes <- c("log_phi_raw", "delta_dna")
  for (prefix in vector_prefixes) {
    cols <- grep(paste0("^", prefix, "\\."), names(previous_opt_row), value = TRUE)
    if (length(cols)) {
      init[[prefix]] <- as.numeric(previous_opt_row[1, cols, drop = TRUE])
    }
  }

  init
}

compute_fixed_sigma_fit <- function(opt_row, sigma_value) {
  spec <- get_extension_registry()[["adjacent_8n_delta_dna"]]
  fit_parts <- compute_fit_table_from_spec(
    input_tbl = input_tbl,
    opt_row = opt_row,
    date_levels = date_levels,
    spec = spec
  )

  fit_tbl <- fit_parts$fit_tbl %>%
    mutate(sigma_delta_dna_fixed = sigma_value)
  param_tbl <- fit_parts$param_tbl %>%
    mutate(
      sigma_delta_dna_fixed = sigma_value,
      objective_lp = opt_row$lp__[[1]]
    )

  list(fit_tbl = fit_tbl, param_tbl = param_tbl)
}

model_path <- file.path(dev_dir, "hypoxia_peak_reduced_adjacent_8n_delta_dna_fixed_sigma.stan")
exe_path <- file.path(sweep_root, "hypoxia_peak_reduced_adjacent_8n_delta_dna_fixed_sigma.exe")
mod <- cmdstanr::cmdstan_model(
  stan_file = model_path,
  exe_file = exe_path,
  compile = TRUE,
  force_recompile = FALSE
)

all_fit_tbl <- list()
all_param_tbl <- list()
previous_opt_row <- NULL

format_sigma_label <- function(x) {
  gsub("\\.", "p", formatC(x, format = "fg", digits = 3, flag = "#"))
}

for (i in seq_along(sigma_grid)) {
  sigma_value <- sigma_grid[[i]]
  sigma_label <- format_sigma_label(sigma_value)
  scale_dir <- file.path(sweep_root, paste0("sigma_", sigma_label))
  dir.create(scale_dir, recursive = TRUE, showWarnings = FALSE)

  stan_data <- build_fixed_sigma_data(prepared, sigma_value)
  init <- build_fixed_sigma_init(input_tbl, date_levels, previous_opt_row = previous_opt_row)

  fit <- mod$optimize(
    data = stan_data,
    init = list(init),
    seed = 123,
    iter = iter,
    algorithm = "lbfgs",
    output_dir = scale_dir,
    output_basename = "hypoxia_peak_reduced_adjacent_8n_delta_dna_fixed_sigma_opt"
  )

  opt_row <- read.csv(fit$output_files(), comment.char = "#", check.names = FALSE)[1, , drop = FALSE]
  parts <- compute_fixed_sigma_fit(opt_row, sigma_value)
  all_fit_tbl[[i]] <- parts$fit_tbl
  all_param_tbl[[i]] <- parts$param_tbl
  previous_opt_row <- opt_row
}

sample_tbl <- bind_rows(all_fit_tbl)
param_tbl <- bind_rows(all_param_tbl)
target_tbl <- sample_tbl %>%
  filter(grepl("^Sample_SUM159_4N_O[12]_A(12|19|22)$", sample_name))

overall_tbl <- sample_tbl %>%
  group_by(sigma_delta_dna_fixed) %>%
  summarise(
    mean_phi1_counterfactual_peak = mean(phi1_counterfactual_peak, na.rm = TRUE),
    mean_fitted_g1 = mean(fitted_g1, na.rm = TRUE),
    mean_abs_delta_dna = mean(abs(delta_dna), na.rm = TRUE),
    n_assigned_4n = sum(assigned_state == "4N", na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    param_tbl %>%
      select(sigma_delta_dna_fixed, p_4n, objective_lp),
    by = "sigma_delta_dna_fixed"
  )

overall_csv <- file.path(sweep_root, "fixed_sigma_sweep_overall.csv")
sample_csv <- file.path(sweep_root, "fixed_sigma_sweep_samples.csv")
target_csv <- file.path(sweep_root, "fixed_sigma_sweep_late_4N_O_samples.csv")
plot_path <- file.path(sweep_root, "fixed_sigma_sweep_summary.png")

write.csv(overall_tbl, overall_csv, row.names = FALSE)
write.csv(sample_tbl, sample_csv, row.names = FALSE)
write.csv(target_tbl, target_csv, row.names = FALSE)

p1 <- ggplot(target_tbl, aes(x = sigma_delta_dna_fixed, y = phi1_counterfactual_peak, color = sample_name)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  labs(
    title = "Late 4N_O* phi1 counterfactual peak",
    x = "Fixed sigma_delta_dna",
    y = "phi1_counterfactual_peak",
    color = "Sample"
  ) +
  theme_minimal(base_size = 11)

p2 <- ggplot(target_tbl, aes(x = sigma_delta_dna_fixed, y = fitted_g1, color = sample_name)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  labs(
    title = "Late 4N_O* fitted lower peak",
    x = "Fixed sigma_delta_dna",
    y = "fitted_g1",
    color = "Sample"
  ) +
  theme_minimal(base_size = 11)

p3 <- ggplot(target_tbl, aes(x = sigma_delta_dna_fixed, y = delta_dna, color = sample_name)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.6) +
  labs(
    title = "Late 4N_O* optimized delta_dna",
    x = "Fixed sigma_delta_dna",
    y = "delta_dna",
    color = "Sample"
  ) +
  theme_minimal(base_size = 11)

png(filename = plot_path, width = 10, height = 13, units = "in", res = 180)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 1, heights = unit(c(5, 5, 5), "null"))))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(p3, vp = viewport(layout.pos.row = 3, layout.pos.col = 1))
dev.off()

cat("Saved fixed-sigma sweep outputs under:", sweep_root, "\n")
cat("Overall summary:", overall_csv, "\n")
cat("Sample summary:", sample_csv, "\n")
cat("Late 4N_O summary:", target_csv, "\n")
cat("Plot:", plot_path, "\n")
