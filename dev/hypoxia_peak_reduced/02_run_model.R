library(cmdstanr)
source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

project_root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[[1]] else ".", winslash = "/", mustWork = TRUE)
dev_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced")
out_dir <- file.path(dev_dir, "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

prepared <- readRDS(file.path(out_dir, "prepared_input.rds"))
input_tbl <- prepared$input_tbl
stan_data <- prepared$stan_data
date_levels <- prepared$date_levels

init <- list(
  alpha = 0,
  beta = 0.6,
  log_M_cen = log(max(input_tbl$cen_peak) * 2),
  log_R_2n = log(8),
  log_R_4n_over_2n = log(2),
  log_rho = 0,
  log_phi_raw = rep(0, length(date_levels)),
  sigma_cen = 0.08,
  sigma_g1 = 0.08,
  p_4n = mean(grepl("^4N", input_tbl$condition))
)

mod <- cmdstan_model(
  stan_file = file.path(dev_dir, "hypoxia_peak_reduced.stan"),
  exe_file = file.path(out_dir, "hypoxia_peak_reduced.exe"),
  compile = TRUE,
  force_recompile = FALSE
)

fit <- mod$optimize(
  data = stan_data,
  init = list(init),
  seed = 123,
  iter = 4000,
  algorithm = "lbfgs",
  output_dir = out_dir,
  output_basename = "hypoxia_peak_reduced_opt"
)

opt_row <- read.csv(fit$output_files(), comment.char = "#", check.names = FALSE)[1, , drop = FALSE]
fit_parts <- compute_fit_table(input_tbl = input_tbl, opt_row = opt_row, date_levels = date_levels)

write.csv(fit_parts$fit_tbl, file.path(out_dir, "sample_predictions.csv"), row.names = FALSE)
write.csv(fit_parts$date_tbl, file.path(out_dir, "date_phi.csv"), row.names = FALSE)
write.csv(fit_parts$param_tbl, file.path(out_dir, "fit_parameters.csv"), row.names = FALSE)

cat("Saved fit outputs under:", out_dir, "\n")
print(fit_parts$param_tbl)
