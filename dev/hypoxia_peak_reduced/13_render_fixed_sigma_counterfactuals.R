library(dplyr)
library(tibble)
source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1) args[[1]] else ".", winslash = "/", mustWork = TRUE)

dev_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced")
sweep_root <- file.path(dev_dir, "output", "fixed_sigma_sweep", "adjacent_8n_delta_dna")
prepared <- readRDS(file.path(dev_dir, "output", "prepared_input.rds"))
spec <- get_extension_registry()[["adjacent_8n_delta_dna"]]

if (!dir.exists(sweep_root)) {
  stop("Fixed-sigma sweep directory not found: ", sweep_root)
}

sigma_dirs <- list.dirs(sweep_root, recursive = FALSE, full.names = TRUE)
sigma_dirs <- sigma_dirs[grepl("sigma_", basename(sigma_dirs))]

if (!length(sigma_dirs)) {
  stop("No sigma-specific fit directories found under: ", sweep_root)
}

render_one_dir <- function(out_dir) {
  opt_files <- list.files(
    out_dir,
    pattern = "^hypoxia_peak_reduced_adjacent_8n_delta_dna_fixed_sigma_opt-.*\\.csv$",
    full.names = TRUE
  )
  if (!length(opt_files)) {
    stop("No optimizer CSV found in ", out_dir)
  }

  opt_path <- opt_files[[1]]
  opt_row <- read.csv(opt_path, comment.char = "#", check.names = FALSE)[1, , drop = FALSE]
  sigma_fixed <- sub("^sigma_", "", basename(out_dir))
  sigma_fixed <- as.numeric(gsub("p", ".", sigma_fixed, fixed = TRUE))

  fit_parts <- compute_fit_table_from_spec(
    input_tbl = prepared$input_tbl,
    opt_row = opt_row,
    date_levels = prepared$date_levels,
    spec = spec
  )

  fit_tbl <- fit_parts$fit_tbl %>%
    mutate(sigma_delta_dna_fixed = sigma_fixed)
  date_tbl <- fit_parts$date_tbl
  param_tbl <- fit_parts$param_tbl %>%
    mutate(
      sigma_delta_dna_fixed = sigma_fixed,
      objective_lp = if ("lp__" %in% names(opt_row)) opt_row$lp__[[1]] else NA_real_
    )

  write.csv(fit_tbl, file.path(out_dir, "sample_predictions.csv"), row.names = FALSE)
  write.csv(date_tbl, file.path(out_dir, "date_phi.csv"), row.names = FALSE)
  write.csv(param_tbl, file.path(out_dir, "fit_parameters.csv"), row.names = FALSE)

  summary_tbl <- render_fit_outputs(out_dir)

  tibble(
    sigma_dir = basename(out_dir),
    sigma_delta_dna_fixed = sigma_fixed,
    phi_counterfactual_png = normalizePath(file.path(out_dir, "phi_counterfactual.png"), winslash = "/", mustWork = TRUE),
    sample_predictions_csv = normalizePath(file.path(out_dir, "sample_predictions.csv"), winslash = "/", mustWork = TRUE),
    sample_peak_summary_csv = normalizePath(file.path(out_dir, "sample_peak_summary.csv"), winslash = "/", mustWork = TRUE),
    n_assigned_4n = sum(summary_tbl$assigned_state == "4N", na.rm = TRUE)
  )
}

manifest_tbl <- bind_rows(lapply(sigma_dirs, render_one_dir)) %>%
  arrange(sigma_delta_dna_fixed)

manifest_path <- file.path(sweep_root, "fixed_sigma_counterfactual_manifest.csv")
write.csv(manifest_tbl, manifest_path, row.names = FALSE)

cat("Rendered full counterfactual figures for", nrow(manifest_tbl), "fixed-sigma fits.\n")
cat("Manifest:", manifest_path, "\n")
