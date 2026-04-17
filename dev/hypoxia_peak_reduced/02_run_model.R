library(cmdstanr)
source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

project_root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[[1]] else ".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output")
prepared <- readRDS(file.path(out_dir, "prepared_input.rds"))
baseline_spec <- get_ablation_registry()[["baseline"]]
run_parts <- run_ablation_fit(project_root = project_root, prepared = prepared, spec = baseline_spec, use_primary_output = TRUE)
write_fit_outputs(run_parts)

cat("Saved fit outputs under:", out_dir, "\n")
print(run_parts$fit_parts$param_tbl)
