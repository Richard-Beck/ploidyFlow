library(ggplot2)
source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1) args[[1]] else ".", winslash = "/", mustWork = TRUE)
spec_name <- if (length(args) >= 2) args[[2]] else "delta_dna"

out <- write_delta_only_counterfactual_outputs(
  project_root = project_root,
  spec_name = spec_name,
  output_category = "extensions"
)

cat("Saved delta-only counterfactual posterior outputs under:", out$nuts_dir, "\n")
cat("Posterior draws:", out$draw_rds_path, "\n")
cat("Posterior summary:", out$summary_csv_path, "\n")
cat("Posterior plot:", out$plot_path, "\n")
