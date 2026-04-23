library(ggplot2)
source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1) args[[1]] else ".", winslash = "/", mustWork = TRUE)
spec_name <- if (length(args) >= 2) args[[2]] else "adjacent_8n_delta_dna"
max_draws <- if (length(args) >= 3) as.integer(args[[3]]) else NA_integer_
if (!is.finite(max_draws)) {
  max_draws <- NULL
}

out <- write_flow_posterior_histogram_outputs(
  project_root = project_root,
  spec_name = spec_name,
  output_category = "extensions",
  max_draws = max_draws
)

cat("Saved flow-only posterior distribution outputs under:", out$out_dir, "\n")
cat("Distribution object:", out$distribution_rds_path, "\n")
cat("Sample summary:", out$sample_summary_csv_path, "\n")
cat("Component summary:", out$component_summary_csv_path, "\n")
cat("Overlay plot:", out$overlay_plot_path, "\n")
