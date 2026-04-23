library(ggplot2)
source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1) args[[1]] else ".", winslash = "/", mustWork = TRUE)
spec_name <- if (length(args) >= 2) args[[2]] else "adjacent_8n_delta_dna"
output_category <- if (length(args) >= 3) args[[3]] else "extensions"

out_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output", output_category, spec_name, "nuts", "flow_only_distribution")
distribution_rds_path <- file.path(out_dir, "flow_only_distribution.rds")
if (!file.exists(distribution_rds_path)) {
  stop(sprintf("Missing distribution object: %s", distribution_rds_path))
}

distribution_object <- readRDS(distribution_rds_path)
state_parts <- build_flow_state_distribution_outputs(distribution_object)
grid_csv_path <- file.path(out_dir, "g0g1_ploidy_density_grid.csv")
write.csv(state_parts$g0g1_density_grid_tbl, grid_csv_path, row.names = FALSE)

cat("Saved compact G0/G1 density grid to:", grid_csv_path, "\n")
