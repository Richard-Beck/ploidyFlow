library(ggplot2)
source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1) args[[1]] else ".", winslash = "/", mustWork = TRUE)
spec_name <- if (length(args) >= 2) args[[2]] else "adjacent_8n_delta_dna"
sample_arg <- if (length(args) >= 3) args[[3]] else "Sample_SUM159_2N_O1_A6,Sample_SUM159_4N_O1_A6"
draw_arg <- if (length(args) >= 4) args[[4]] else "1,2,3"

sample_names <- strsplit(sample_arg, ",", fixed = TRUE)[[1]]
sample_names <- trimws(sample_names)
sample_names <- sample_names[nzchar(sample_names)]

draw_ids <- as.integer(trimws(strsplit(draw_arg, ",", fixed = TRUE)[[1]]))
draw_ids <- draw_ids[is.finite(draw_ids)]

out <- write_flow_posterior_histogram_trial_outputs(
  project_root = project_root,
  spec_name = spec_name,
  output_category = "extensions",
  sample_names = sample_names,
  draw_ids = draw_ids
)

cat("Saved trial flow-only outputs under:", out$out_dir, "\n")
cat("Overlay plot:", out$overlay_plot_path, "\n")
cat("Draw overlay table:", out$draw_overlay_csv_path, "\n")
cat("Component draws:", out$component_draw_csv_path, "\n")
