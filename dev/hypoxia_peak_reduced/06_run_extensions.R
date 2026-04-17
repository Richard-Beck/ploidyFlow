library(cmdstanr)
library(dplyr)
source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1) args[[1]] else ".", winslash = "/", mustWork = TRUE)
iter <- if (length(args) >= 2) as.integer(args[[2]]) else 4000
if (!is.finite(iter) || is.na(iter) || iter < 1) {
  stop("`iter` must be a positive integer.")
}
dev_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced")
out_dir <- file.path(dev_dir, "output")
extension_root <- file.path(out_dir, "extensions")
dir.create(extension_root, recursive = TRUE, showWarnings = FALSE)

prepared <- readRDS(file.path(out_dir, "prepared_input.rds"))
extension_registry <- get_extension_registry()

run_parts_list <- lapply(extension_registry, function(spec) {
  cat("Running extension:", spec$name, "\n")
  run_parts <- run_ablation_fit(
    project_root = project_root,
    prepared = prepared,
    spec = spec,
    iter = iter,
    output_category = "extensions"
  )
  write_fit_outputs(run_parts)
  render_fit_outputs(run_parts$out_dir)
  run_parts
})
names(run_parts_list) <- names(extension_registry)

comparison_tbl <- bind_rows(lapply(run_parts_list, summarise_ablation_run))
write.csv(comparison_tbl, file.path(extension_root, "extension_comparison.csv"), row.names = FALSE)

cat("Saved extension outputs under:", extension_root, "\n")
cat("Optimizer iterations per run:", iter, "\n")
print(comparison_tbl)
