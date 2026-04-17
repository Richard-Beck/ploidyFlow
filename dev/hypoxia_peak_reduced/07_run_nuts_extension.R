library(cmdstanr)
source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1) args[[1]] else ".", winslash = "/", mustWork = TRUE)
spec_name <- if (length(args) >= 2) args[[2]] else "delta_dna"
chains <- if (length(args) >= 3) as.integer(args[[3]]) else 4L
iter_warmup <- if (length(args) >= 4) as.integer(args[[4]]) else 1000L
iter_sampling <- if (length(args) >= 5) as.integer(args[[5]]) else 1000L
adapt_delta <- if (length(args) >= 6) as.numeric(args[[6]]) else 0.9
max_treedepth <- if (length(args) >= 7) as.integer(args[[7]]) else 12L

if (!is.finite(chains) || is.na(chains) || chains < 1) {
  stop("`chains` must be a positive integer.")
}
if (!is.finite(iter_warmup) || is.na(iter_warmup) || iter_warmup < 0) {
  stop("`iter_warmup` must be a non-negative integer.")
}
if (!is.finite(iter_sampling) || is.na(iter_sampling) || iter_sampling < 1) {
  stop("`iter_sampling` must be a positive integer.")
}
if (!is.finite(adapt_delta) || is.na(adapt_delta) || adapt_delta <= 0 || adapt_delta >= 1) {
  stop("`adapt_delta` must be between 0 and 1.")
}
if (!is.finite(max_treedepth) || is.na(max_treedepth) || max_treedepth < 1) {
  stop("`max_treedepth` must be a positive integer.")
}

out_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output")
prepared <- readRDS(file.path(out_dir, "prepared_input.rds"))
extension_registry <- get_extension_registry()

if (!spec_name %in% names(extension_registry)) {
  stop(sprintf(
    "Unknown extension `%s`. Available extensions: %s",
    spec_name,
    paste(names(extension_registry), collapse = ", ")
  ))
}

run_parts <- run_ablation_nuts(
  project_root = project_root,
  prepared = prepared,
  spec = extension_registry[[spec_name]],
  chains = chains,
  parallel_chains = chains,
  iter_warmup = iter_warmup,
  iter_sampling = iter_sampling,
  adapt_delta = adapt_delta,
  max_treedepth = max_treedepth
)

cat("Saved NUTS outputs under:", run_parts$nuts_out_dir, "\n")
cat("CmdStan CSV files:\n")
cat(paste(" -", run_parts$metadata$stan_csv_files), sep = "\n")
cat("\nCombined draws matrix:", file.path(run_parts$nuts_out_dir, "draws_matrix.rds"), "\n")
cat("Draw summary:", file.path(run_parts$nuts_out_dir, "draws_summary.csv"), "\n")
