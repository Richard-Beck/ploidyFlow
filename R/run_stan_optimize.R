library(cmdstanr)

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
stan_file <- file.path(project_root, "stan", "ploidy_histogram_first_pass.stan")
stan_data_path <- file.path(project_root, "processed_data", "stan_data.Rds")
fit_root <- file.path(project_root, "processed_data", "stan_optimize_first_pass")
seed <- as.integer(Sys.getenv("PLOIDYFLOW_STAN_SEED", "123"))
run_id <- sprintf("seed_%d", seed)
manifest_path <- file.path(fit_root, sprintf("optimize_manifest_%s.rds", run_id))

dir.create(fit_root, recursive = TRUE, showWarnings = FALSE)

stan_bundle <- readRDS(stan_data_path)
stan_data <- stan_bundle$stan_data
if (is.null(stan_data$N_s_phase)) {
  stan_data$N_s_phase <- 5L
}

mod <- cmdstan_model(stan_file, compile = TRUE)

fit <- mod$optimize(
  data = stan_data,
  seed = seed,
  output_dir = fit_root,
  output_basename = sprintf("optimize_output_%s", run_id)
)

manifest <- list(
  generated_at = Sys.time(),
  seed = seed,
  stan_file = normalizePath(stan_file, winslash = "/", mustWork = TRUE),
  stan_data_path = normalizePath(stan_data_path, winslash = "/", mustWork = TRUE),
  output_file = normalizePath(fit$output_files(), winslash = "/", mustWork = TRUE),
  metadata = fit$metadata(),
  mle = fit$draws()
)

saveRDS(manifest, manifest_path)

cat(sprintf("Saved optimize manifest to %s\n", manifest_path))
