library(cmdstanr)

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
stan_file <- file.path(project_root, "stan", "ploidy_histogram_first_pass.stan")
stan_data_path <- file.path(project_root, "processed_data", "stan_data.Rds")
fit_root <- file.path(project_root, "processed_data", "stan_fit_first_pass")
output_dir <- file.path(fit_root, "csv")
manifest_path <- file.path(fit_root, "fit_manifest.rds")

dir.create(fit_root, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stan_bundle <- readRDS(stan_data_path)
stan_data <- stan_bundle$stan_data
if (is.null(stan_data$N_s_phase)) {
  stan_data$N_s_phase <- 5L
}

mod <- cmdstan_model(stan_file, compile = TRUE)

fit <- mod$sample(
  data = stan_data,
  chains = 3,
  parallel_chains = 3,
  iter_warmup = 500,
  iter_sampling = 1000,
  refresh = 100,
  seed = 123,
  output_dir = output_dir
)

manifest <- list(
  generated_at = Sys.time(),
  stan_file = normalizePath(stan_file, winslash = "/", mustWork = TRUE),
  stan_data_path = normalizePath(stan_data_path, winslash = "/", mustWork = TRUE),
  output_dir = normalizePath(output_dir, winslash = "/", mustWork = TRUE),
  output_files = normalizePath(fit$output_files(), winslash = "/", mustWork = TRUE),
  fit_summary = fit$summary()
)

saveRDS(manifest, manifest_path)

cat(sprintf("Saved fit manifest to %s\n", manifest_path))
