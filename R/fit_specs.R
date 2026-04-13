fit_specs <- function(project_root = ".") {
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)

  list(
    first_pass = list(
      spec_name = "first_pass",
      model_variant = "first_pass",
      experiment_filter = "anoxia",
      exclude_sample_patterns = character(0),
      stan_file = file.path(project_root, "stan", "ploidy_histogram_first_pass.stan"),
      fit_root = file.path(project_root, "processed_data", "stan_optimize_fixed_init_anoxia"),
      plot_dir = file.path(project_root, "figure", "stan_optimize_fixed_init_anoxia"),
      exe_file = file.path(project_root, "processed_data", "stan_optimize_fixed_init_anoxia", "ploidy_histogram_first_pass_fixed_init.exe"),
      output_basename = "fixed_init_optimize",
      manifest_path = file.path(project_root, "processed_data", "stan_optimize_fixed_init_anoxia", "fixed_init_optimize_manifest.rds"),
      archive_dir = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia"),
      result_csv_path = file.path(project_root, "processed_data", "stan_optimize_fixed_init_anoxia", "fixed_init_optimize-1.csv"),
      report_csv_path = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia", "fits", "first_pass", "fixed_init_optimize-1.csv"),
      report_manifest_path = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia", "fits", "first_pass", "fixed_init_optimize_manifest.rds"),
      report_stan_data_path = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia", "inputs", "stan_data.Rds"),
      report_sample_metadata_path = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia", "inputs", "sample_metadata.csv"),
      checkpoint_dir = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia"),
      plot_label_prefix = "fixed_init",
      title_prefix = "Fixed-init optimization checkpoint overlay",
      include_group_effects = FALSE,
      cen_threshold = NULL
    ),
    dye_ratio = list(
      spec_name = "dye_ratio",
      model_variant = "original",
      experiment_filter = "anoxia",
      exclude_sample_patterns = character(0),
      stan_file = file.path(project_root, "stan", "ploidy_histogram_dye_ratio.stan"),
      fit_root = file.path(project_root, "processed_data", "stan_optimize_fixed_init_dye_ratio_anoxia"),
      plot_dir = file.path(project_root, "figure", "stan_optimize_fixed_init_dye_ratio_anoxia"),
      exe_file = file.path(project_root, "processed_data", "stan_optimize_fixed_init_dye_ratio_anoxia", "ploidy_histogram_dye_ratio_fixed_init.exe"),
      output_basename = "fixed_init_dye_ratio_optimize",
      manifest_path = file.path(project_root, "processed_data", "stan_optimize_fixed_init_dye_ratio_anoxia", "fixed_init_dye_ratio_optimize_manifest.rds"),
      archive_dir = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia"),
      result_csv_path = file.path(project_root, "processed_data", "stan_optimize_fixed_init_dye_ratio_anoxia", "fixed_init_dye_ratio_optimize-1.csv"),
      report_csv_path = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia", "fits", "dye_ratio", "fixed_init_dye_ratio_optimize-1.csv"),
      report_manifest_path = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia", "fits", "dye_ratio", "fixed_init_dye_ratio_optimize_manifest.rds"),
      report_stan_data_path = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia", "inputs", "stan_data.Rds"),
      report_sample_metadata_path = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia", "inputs", "sample_metadata.csv"),
      checkpoint_dir = NULL,
      plot_label_prefix = "fixed_init_dye_ratio",
      title_prefix = "Dye-ratio fixed-init optimization overlay",
      include_group_effects = FALSE,
      cen_threshold = NULL
    ),
    observed_cen_scaled = list(
      spec_name = "observed_cen_scaled",
      model_variant = "observed_cen_scaled",
      experiment_filter = "anoxia",
      exclude_sample_patterns = c("Sample_CEN", "2N_Control", "4N_Control"),
      stan_file = file.path(project_root, "stan", "ploidy_histogram_dye_ratio_observed_cen_scaled.stan"),
      fit_root = file.path(project_root, "processed_data", "stan_optimize_fixed_init_dye_ratio_observed_cen_scaled_anoxia_exclude_2n_4n_controls"),
      plot_dir = file.path(project_root, "figure", "stan_optimize_fixed_init_dye_ratio_observed_cen_scaled_anoxia_exclude_2n_4n_controls"),
      exe_file = file.path(project_root, "processed_data", "stan_optimize_fixed_init_dye_ratio_observed_cen_scaled_anoxia_exclude_2n_4n_controls", "ploidy_histogram_dye_ratio_observed_cen_scaled_fixed_init.exe"),
      output_basename = "fixed_init_dye_ratio_observed_cen_scaled_optimize",
      manifest_path = file.path(project_root, "processed_data", "stan_optimize_fixed_init_dye_ratio_observed_cen_scaled_anoxia_exclude_2n_4n_controls", "fixed_init_dye_ratio_observed_cen_scaled_optimize_manifest.rds"),
      archive_dir = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia"),
      result_csv_path = file.path(project_root, "processed_data", "stan_optimize_fixed_init_dye_ratio_observed_cen_scaled_anoxia_exclude_2n_4n_controls", "fixed_init_dye_ratio_observed_cen_scaled_optimize-1.csv"),
      report_csv_path = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia", "fits", "observed_cen_scaled", "fixed_init_dye_ratio_observed_cen_scaled_optimize-1.csv"),
      report_manifest_path = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia", "fits", "observed_cen_scaled", "fixed_init_dye_ratio_observed_cen_scaled_optimize_manifest.rds"),
      report_stan_data_path = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia", "inputs", "stan_data.Rds"),
      report_sample_metadata_path = file.path(project_root, "results", "2026-04-13-fixed-init-anoxia", "inputs", "sample_metadata.csv"),
      checkpoint_dir = NULL,
      plot_label_prefix = "fixed_init_dye_ratio_observed_cen_scaled",
      title_prefix = "Observed-CEN-scaled optimization overlay",
      include_group_effects = TRUE,
      cen_threshold = 20000
    )
  )
}

get_fit_spec <- function(spec_name, project_root = ".") {
  specs <- fit_specs(project_root)
  if (!spec_name %in% names(specs)) {
    stop(sprintf("Unknown fit spec '%s'.", spec_name))
  }
  specs[[spec_name]]
}

derive_sample_group <- function(sample_name) {
  sample_id <- tools::file_path_sans_ext(basename(sample_name))
  if (identical(sample_id, "Sample_CEN")) {
    return("CEN")
  }
  if (grepl("^Sample_SUM159_", sample_id)) {
    return(sub("^Sample_SUM159_([0-9]+_)?", "", sample_id))
  }
  sample_id
}

load_fit_stan_bundle <- function(project_root = ".", n_s_phase = 3L, stan_data_path = NULL) {
  project_root <- normalizePath(project_root, winslash = "/", mustWork = TRUE)
  if (is.null(stan_data_path)) {
    stan_data_path <- file.path(project_root, "processed_data", "stan_data.Rds")
  }
  stan_bundle <- readRDS(stan_data_path)
  if (is.null(stan_bundle$stan_data$N_s_phase)) {
    stan_bundle$stan_data$N_s_phase <- as.integer(n_s_phase)
  }
  stan_bundle
}

subset_stan_bundle_by_experiment <- function(stan_data, lut, experiment_pattern) {
  keep_experiments <- grepl(experiment_pattern, lut$experiment$experiment_name, ignore.case = TRUE)
  if (!any(keep_experiments)) {
    stop(sprintf("No experiments matched pattern '%s'.", experiment_pattern))
  }

  old_exp_ids <- lut$experiment$experiment_id[keep_experiments]
  keep_samples <- lut$sample$experiment_id %in% old_exp_ids
  old_sample_ids <- lut$sample$sample_id[keep_samples]
  keep_bins <- stan_data$sampleID %in% old_sample_ids

  exp_map <- setNames(seq_along(old_exp_ids), old_exp_ids)
  sample_map <- setNames(seq_along(old_sample_ids), old_sample_ids)

  experiment_lut <- lut$experiment[keep_experiments, , drop = FALSE]
  experiment_lut$experiment_id <- seq_len(nrow(experiment_lut))

  sample_lut <- lut$sample[keep_samples, , drop = FALSE]
  sample_lut$sample_id <- seq_len(nrow(sample_lut))
  sample_lut$experiment_id <- unname(exp_map[as.character(sample_lut$experiment_id)])

  stan_data$areaDNA_bin_center <- stan_data$areaDNA_bin_center[keep_bins]
  stan_data$count <- stan_data$count[keep_bins]
  stan_data$expID <- unname(exp_map[as.character(stan_data$expID[keep_bins])])
  stan_data$sampleID <- unname(sample_map[as.character(stan_data$sampleID[keep_bins])])
  stan_data$sampleExpID <- unname(exp_map[as.character(lut$sample$experiment_id[keep_samples])])
  stan_data$N <- length(stan_data$count)
  stan_data$N_exp <- nrow(experiment_lut)
  stan_data$N_sample <- nrow(sample_lut)

  list(stan_data = stan_data, lut = list(experiment = experiment_lut, sample = sample_lut))
}

exclude_samples_from_stan_bundle <- function(stan_data, lut, sample_patterns) {
  if (length(sample_patterns) == 0) {
    return(list(stan_data = stan_data, lut = lut, excluded_samples = lut$sample[FALSE, , drop = FALSE]))
  }

  sample_name <- tools::file_path_sans_ext(basename(lut$sample$sample_name))
  exclude_samples <- Reduce(`|`, lapply(sample_patterns, function(pattern) grepl(pattern, sample_name, ignore.case = TRUE)))

  if (!any(exclude_samples)) {
    return(list(stan_data = stan_data, lut = lut, excluded_samples = lut$sample[FALSE, , drop = FALSE]))
  }

  keep_samples <- !exclude_samples
  old_sample_ids <- lut$sample$sample_id[keep_samples]
  keep_bins <- stan_data$sampleID %in% old_sample_ids
  sample_map <- setNames(seq_along(old_sample_ids), old_sample_ids)

  excluded_samples <- lut$sample[exclude_samples, , drop = FALSE]
  sample_lut <- lut$sample[keep_samples, , drop = FALSE]
  sample_lut$sample_id <- seq_len(nrow(sample_lut))

  stan_data$areaDNA_bin_center <- stan_data$areaDNA_bin_center[keep_bins]
  stan_data$count <- stan_data$count[keep_bins]
  stan_data$expID <- stan_data$expID[keep_bins]
  stan_data$sampleID <- unname(sample_map[as.character(stan_data$sampleID[keep_bins])])
  stan_data$sampleExpID <- sample_lut$experiment_id
  stan_data$N <- length(stan_data$count)
  stan_data$N_sample <- nrow(sample_lut)

  list(stan_data = stan_data, lut = list(experiment = lut$experiment, sample = sample_lut), excluded_samples = excluded_samples)
}

augment_stan_bundle_with_groups <- function(stan_data, lut) {
  sample_groups <- vapply(lut$sample$sample_name, derive_sample_group, character(1))
  group_levels <- unique(sample_groups)
  lut$sample$sample_group <- sample_groups
  lut$sample$group_id <- match(sample_groups, group_levels)
  lut$group <- data.frame(group_id = seq_along(group_levels), group_name = group_levels, stringsAsFactors = FALSE)
  stan_data$N_group <- length(group_levels)
  stan_data$sampleGroupID <- as.integer(lut$sample$group_id)
  list(stan_data = stan_data, lut = lut)
}

prepare_fit_inputs <- function(spec, project_root = ".", stan_data_path = NULL) {
  stan_bundle <- load_fit_stan_bundle(project_root = project_root, stan_data_path = stan_data_path)
  prepared <- subset_stan_bundle_by_experiment(stan_bundle$stan_data, stan_bundle$lut, spec$experiment_filter)
  stan_data <- prepared$stan_data
  lut <- prepared$lut
  excluded_samples <- lut$sample[FALSE, , drop = FALSE]

  if (length(spec$exclude_sample_patterns) > 0) {
    prepared <- exclude_samples_from_stan_bundle(stan_data, lut, spec$exclude_sample_patterns)
    stan_data <- prepared$stan_data
    lut <- prepared$lut
    excluded_samples <- prepared$excluded_samples
  }

  if (isTRUE(spec$include_group_effects)) {
    prepared <- augment_stan_bundle_with_groups(stan_data, lut)
    stan_data <- prepared$stan_data
    lut <- prepared$lut
  }

  if (!is.null(spec$cen_threshold)) {
    stan_data$cen_threshold <- as.numeric(spec$cen_threshold)
  }

  list(stan_data = stan_data, lut = lut, excluded_samples = excluded_samples)
}
