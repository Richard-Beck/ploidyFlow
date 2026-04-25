args <- commandArgs(trailingOnly = TRUE)

project_root <- normalizePath(
  if (length(args) >= 1L) args[[1]] else ".",
  winslash = "/",
  mustWork = TRUE
)

job_table_path <- normalizePath(
  file.path(project_root, if (length(args) >= 2L) args[[2]] else "dev/opencyto_minimal_example/config/jobs.csv"),
  winslash = "/",
  mustWork = TRUE
)

example_dir <- file.path(project_root, "dev", "opencyto_minimal_example")

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggcyto)
  library(ggplot2)
  library(graph)
  library(grid)
  library(openCyto)
  library(patchwork)
  library(png)
})

source(file.path(example_dir, "R", "custom_gates.R"))
source(file.path(example_dir, "R", "channel_mapping.R"))
source(file.path(example_dir, "R", "transformations.R"))
source(file.path(example_dir, "R", "template_validation.R"))
source(file.path(example_dir, "R", "plotting.R"))

resolve_project_path <- function(project_root, path, base_dir = project_root, must_work = FALSE) {
  if (grepl("^[A-Za-z]:|^/", path)) {
    return(normalizePath(path, winslash = "/", mustWork = must_work))
  }
  normalizePath(file.path(base_dir, path), winslash = "/", mustWork = must_work)
}

read_dataset_config <- function(project_root, config_dir) {
  config_dir <- resolve_project_path(project_root, config_dir, must_work = TRUE)
  config <- jsonlite::read_json(file.path(config_dir, "config.json"), simplifyVector = TRUE)

  config$config_dir <- config_dir
  config$fcs_dir <- resolve_project_path(project_root, config$fcs_dir, must_work = TRUE)
  config$template_path <- resolve_project_path(project_root, config$template_path, base_dir = config_dir, must_work = TRUE)
  config$channel_map_path <- resolve_project_path(project_root, config$channel_map_path, base_dir = config_dir, must_work = TRUE)
  config$outputs$gated_dir <- resolve_project_path(project_root, config$outputs$gated_dir)
  config$outputs$debug_dir <- resolve_project_path(project_root, config$outputs$debug_dir)
  config$plots$output_dir <- resolve_project_path(project_root, config$plots$output_dir)
  config$plots$debug_dir <- config$outputs$debug_dir

  config
}

sample_summary_file <- function(sample_name) {
  paste0(sanitize_filename(sample_name), "_summary.png")
}

gate_range_for_dim <- function(gate, dim_name) {
  gate_slots <- methods::slotNames(gate)

  if ("boundaries" %in% gate_slots) {
    values <- gate@boundaries[, dim_name]
  } else if (all(c("min", "max") %in% gate_slots)) {
    values <- c(gate@min[[dim_name]], gate@max[[dim_name]])
  } else if (all(c("mean", "cov", "distance") %in% gate_slots)) {
    if (!dim_name %in% names(gate@mean)) {
      stop("Gate does not contain dimension '", dim_name, "'.")
    }
    if (!dim_name %in% rownames(gate@cov)) {
      stop("Gate covariance does not contain dimension '", dim_name, "'.")
    }

    radius <- sqrt(max(gate@distance, na.rm = TRUE) * gate@cov[dim_name, dim_name])
    values <- gate@mean[[dim_name]] + c(-radius, radius)
  } else {
    stop("Unsupported gate type for QC geometry: ", class(gate)[[1]])
  }

  range(values, finite = TRUE)
}

normalized_margin_distance <- function(gate_bounds, parent_bounds) {
  x_span <- parent_bounds$x_max - parent_bounds$x_min
  y_span <- parent_bounds$y_max - parent_bounds$y_min

  if (!is.finite(x_span) || !is.finite(y_span) || x_span <= 0 || y_span <= 0) {
    return(NA_real_)
  }

  min(
    abs(gate_bounds$x_min - parent_bounds$x_min) / x_span,
    abs(parent_bounds$x_max - gate_bounds$x_max) / x_span,
    abs(gate_bounds$y_min - parent_bounds$y_min) / y_span,
    abs(parent_bounds$y_max - gate_bounds$y_max) / y_span,
    na.rm = TRUE
  )
}

render_fsc_singlet_qc <- function(gs, margin_threshold = 0.02) {
  qc_rows <- list()
  population_paths <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")

  if (!"tumor_fsc_singlets" %in% basename(population_paths)) {
    return(data.frame())
  }

  for (sample_name in flowCore::sampleNames(gs)) {
    gh <- gs[[sample_name]]
    parent <- flowWorkspace::gh_pop_get_data(gh, "tumor")
    child <- flowWorkspace::gh_pop_get_data(gh, "tumor_fsc_singlets")
    gate <- flowWorkspace::gh_pop_get_gate(gh, "tumor_fsc_singlets")
    parent_expr <- flowCore::exprs(parent)

    required_dims <- c("FSC-A", "FSC-H")
    missing_dims <- setdiff(required_dims, colnames(parent_expr))
    if (length(missing_dims) > 0L) {
      stop("FSC singlet QC missing dimensions: ", paste(missing_dims, collapse = ", "))
    }

    parent_bounds <- list(
      x_min = min(parent_expr[, "FSC-A"], finite = TRUE),
      x_max = max(parent_expr[, "FSC-A"], finite = TRUE),
      y_min = min(parent_expr[, "FSC-H"], finite = TRUE),
      y_max = max(parent_expr[, "FSC-H"], finite = TRUE)
    )

    x_gate <- gate_range_for_dim(gate, "FSC-A")
    y_gate <- gate_range_for_dim(gate, "FSC-H")
    gate_bounds <- list(
      x_min = x_gate[[1]],
      x_max = x_gate[[2]],
      y_min = y_gate[[1]],
      y_max = y_gate[[2]]
    )

    tumor_n <- nrow(parent_expr)
    fsc_singlet_n <- nrow(flowCore::exprs(child))
    margin_distance <- normalized_margin_distance(gate_bounds, parent_bounds)
    fsc_gate_qc <- if (fsc_singlet_n == 0L) {
      "empty"
    } else if (is.finite(margin_distance) && margin_distance <= margin_threshold) {
      "near_margin"
    } else {
      "ok"
    }

    qc_rows[[length(qc_rows) + 1L]] <- data.frame(
      sample = sample_name,
      tumor_n = tumor_n,
      fsc_singlet_n = fsc_singlet_n,
      fsc_singlet_frac_of_tumor = fsc_singlet_n / tumor_n,
      fsc_gate_x_min = gate_bounds$x_min,
      fsc_gate_x_max = gate_bounds$x_max,
      fsc_gate_y_min = gate_bounds$y_min,
      fsc_gate_y_max = gate_bounds$y_max,
      parent_x_min = parent_bounds$x_min,
      parent_x_max = parent_bounds$x_max,
      parent_y_min = parent_bounds$y_min,
      parent_y_max = parent_bounds$y_max,
      fsc_gate_margin_distance = margin_distance,
      fsc_gate_qc = fsc_gate_qc,
      stringsAsFactors = FALSE
    )
  }

  do.call(rbind, qc_rows)
}

resolve_job_fcs_files <- function(project_root, config, fcs_files = "") {
  if (!is.null(fcs_files) && !is.na(fcs_files) && nzchar(trimws(fcs_files))) {
    files <- trimws(strsplit(fcs_files, ";", fixed = TRUE)[[1]])
    files <- files[nzchar(files)]
    return(sort(vapply(
      files,
      resolve_project_path,
      character(1),
      project_root = project_root,
      must_work = TRUE
    )))
  }

  fcs_files <- list.files(
    config$fcs_dir,
    pattern = config$fcs_pattern,
    full.names = TRUE,
    ignore.case = TRUE
  )
  sort(normalizePath(fcs_files, winslash = "/", mustWork = TRUE))
}

run_dataset_job <- function(project_root, config_dir, fcs_files = "") {
  config <- read_dataset_config(project_root, config_dir)
  set.seed(config$random_seed)

  if (dir.exists(config$outputs$gated_dir)) {
    unlink(config$outputs$gated_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(config$outputs$gated_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$outputs$debug_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(config$plots$output_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(config$plots$output_dir, full.names = TRUE), recursive = TRUE, force = TRUE)
  unlink(file.path(config$outputs$debug_dir, "panels"), recursive = TRUE, force = TRUE)

  fcs_files <- resolve_job_fcs_files(project_root, config, fcs_files)
  if (length(fcs_files) == 0L) {
    stop("No FCS files found for dataset '", config$dataset, "' in ", config$fcs_dir)
  }

  cat("\nRunning dataset: ", config$dataset, "\n", sep = "")
  cat("  FCS files: ", length(fcs_files), "\n", sep = "")
  cat("  Config: ", config$config_dir, "\n", sep = "")

  channel_map <- load_channel_map(config$channel_map_path)
  fs <- flowCore::read.flowSet(
    files = fcs_files,
    transformation = FALSE,
    truncate_max_range = FALSE
  )
  fs <- standardize_channels(fs, channel_map)
  fs <- preprocess_flow_data(fs, config)

  template_parts <- load_and_validate_template(config$template_path, fs, config)
  fs <- add_template_transform_columns(fs, template_parts$table, config)

  gs <- flowWorkspace::GatingSet(fs)
  openCyto::gt_gating(template_parts$gt, gs)

  stats <- flowWorkspace::gs_pop_get_stats(gs)
  stats_file <- file.path(config$outputs$debug_dir, "population_stats.csv")
  utils::write.csv(stats, stats_file, row.names = FALSE)

  qc_table <- render_fsc_singlet_qc(gs)
  if (nrow(qc_table) > 0L) {
    qc_table$dataset <- config$dataset
    qc_table <- qc_table[, c("dataset", setdiff(colnames(qc_table), "dataset"))]
    qc_file <- file.path(config$outputs$debug_dir, "qc_table.csv")
    utils::write.csv(qc_table, qc_file, row.names = FALSE)
  }

  flowWorkspace::save_gs(gs, path = config$outputs$gated_dir, backend_opt = "copy")

  figure_files <- character()
  sample_names <- flowCore::sampleNames(gs)
  for (sample_index in seq_along(sample_names)) {
    sample_plots <- config$plots
    sample_plots$summary_file <- sample_summary_file(sample_names[[sample_index]])
    plot_files <- render_configured_plots(
      gs,
      template_parts$gt,
      sample_plots,
      sample_index = sample_index
    )
    figure_files <- c(figure_files, plot_files$summary_file)
  }

  manifest <- data.frame(
    dataset = config$dataset,
    sample_name = sample_names,
    fcs_file = fcs_files,
    figure_file = figure_files,
    stringsAsFactors = FALSE
  )
  manifest_file <- file.path(config$outputs$debug_dir, "manifest.csv")
  utils::write.csv(manifest, manifest_file, row.names = FALSE)

  cat("  Gated flow dataset: ", config$outputs$gated_dir, "\n", sep = "")
  cat("  Combined figures: ", config$plots$output_dir, "\n", sep = "")
  cat("  Debug tables: ", config$outputs$debug_dir, "\n", sep = "")

  invisible(list(config = config, stats = stats, manifest = manifest))
}

register_custom_gates()

jobs <- utils::read.csv(job_table_path, stringsAsFactors = FALSE)
required_cols <- c("dataset", "config_dir")
missing_cols <- setdiff(required_cols, colnames(jobs))
if (length(missing_cols) > 0L) {
  stop("Job table is missing columns: ", paste(missing_cols, collapse = ", "))
}

for (i in seq_len(nrow(jobs))) {
  job_fcs_files <- if ("fcs_files" %in% colnames(jobs)) jobs$fcs_files[[i]] else ""
  run_dataset_job(project_root, jobs$config_dir[[i]], job_fcs_files)
}
