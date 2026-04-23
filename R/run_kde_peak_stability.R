args <- commandArgs(trailingOnly = TRUE)

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
source(file.path(project_root, "R", "kde_peak_stability.R"))

parse_cli_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected positional argument '%s'. Use --key value arguments.", key))
    }
    key <- sub("^--", "", key)
    if (key %in% c("no-plot", "help")) {
      out[[key]] <- TRUE
      i <- i + 1L
    } else {
      if (i == length(args)) {
        stop(sprintf("Missing value for --%s.", key))
      }
      out[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  out
}

print_usage <- function() {
  cat(
    paste(
      "Usage:",
      "  Rscript R/run_kde_peak_stability.R --output-root OUT_DIR [options]",
      "",
      "Inputs:",
      "  --dataset DATASET_ID       run one processed_data dataset folder",
      "  --datasets A,B,C           run multiple processed_data dataset folders",
      "  --export-path PATH         run one direct filtered_dna_area_vectors.rds export",
      "                             if neither is supplied, run every export under processed_data/",
      "",
      "Options:",
      "  --scale raw|transformed    DNA-A vector scale for input; default raw",
      "  --transform log10|raw      KDE scale; default log10",
      "  --bandwidth-method SJ|nrd0 base bandwidth estimator; default SJ",
      "  --global-bandwidth MODE    none|median-sample|pooled; default none",
      "  --bandwidths CSV           bandwidth multipliers; default 0.70,0.85,1.00,1.15,1.30",
      "  --density-n N              KDE grid size; default 2048",
      "  --min-relative-height X    candidate height cutoff; default 0.01",
      "  --min-relative-prominence X candidate prominence cutoff; default 0.008",
      "  --cluster-tolerance X      peak clustering tolerance on transform scale; default 0.035",
      "  --min-bandwidth-support N  bandwidths required for stable peak; default 4",
      "  --min-n N                  minimum finite positive observations; default 20",
      "  --plot-bins N              raw DNA histogram bins; default 180",
      "  --no-plot                  skip raw histogram PNGs",
      sep = "\n"
    ),
    "\n"
  )
}

parse_numeric_csv <- function(x, arg_name) {
  values <- as.numeric(trimws(strsplit(x, ",", fixed = TRUE)[[1]]))
  if (length(values) == 0L || any(!is.finite(values))) {
    stop(sprintf("--%s must be a comma-separated numeric list.", arg_name))
  }
  values
}

parse_positive_integer <- function(x, arg_name) {
  value <- as.integer(x)
  if (!is.finite(value) || value < 1L) {
    stop(sprintf("--%s must be a positive integer.", arg_name))
  }
  value
}

prefix_columns <- function(tbl, dataset_id, sample_name) {
  if (!nrow(tbl)) {
    return(tbl)
  }
  tbl$dataset_id <- dataset_id
  tbl$sample_name <- sample_name
  tbl[, c("dataset_id", "sample_name", setdiff(names(tbl), c("dataset_id", "sample_name"))), drop = FALSE]
}

empty_peaks_table <- function() {
  data.frame(
    dataset_id = character(),
    sample_name = character(),
    cluster_id = integer(),
    stable_rank = integer(),
    bandwidth_support = integer(),
    n_peak_calls = integer(),
    peak_fit_median = numeric(),
    peak_raw_median = numeric(),
    peak_raw_min = numeric(),
    peak_raw_max = numeric(),
    mean_relative_prominence = numeric(),
    stringsAsFactors = FALSE
  )
}

empty_diagnostics_table <- function() {
  data.frame(
    dataset_id = character(),
    sample_name = character(),
    bandwidth_multiplier = numeric(),
    bandwidth_method = character(),
    base_bandwidth = numeric(),
    bandwidth = numeric(),
    peak_fit = numeric(),
    peak_raw = numeric(),
    peak_height = numeric(),
    relative_height = numeric(),
    prominence = numeric(),
    relative_prominence = numeric(),
    stringsAsFactors = FALSE
  )
}

cli <- parse_cli_args(args)
if (isTRUE(cli$help) || length(args) == 0L) {
  print_usage()
  quit(status = if (isTRUE(cli$help)) 0L else 1L)
}

if (is.null(cli[["output-root"]])) {
  print_usage()
  stop("--output-root is required.")
}

scale <- cli[["scale"]] %||% "raw"
transform <- cli[["transform"]] %||% "log10"
bandwidth_method <- cli[["bandwidth-method"]] %||% "SJ"
global_bandwidth <- cli[["global-bandwidth"]] %||% "none"
bandwidth_multipliers <- parse_numeric_csv(cli[["bandwidths"]] %||% "0.70,0.85,1.00,1.15,1.30", "bandwidths")
density_n <- parse_positive_integer(cli[["density-n"]] %||% "2048", "density-n")
min_bandwidth_support <- parse_positive_integer(cli[["min-bandwidth-support"]] %||% "4", "min-bandwidth-support")
min_n <- parse_positive_integer(cli[["min-n"]] %||% "20", "min-n")
plot_bins <- parse_positive_integer(cli[["plot-bins"]] %||% "180", "plot-bins")
min_relative_height <- as.numeric(cli[["min-relative-height"]] %||% "0.01")
min_relative_prominence <- as.numeric(cli[["min-relative-prominence"]] %||% "0.008")
cluster_tolerance <- as.numeric(cli[["cluster-tolerance"]] %||% "0.035")

if (!scale %in% c("raw", "transformed")) {
  stop("--scale must be raw or transformed.")
}
if (!transform %in% c("log10", "raw")) {
  stop("--transform must be log10 or raw.")
}
if (!bandwidth_method %in% c("SJ", "nrd0")) {
  stop("--bandwidth-method must be SJ or nrd0.")
}
if (!global_bandwidth %in% c("none", "median-sample", "pooled")) {
  stop("--global-bandwidth must be none, median-sample, or pooled.")
}
if (!is.finite(min_relative_height) || min_relative_height < 0) {
  stop("--min-relative-height must be a non-negative number.")
}
if (!is.finite(min_relative_prominence) || min_relative_prominence < 0) {
  stop("--min-relative-prominence must be a non-negative number.")
}
if (!is.finite(cluster_tolerance) || cluster_tolerance <= 0) {
  stop("--cluster-tolerance must be a positive number.")
}
if (min_bandwidth_support > length(unique(bandwidth_multipliers))) {
  stop("--min-bandwidth-support cannot exceed the number of unique bandwidth multipliers.")
}

if (!is.null(cli[["export-path"]]) && (!is.null(cli[["dataset"]]) || !is.null(cli[["datasets"]]))) {
  stop("Use only one of --export-path, --dataset, or --datasets.")
}
if (!is.null(cli[["dataset"]]) && !is.null(cli[["datasets"]])) {
  stop("Use only one of --dataset or --datasets.")
}

if (!is.null(cli[["export-path"]])) {
  export_files <- normalizePath(cli[["export-path"]], winslash = "/", mustWork = TRUE)
} else if (!is.null(cli[["datasets"]])) {
  dataset_ids <- trimws(strsplit(cli[["datasets"]], ",", fixed = TRUE)[[1]])
  dataset_ids <- dataset_ids[nzchar(dataset_ids)]
  if (length(dataset_ids) == 0L) {
    stop("--datasets must contain at least one dataset id.")
  }
  export_files <- normalizePath(
    file.path(project_root, "processed_data", dataset_ids, "filtered_dna_area_vectors.rds"),
    winslash = "/",
    mustWork = TRUE
  )
} else if (!is.null(cli[["dataset"]])) {
  export_files <- normalizePath(
    file.path(project_root, "processed_data", cli[["dataset"]], "filtered_dna_area_vectors.rds"),
    winslash = "/",
    mustWork = TRUE
  )
} else {
  export_files <- find_filtered_dna_exports(file.path(project_root, "processed_data"))
}

if (length(export_files) == 0L) {
  stop("No filtered_dna_area_vectors.rds files found.")
}

output_root <- cli[["output-root"]]
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

all_dataset_summaries <- list()

for (export_path in export_files) {
  dataset_id <- dataset_id_from_export(export_path)
  payload <- read_filtered_dna_export(export_path)
  sample_names <- names(payload[[scale]])
  dataset_output_dir <- file.path(output_root, filesystem_safe(dataset_id))
  overlay_dir <- file.path(dataset_output_dir, "overlays")
  dir.create(overlay_dir, recursive = TRUE, showWarnings = FALSE)

  global_base_bandwidth <- NULL
  if (!identical(global_bandwidth, "none")) {
    message(sprintf("Estimating %s / %s bandwidth with mode %s", dataset_id, bandwidth_method, global_bandwidth))
    global_base_bandwidth <- estimate_global_kde_base_bandwidth(
      value_list = payload[[scale]],
      bandwidth_method = bandwidth_method,
      transform = transform,
      mode = global_bandwidth,
      min_n = min_n
    )
    message(sprintf("Using %s base bandwidth %.6g on the %s scale", dataset_id, global_base_bandwidth, transform))
  }

  dataset_peaks <- list()
  dataset_diagnostics <- list()
  dataset_summary <- list()

  for (sample_name in sample_names) {
    message(sprintf("Detecting KDE-stable peaks for %s / %s", dataset_id, sample_name))
    values <- payload[[scale]][[sample_name]]
    clean_values <- clean_positive_values(values, min_n = min_n)
    sample_base_bandwidth <- global_base_bandwidth %||% estimate_kde_base_bandwidth(
      values = values,
      bandwidth_method = bandwidth_method,
      transform = transform,
      min_n = min_n
    )

    result <- detect_kde_stable_peaks(
      values = values,
      bandwidth_multipliers = bandwidth_multipliers,
      bandwidth_method = bandwidth_method,
      base_bandwidth = sample_base_bandwidth,
      transform = transform,
      density_n = density_n,
      min_relative_height = min_relative_height,
      min_relative_prominence = min_relative_prominence,
      cluster_tolerance = cluster_tolerance,
      min_bandwidth_support = min_bandwidth_support,
      min_n = min_n
    )

    stable_peaks <- result$stable_peaks[result$stable_peaks$stable, , drop = FALSE]
    stable_peaks <- stable_peaks[, setdiff(names(stable_peaks), "stable"), drop = FALSE]
    if (nrow(stable_peaks)) {
      dataset_peaks[[length(dataset_peaks) + 1L]] <- prefix_columns(stable_peaks, dataset_id, sample_name)
    }
    if (nrow(result$peak_calls)) {
      dataset_diagnostics[[length(dataset_diagnostics) + 1L]] <- prefix_columns(result$peak_calls, dataset_id, sample_name)
    }

    if (!isTRUE(cli[["no-plot"]])) {
      plot_kde_stable_peak_histogram(
        values = values,
        stable_peaks = result$stable_peaks,
        output_path = file.path(overlay_dir, sprintf("%s_peak_overlay.png", filesystem_safe(sample_name))),
        sample_name = sample_name,
        dataset_id = dataset_id,
        bins = plot_bins
      )
    }

    dataset_summary[[length(dataset_summary) + 1L]] <- data.frame(
      dataset_id = dataset_id,
      sample_name = sample_name,
      scale = scale,
      transform = transform,
      bandwidth_method = bandwidth_method,
      global_bandwidth = global_bandwidth,
      base_bandwidth = sample_base_bandwidth,
      bandwidth_multipliers = paste(bandwidth_multipliers, collapse = ","),
      density_n = density_n,
      min_relative_height = min_relative_height,
      min_relative_prominence = min_relative_prominence,
      cluster_tolerance = cluster_tolerance,
      min_bandwidth_support = min_bandwidth_support,
      n_events = length(clean_values),
      n_peak_calls = nrow(result$peak_calls),
      n_peak_clusters = nrow(result$stable_peaks),
      n_stable_peaks = nrow(stable_peaks),
      stable_peak_raw_locations = paste(round(stable_peaks$peak_raw_median, 3), collapse = ";"),
      stringsAsFactors = FALSE
    )
  }

  peaks_tbl <- if (length(dataset_peaks)) do.call(rbind, dataset_peaks) else empty_peaks_table()
  diagnostics_tbl <- if (length(dataset_diagnostics)) do.call(rbind, dataset_diagnostics) else empty_diagnostics_table()
  summary_tbl <- do.call(rbind, dataset_summary)

  utils::write.csv(peaks_tbl, file.path(dataset_output_dir, "peaks.csv"), row.names = FALSE)
  utils::write.csv(diagnostics_tbl, file.path(dataset_output_dir, "diagnostics.csv"), row.names = FALSE)
  utils::write.csv(summary_tbl, file.path(dataset_output_dir, "sample_summary.csv"), row.names = FALSE)

  all_dataset_summaries[[length(all_dataset_summaries) + 1L]] <- summary_tbl
}

combined_summary <- if (length(all_dataset_summaries)) do.call(rbind, all_dataset_summaries) else data.frame()
utils::write.csv(combined_summary, file.path(output_root, "combined_sample_summary.csv"), row.names = FALSE)

cat(sprintf(
  "Wrote KDE peak stability outputs for %d samples across %d datasets to %s\n",
  nrow(combined_summary),
  length(export_files),
  normalizePath(output_root, winslash = "/", mustWork = TRUE)
))
