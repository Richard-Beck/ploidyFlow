args <- commandArgs(trailingOnly = TRUE)

project_root <- normalizePath(".", winslash = "/", mustWork = TRUE)
source(file.path(project_root, "R", "gmm_fit.R"))

parse_cli_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected positional argument '%s'. Use --key value arguments.", key))
    }
    key <- sub("^--", "", key)
    if (key %in% c("classification", "no-plot", "verbose", "help")) {
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
      "  Rscript R/fit_gmm_all_samples.R --output-root OUT_DIR [options]",
      "",
      "Inputs:",
      "  --dataset DATASET_ID       fit one processed_data dataset folder",
      "  --datasets A,B,C           fit multiple processed_data dataset folders",
      "  --export-path PATH         fit one direct filtered_dna_area_vectors.rds export",
      "                             if neither is supplied, fit every export under processed_data/",
      "",
      "Options:",
      "  --scale raw|transformed    DNA-A vector scale; default raw",
      "  --gmax N                   fit G = 1:N; default 16",
      "  --models E|V|E,V           variance models to search; default E",
      "  --criterion ICL|BIC        model-selection criterion; default ICL",
      "  --hc-subset N              deterministic max events for hcPairs; default 2000",
      "  --hc-model E|V             model used to build hcPairs; default E",
      "  --min-n N                  minimum finite observations; default 10",
      "  --classification           write per-event classification.csv for every sample",
      "  --no-plot                  skip fit_plot.png",
      "  --verbose                  pass verbose=TRUE to mclustBIC",
      sep = "\n"
    ),
    "\n"
  )
}

filesystem_safe <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  gsub("(^_+|_+$)", "", x)
}

dataset_id_from_export <- function(export_path) {
  basename(dirname(normalizePath(export_path, winslash = "/", mustWork = TRUE)))
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
gmax <- as.integer(cli[["gmax"]] %||% "16")
min_n <- as.integer(cli[["min-n"]] %||% "10")
hc_subset <- as.integer(cli[["hc-subset"]] %||% "2000")
hc_model <- cli[["hc-model"]] %||% "E"
model_names <- strsplit(cli[["models"]] %||% "E", ",", fixed = TRUE)[[1]]
model_names <- trimws(model_names)
criterion <- cli[["criterion"]] %||% "ICL"

if (!is.finite(gmax) || gmax < 1L) {
  stop("--gmax must be a positive integer.")
}
if (!is.finite(min_n) || min_n < 1L) {
  stop("--min-n must be a positive integer.")
}
if (!is.finite(hc_subset) || hc_subset < 1L) {
  stop("--hc-subset must be a positive integer.")
}
if (!hc_model %in% c("E", "V")) {
  stop("--hc-model must be E or V.")
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
summary_parts <- list()

for (export_path in export_files) {
  dataset_id <- dataset_id_from_export(export_path)
  payload <- read_filtered_dna_export(export_path)
  sample_names <- names(payload[[scale]])

  for (sample_name in sample_names) {
    message(sprintf("Fitting %s / %s", dataset_id, sample_name))
    fit <- fit_univariate_gmm_bic(
      values = payload[[scale]][[sample_name]],
      G = seq_len(gmax),
      modelNames = model_names,
      criterion = criterion,
      hc_subset = hc_subset,
      hc_modelName = hc_model,
      min_n = min_n,
      verbose = isTRUE(cli$verbose)
    )
    fit$source <- list(
      export_path = export_path,
      sample_name = sample_name,
      scale = scale
    )

    sample_output_dir <- file.path(output_root, filesystem_safe(dataset_id), filesystem_safe(sample_name))
    write_gmm_fit_outputs(
      fit = fit,
      output_dir = sample_output_dir,
      sample_name = sample_name,
      dataset_id = dataset_id,
      scale = scale,
      include_classification = isTRUE(cli$classification),
      include_plot = !isTRUE(cli[["no-plot"]])
    )
    summary_parts[[length(summary_parts) + 1L]] <- gmm_fit_summary(
      fit = fit,
      sample_name = sample_name,
      dataset_id = dataset_id,
      scale = scale
    )
  }
}

summary_tbl <- do.call(rbind, summary_parts)
utils::write.csv(summary_tbl, file.path(output_root, "combined_summary.csv"), row.names = FALSE)
cat(sprintf("Wrote %d sample fits to %s\n", nrow(summary_tbl), normalizePath(output_root, winslash = "/", mustWork = TRUE)))
