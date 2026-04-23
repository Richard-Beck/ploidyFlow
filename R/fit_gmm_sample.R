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
    if (key %in% c("no-classification", "no-plot", "verbose", "help")) {
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
      "  Rscript R/fit_gmm_sample.R --dataset DATASET_ID --sample SAMPLE_NAME --output-dir OUT_DIR [options]",
      "",
      "Required unless --export-path is supplied:",
      "  --dataset DATASET_ID       processed_data dataset folder, e.g. hypoxia-sum159",
      "",
      "Alternative input:",
      "  --export-path PATH         direct filtered_dna_area_vectors.rds path",
      "",
      "Required:",
      "  --sample SAMPLE_NAME       sample name in the filtered DNA export",
      "  --output-dir OUT_DIR       output folder for gmm_fit.rds and CSV summaries",
      "",
      "Options:",
      "  --scale raw|transformed    DNA-A vector scale; default raw",
        "  --gmax N                   fit G = 1:N; default 16",
        "  --models E|V|E,V           variance models to search; default E",
        "  --criterion ICL|BIC        model-selection criterion; default ICL",
        "  --hc-subset N              deterministic max events for hcPairs; default 2000",
        "  --hc-model E|V             model used to build hcPairs; default E",
        "  --min-n N                  minimum finite observations; default 10",
        "  --no-classification        skip per-event classification.csv",
        "  --no-plot                  skip fit_plot.png",
        "  --verbose                  pass verbose=TRUE to mclustBIC",
      sep = "\n"
    ),
    "\n"
  )
}

safe_dataset_id_from_export <- function(export_path) {
  basename(dirname(normalizePath(export_path, winslash = "/", mustWork = TRUE)))
}

cli <- parse_cli_args(args)
if (isTRUE(cli$help) || length(args) == 0L) {
  print_usage()
  quit(status = if (isTRUE(cli$help)) 0L else 1L)
}

scale <- cli[["scale"]] %||% "raw"
gmax <- as.integer(cli[["gmax"]] %||% "16")
min_n <- as.integer(cli[["min-n"]] %||% "10")
hc_subset <- as.integer(cli[["hc-subset"]] %||% "2000")
hc_model <- cli[["hc-model"]] %||% "E"
model_names <- strsplit(cli[["models"]] %||% "E", ",", fixed = TRUE)[[1]]
model_names <- trimws(model_names)
criterion <- cli[["criterion"]] %||% "ICL"

if (is.null(cli[["sample"]]) || is.null(cli[["output-dir"]])) {
  print_usage()
  stop("--sample and --output-dir are required.")
}
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

if (!is.null(cli[["export-path"]])) {
  export_path <- cli[["export-path"]]
  dataset_id <- cli[["dataset"]] %||% safe_dataset_id_from_export(export_path)
} else {
  if (is.null(cli[["dataset"]])) {
    print_usage()
    stop("Either --dataset or --export-path is required.")
  }
  dataset_id <- cli[["dataset"]]
  export_path <- file.path(project_root, "processed_data", dataset_id, "filtered_dna_area_vectors.rds")
}

fit <- fit_filtered_dna_sample_gmm(
  export_path = export_path,
  sample_name = cli[["sample"]],
  scale = scale,
  G = seq_len(gmax),
  modelNames = model_names,
  criterion = criterion,
  hc_subset = hc_subset,
  hc_modelName = hc_model,
  min_n = min_n,
  verbose = isTRUE(cli$verbose)
)

written <- write_gmm_fit_outputs(
  fit = fit,
  output_dir = cli[["output-dir"]],
  sample_name = cli[["sample"]],
  dataset_id = dataset_id,
  scale = scale,
  include_classification = !isTRUE(cli[["no-classification"]]),
  include_plot = !isTRUE(cli[["no-plot"]])
)

summary_tbl <- gmm_fit_summary(fit, sample_name = cli[["sample"]], dataset_id = dataset_id, scale = scale)
cat(sprintf(
  "Selected %s model with G=%d by global %s %.3f across G=1:%d and modelNames=%s.\n",
  summary_tbl$modelName,
  summary_tbl$G,
  summary_tbl$selection_criterion,
  summary_tbl$selection_value,
  gmax,
  paste(model_names, collapse = ",")
))
cat(sprintf("Wrote outputs to %s\n", normalizePath(cli[["output-dir"]], winslash = "/", mustWork = TRUE)))
