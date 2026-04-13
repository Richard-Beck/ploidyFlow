args <- commandArgs(trailingOnly = TRUE)

project_root <- normalizePath(
  if (length(args) >= 1) args[[1]] else ".",
  winslash = "/",
  mustWork = TRUE
)
processed_root <- file.path(project_root, "processed_data")
output_path <- file.path(processed_root, "stan_data.Rds")
n_bins <- 512L
n_s_phase <- 5L

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

find_export_files <- function(root_dir) {
  files <- list.files(
    root_dir,
    pattern = "^filtered_dna_area_vectors\\.rds$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  sort(normalizePath(files, winslash = "/", mustWork = TRUE))
}

normalize_export_payload <- function(payload, export_path) {
  if (!is.list(payload) || is.null(payload$raw)) {
    stop(sprintf(
      "Expected '%s' to contain a list with at least a named '$raw' entry. Rerun workflow.Rmd first.",
      export_path
    ))
  }

  if (is.null(names(payload$raw)) || anyDuplicated(names(payload$raw))) {
    stop(sprintf(
      "Invalid sample names in '$raw' for '%s'.",
      export_path
    ))
  }

  transformed <- payload$transformed
  if (is.null(transformed)) {
    transformed <- payload$raw
  }

  if (!identical(names(payload$raw), names(transformed))) {
    stop(sprintf(
      "Raw and transformed sample names do not align in '%s'.",
      export_path
    ))
  }

  list(
    raw = payload$raw,
    transformed = transformed
  )
}

coerce_sample_values <- function(raw_values, transformed_values, sample_name) {
  raw_values <- as.numeric(raw_values)
  transformed_values <- as.numeric(transformed_values)

  keep <- is.finite(raw_values) & is.finite(transformed_values)
  raw_values <- raw_values[keep]

  if (length(raw_values) == 0L) {
    stop(sprintf("Sample '%s' has no finite raw/transformed values after filtering.", sample_name))
  }

  raw_values
}

export_files <- find_export_files(processed_root)
if (length(export_files) == 0) {
  stop("No filtered_dna_area_vectors.rds files found under processed_data/.")
}

experiment_lut <- data.frame(
  experiment_id = seq_along(export_files),
  experiment_name = basename(dirname(export_files)),
  export_path = export_files,
  stringsAsFactors = FALSE
)

sample_lut_parts <- vector("list", length(export_files))
sample_parts <- vector("list", length(export_files))
sample_id_counter <- 1L

for (i in seq_along(export_files)) {
  experiment_id <- experiment_lut$experiment_id[[i]]
  experiment_name <- experiment_lut$experiment_name[[i]]
  export_payload <- normalize_export_payload(readRDS(export_files[[i]]), export_files[[i]])

  sample_names <- names(export_payload$raw)
  sample_ids <- seq.int(sample_id_counter, length.out = length(sample_names))
  sample_id_counter <- max(sample_ids) + 1L

  sample_lut_parts[[i]] <- data.frame(
    sample_id = sample_ids,
    sample_name = sample_names,
    experiment_id = experiment_id,
    experiment_name = experiment_name,
    stringsAsFactors = FALSE
  )

  sample_parts[[i]] <- Map(
    f = function(raw_values, transformed_values, sample_id, sample_name) {
      clean_raw <- coerce_sample_values(raw_values, transformed_values, sample_name)
      data.frame(
        areaDNA = clean_raw,
        expID = rep.int(experiment_id, length(clean_raw)),
        sampleID = rep.int(sample_id, length(clean_raw)),
        sample_name = rep.int(sample_name, length(clean_raw)),
        stringsAsFactors = FALSE
      )
    },
    raw_values = export_payload$raw,
    transformed_values = export_payload$transformed,
    sample_id = sample_ids,
    sample_name = sample_names
  )
}

sample_lut <- do.call(rbind, sample_lut_parts)
sample_lut$sample_group <- vapply(sample_lut$sample_name, derive_sample_group, character(1))
group_levels <- unique(sample_lut$sample_group)
group_lut <- data.frame(
  group_id = seq_along(group_levels),
  group_name = group_levels,
  stringsAsFactors = FALSE
)
sample_lut$group_id <- match(sample_lut$sample_group, group_lut$group_name)
all_samples <- unlist(sample_parts, recursive = FALSE)
combined <- do.call(rbind, all_samples)

global_min <- min(combined$areaDNA)
global_max <- max(combined$areaDNA)

if (!is.finite(global_min) || !is.finite(global_max)) {
  stop("Global raw DNA area range is not finite.")
}

if (global_min == global_max) {
  stop("Global raw DNA area range has zero width; cannot construct shared bins.")
}

global_breaks <- seq(global_min, global_max, length.out = n_bins + 1L)
global_centers <- (head(global_breaks, -1L) + tail(global_breaks, -1L)) / 2

bin_parts <- lapply(split(combined, combined$sampleID), function(sample_df) {
  hist_out <- hist(
    sample_df$areaDNA,
    breaks = global_breaks,
    plot = FALSE,
    include.lowest = TRUE,
    right = TRUE
  )

  data.frame(
    areaDNA_bin_center = global_centers,
    count = as.integer(hist_out$counts),
    expID = rep.int(sample_df$expID[[1]], n_bins),
    sampleID = rep.int(sample_df$sampleID[[1]], n_bins),
    stringsAsFactors = FALSE
  )
})

binned <- do.call(rbind, bin_parts)
rownames(binned) <- NULL

stan_data <- list(
  N = nrow(binned),
  N_exp = nrow(experiment_lut),
  N_sample = nrow(sample_lut),
  N_group = nrow(group_lut),
  N_bin = n_bins,
  N_s_phase = n_s_phase,
  areaDNA_bin_center = as.numeric(binned$areaDNA_bin_center),
  count = as.integer(binned$count),
  sampleExpID = as.integer(sample_lut$experiment_id),
  sampleGroupID = as.integer(sample_lut$group_id),
  expID = as.integer(binned$expID),
  sampleID = as.integer(binned$sampleID)
)

lut <- list(
  experiment = experiment_lut,
  sample = sample_lut,
  group = group_lut
)

dir.create(processed_root, recursive = TRUE, showWarnings = FALSE)
saveRDS(
  list(
    stan_data = stan_data,
    lut = lut
  ),
  file = output_path
)

cat(sprintf("Wrote %s\n", output_path))
