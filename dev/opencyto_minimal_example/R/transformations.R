split_transformed_dim <- function(dim_name, allowed_transforms) {
  suffix <- sub("^.*_", "", dim_name)

  if (!suffix %in% allowed_transforms) {
    return(list(base_dim = dim_name, transform = NA_character_, transformed = FALSE))
  }

  list(
    base_dim = sub(paste0("_", suffix, "$"), "", dim_name),
    transform = suffix,
    transformed = TRUE
  )
}

apply_allowed_transform <- function(x, transform, config) {
  if (transform == "ASINH") {
    return(asinh(x / config$transforms$ASINH$cofactor))
  }

  stop("Unsupported transform: ", transform)
}

preprocess_flow_frame_peacoqc <- function(fr, config) {
  peacoqc_config <- config$preprocessing
  missing_channels <- setdiff(peacoqc_config$channels, colnames(fr))

  if (length(missing_channels) > 0L) {
    stop("PeacoQC channels are missing from flowFrame: ", paste(missing_channels, collapse = ", "))
  }

  peacoqc_args <- list(
    ff = fr,
    channels = peacoqc_config$channels,
    determine_good_cells = peacoqc_config$determine_good_cells,
    plot = peacoqc_config$plot,
    save_fcs = peacoqc_config$save_fcs,
    report = peacoqc_config$report
  )
  peacoqc_result <- do.call(PeacoQC::PeacoQC, peacoqc_args)
  peacoqc_result$FinalFF
}

preprocess_flow_data <- function(fs, config) {
  preprocessing <- config$preprocessing
  if (is.null(preprocessing$enabled) || !isTRUE(preprocessing$enabled)) {
    return(fs)
  }

  if (preprocessing$method == "PeacoQC") {
    return(flowCore::fsApply(
      fs,
      function(fr) preprocess_flow_frame_peacoqc(fr, config),
      use.exprs = FALSE
    ))
  }

  stop("Unsupported preprocessing method: ", preprocessing$method)
}

add_transform_columns_to_flow_frame <- function(fr, dims, config) {
  expr <- flowCore::exprs(fr)
  new_cols <- list()

  for (dim_name in dims) {
    parsed <- split_transformed_dim(dim_name, config$allowed_transforms)
    if (!isTRUE(parsed$transformed)) {
      next
    }
    if (!parsed$base_dim %in% colnames(expr)) {
      stop("Cannot transform missing base channel: ", parsed$base_dim)
    }
    if (dim_name %in% colnames(expr)) {
      next
    }

    new_cols[[dim_name]] <- apply_allowed_transform(expr[, parsed$base_dim], parsed$transform, config)
  }

  if (length(new_cols) > 0L) {
    fr <- flowCore::fr_append_cols(fr, as.matrix(as.data.frame(new_cols, check.names = FALSE)))
  }

  fr
}

add_template_transform_columns <- function(fs, template_table, config) {
  dims <- unique(unlist(strsplit(template_table$dims[template_table$dims != ""], ",", fixed = TRUE)))
  dims <- trimws(dims)
  dims <- dims[nzchar(dims)]

  flowCore::fsApply(
    fs,
    function(fr) add_transform_columns_to_flow_frame(fr, dims, config),
    use.exprs = FALSE
  )
}
