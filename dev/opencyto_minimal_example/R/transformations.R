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

template_preprocessing_channels <- function(template_table, config, always_include = c("DNA_AREA", "DNA_HEIGHT")) {
  dims <- unique(unlist(strsplit(template_table$dims[template_table$dims != ""], ",", fixed = TRUE)))
  dims <- trimws(dims)
  dims <- dims[nzchar(dims)]

  base_dims <- vapply(
    dims,
    function(dim_name) split_transformed_dim(dim_name, config$allowed_transforms)$base_dim,
    character(1)
  )

  unique(c(base_dims, always_include))
}

standard_to_raw_channels <- function(standard_channels, channel_map) {
  missing_channels <- setdiff(standard_channels, channel_map$standard_channel)
  if (length(missing_channels) > 0L) {
    stop("Channel map is missing standard channels: ", paste(missing_channels, collapse = ", "))
  }

  channel_map$raw_channel[match(standard_channels, channel_map$standard_channel)]
}

remove_margins_flow_frame <- function(fr, channels) {
  missing_channels <- setdiff(channels, colnames(fr))
  if (length(missing_channels) > 0L) {
    stop("RemoveMargins channels are missing from flowFrame: ", paste(missing_channels, collapse = ", "))
  }

  n_before <- nrow(flowCore::exprs(fr))
  margin_result <- PeacoQC::RemoveMargins(
    ff = fr,
    channels = channels,
    output = "frame"
  )
  n_after <- nrow(flowCore::exprs(margin_result))

  list(
    frame = margin_result,
    metrics = data.frame(
      margin_events_before = n_before,
      margin_events_after = n_after,
      margin_removed_events = n_before - n_after,
      margin_removed_percent = if (n_before > 0L) 100 * (n_before - n_after) / n_before else NA_real_,
      margin_channels = paste(channels, collapse = ";"),
      stringsAsFactors = FALSE
    )
  )
}

remove_margins_flow_data <- function(fs, channels) {
  sample_names <- flowCore::sampleNames(fs)
  frames <- vector("list", length(sample_names))
  metrics <- vector("list", length(sample_names))

  for (sample_index in seq_along(sample_names)) {
    sample_name <- sample_names[[sample_index]]
    parts <- remove_margins_flow_frame(fs[[sample_name]], channels)
    flowCore::identifier(parts$frame) <- sample_name
    frames[[sample_index]] <- parts$frame
    metrics[[sample_index]] <- cbind(
      data.frame(sample = sample_name, stringsAsFactors = FALSE),
      parts$metrics
    )
  }

  names(frames) <- sample_names

  list(
    fs = flowCore::flowSet(frames),
    qc = do.call(rbind, metrics)
  )
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

  n_before <- nrow(flowCore::exprs(fr))
  n_after <- nrow(flowCore::exprs(peacoqc_result$FinalFF))
  increasing_channels <- unlist(peacoqc_result$WeirdChannels$Increasing, use.names = FALSE)
  decreasing_channels <- unlist(peacoqc_result$WeirdChannels$Decreasing, use.names = FALSE)

  list(
    frame = peacoqc_result$FinalFF,
    metrics = data.frame(
      peacoqc_events_before = n_before,
      peacoqc_events_after = n_after,
      peacoqc_removed_events = n_before - n_after,
      peacoqc_removed_percent = as.numeric(peacoqc_result$PercentageRemoved),
      peacoqc_it_removed_percent = as.numeric(peacoqc_result$ITPercentage),
      peacoqc_mad_removed_percent = as.numeric(peacoqc_result$MADPercentage),
      peacoqc_consecutive_removed_percent = as.numeric(peacoqc_result$ConsecutiveCellsPercentage),
      peacoqc_channels = paste(peacoqc_config$channels, collapse = ";"),
      peacoqc_trend_flag = peacoqc_result$WeirdChannels$Changing_channel,
      peacoqc_increasing_channels = paste(increasing_channels, collapse = ";"),
      peacoqc_decreasing_channels = paste(decreasing_channels, collapse = ";"),
      stringsAsFactors = FALSE
    )
  )
}

preprocess_flow_data <- function(fs, config) {
  preprocessing <- config$preprocessing
  if (is.null(preprocessing$enabled) || !isTRUE(preprocessing$enabled)) {
    return(list(fs = fs, qc = data.frame()))
  }

  if (preprocessing$method == "PeacoQC") {
    sample_names <- flowCore::sampleNames(fs)
    frames <- vector("list", length(sample_names))
    metrics <- vector("list", length(sample_names))

    for (sample_index in seq_along(sample_names)) {
      sample_name <- sample_names[[sample_index]]
      parts <- preprocess_flow_frame_peacoqc(fs[[sample_name]], config)
      flowCore::identifier(parts$frame) <- sample_name
      frames[[sample_index]] <- parts$frame
      metrics[[sample_index]] <- cbind(
        data.frame(sample = sample_name, stringsAsFactors = FALSE),
        parts$metrics
      )
    }

    names(frames) <- sample_names

    return(list(
      fs = flowCore::flowSet(frames),
      qc = do.call(rbind, metrics)
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
