interpolate_half_height_crossing <- function(x_left, y_left, x_right, y_right, target_y) {
  if (!all(is.finite(c(x_left, y_left, x_right, y_right, target_y)))) {
    return(NA_real_)
  }
  if (y_left == y_right) {
    return(mean(c(x_left, x_right)))
  }
  x_left + (target_y - y_left) * (x_right - x_left) / (y_right - y_left)
}

compute_peak_prominence <- function(x, y, peak_idx) {
  peak_y <- y[[peak_idx]]
  left_higher <- which(y[seq_len(peak_idx)] > peak_y)
  left_boundary <- if (length(left_higher) == 0) 1 else max(left_higher)
  right_higher_rel <- which(y[peak_idx:length(y)] > peak_y)
  right_boundary <- if (length(right_higher_rel) == 0) length(y) else peak_idx + min(right_higher_rel) - 1

  left_base <- min(y[left_boundary:peak_idx], na.rm = TRUE)
  right_base <- min(y[peak_idx:right_boundary], na.rm = TRUE)
  peak_y - max(left_base, right_base)
}

summarise_threshold_peak <- function(values, side = c("below", "above"), threshold = 12000, density_n = 1024) {
  side <- match.arg(side)
  values <- values[is.finite(values)]
  values <- if (identical(side, "below")) values[values < threshold] else values[values > threshold]

  if (length(values) < 10 || length(unique(values)) < 2) {
    return(tibble::tibble(
      peak_location = NA_real_,
      peak_height = NA_real_,
      half_height_left = NA_real_,
      half_height_right = NA_real_,
      peak_width_fwhm = NA_real_,
      peak_n = length(values)
    ))
  }

  density_obj <- stats::density(values, n = density_n, from = min(values), to = max(values), na.rm = TRUE)
  x <- density_obj$x
  y <- density_obj$y

  if (length(x) < 3 || all(!is.finite(y))) {
    return(tibble::tibble(
      peak_location = NA_real_,
      peak_height = NA_real_,
      half_height_left = NA_real_,
      half_height_right = NA_real_,
      peak_width_fwhm = NA_real_,
      peak_n = length(values)
    ))
  }

  local_max_idx <- which(y >= dplyr::lag(y, default = -Inf) & y >= dplyr::lead(y, default = -Inf))
  if (length(local_max_idx) == 0) {
    local_max_idx <- which.max(y)
  }
  peak_idx <- local_max_idx[[which.max(y[local_max_idx])]]
  peak_x <- x[[peak_idx]]
  peak_y <- y[[peak_idx]]
  half_height <- peak_y / 2

  left_idx <- which(y[seq_len(peak_idx)] < half_height)
  if (length(left_idx) == 0) {
    left_cross <- x[[1]]
  } else {
    left_low_idx <- max(left_idx)
    left_high_idx <- min(left_low_idx + 1, peak_idx)
    left_cross <- interpolate_half_height_crossing(
      x[[left_low_idx]], y[[left_low_idx]],
      x[[left_high_idx]], y[[left_high_idx]],
      half_height
    )
  }

  right_range <- peak_idx:length(y)
  right_idx_rel <- which(y[right_range] < half_height)
  if (length(right_idx_rel) == 0) {
    right_cross <- x[[length(x)]]
  } else {
    right_low_idx <- right_range[[min(right_idx_rel)]] - 1
    right_low_idx <- max(right_low_idx, peak_idx)
    right_high_idx <- min(right_low_idx + 1, length(x))
    right_cross <- interpolate_half_height_crossing(
      x[[right_low_idx]], y[[right_low_idx]],
      x[[right_high_idx]], y[[right_high_idx]],
      half_height
    )
  }

  tibble::tibble(
    peak_location = peak_x,
    peak_height = peak_y,
    half_height_left = left_cross,
    half_height_right = right_cross,
    peak_width_fwhm = right_cross - left_cross,
    peak_n = length(values)
  )
}

find_threshold_peak_candidates <- function(values, side = c("below", "above"), threshold = 12000, density_n = 2048) {
  side <- match.arg(side)
  values <- values[is.finite(values)]
  values <- if (identical(side, "below")) values[values < threshold] else values[values > threshold]

  if (length(values) < 10 || length(unique(values)) < 2) {
    return(tibble::tibble(
      candidate_rank = integer(),
      peak_location = numeric(),
      peak_height = numeric(),
      prominence = numeric(),
      half_height_left = numeric(),
      half_height_right = numeric(),
      peak_width_fwhm = numeric(),
      peak_n = integer()
    ))
  }

  density_obj <- stats::density(values, n = density_n, from = min(values), to = max(values), na.rm = TRUE)
  x <- density_obj$x
  y <- density_obj$y

  if (length(x) < 3 || all(!is.finite(y))) {
    return(tibble::tibble(
      candidate_rank = integer(),
      peak_location = numeric(),
      peak_height = numeric(),
      prominence = numeric(),
      half_height_left = numeric(),
      half_height_right = numeric(),
      peak_width_fwhm = numeric(),
      peak_n = integer()
    ))
  }

  local_max_idx <- which(y >= dplyr::lag(y, default = -Inf) & y >= dplyr::lead(y, default = -Inf))
  if (length(local_max_idx) == 0) {
    local_max_idx <- which.max(y)
  }

  cand_tbl <- dplyr::bind_rows(lapply(local_max_idx, function(peak_idx) {
    peak_x <- x[[peak_idx]]
    peak_y <- y[[peak_idx]]
    half_height <- peak_y / 2
    prominence <- compute_peak_prominence(x, y, peak_idx)

    left_idx <- which(y[seq_len(peak_idx)] < half_height)
    if (length(left_idx) == 0) {
      left_cross <- x[[1]]
    } else {
      left_low_idx <- max(left_idx)
      left_high_idx <- min(left_low_idx + 1, peak_idx)
      left_cross <- interpolate_half_height_crossing(
        x[[left_low_idx]], y[[left_low_idx]],
        x[[left_high_idx]], y[[left_high_idx]],
        half_height
      )
    }

    right_range <- peak_idx:length(y)
    right_idx_rel <- which(y[right_range] < half_height)
    if (length(right_idx_rel) == 0) {
      right_cross <- x[[length(x)]]
    } else {
      right_low_idx <- right_range[[min(right_idx_rel)]] - 1
      right_low_idx <- max(right_low_idx, peak_idx)
      right_high_idx <- min(right_low_idx + 1, length(x))
      right_cross <- interpolate_half_height_crossing(
        x[[right_low_idx]], y[[right_low_idx]],
        x[[right_high_idx]], y[[right_high_idx]],
        half_height
      )
    }

    tibble::tibble(
      peak_location = peak_x,
      peak_height = peak_y,
      prominence = prominence,
      half_height_left = left_cross,
      half_height_right = right_cross,
      peak_width_fwhm = right_cross - left_cross,
      peak_n = length(values)
    )
  }))

  cand_tbl %>%
    dplyr::arrange(dplyr::desc(prominence), dplyr::desc(peak_height), peak_location) %>%
    dplyr::mutate(candidate_rank = dplyr::row_number()) %>%
    dplyr::select(candidate_rank, dplyr::everything())
}

select_separated_peaks <- function(candidate_tbl, n_keep = 2, min_separation = 35000, min_prominence = 0) {
  if (nrow(candidate_tbl) == 0) {
    return(candidate_tbl)
  }

  filtered_tbl <- candidate_tbl %>%
    dplyr::filter(is.finite(prominence), prominence >= min_prominence) %>%
    dplyr::arrange(dplyr::desc(prominence), dplyr::desc(peak_height), peak_location)

  if (nrow(filtered_tbl) == 0) {
    return(filtered_tbl)
  }

  keep_idx <- integer()
  for (i in seq_len(nrow(filtered_tbl))) {
    loc_i <- filtered_tbl$peak_location[[i]]
    if (!length(keep_idx) || all(abs(loc_i - filtered_tbl$peak_location[keep_idx]) >= min_separation)) {
      keep_idx <- c(keep_idx, i)
    }
    if (length(keep_idx) >= n_keep) {
      break
    }
  }

  filtered_tbl[keep_idx, , drop = FALSE] %>%
    dplyr::arrange(peak_location) %>%
    dplyr::mutate(selected_rank = dplyr::row_number())
}

detect_sample_peaks <- function(
    values,
    threshold = 12000,
    tumor_peak_min_separation = 35000,
    tumor_peak_min_prominence = 0,
    density_n = 2048) {
  values <- as.numeric(values)
  values <- values[is.finite(values)]

  below_peak <- summarise_threshold_peak(values, side = "below", threshold = threshold, density_n = density_n)
  above_peak <- summarise_threshold_peak(values, side = "above", threshold = threshold, density_n = density_n)
  tumor_candidates <- find_threshold_peak_candidates(values, side = "above", threshold = threshold, density_n = density_n)
  tumor_selected <- select_separated_peaks(
    tumor_candidates,
    n_keep = 2,
    min_separation = tumor_peak_min_separation,
    min_prominence = tumor_peak_min_prominence
  )

  sum_below <- sum(values[values < threshold], na.rm = TRUE)
  sum_above <- sum(values[values > threshold], na.rm = TRUE)
  n_below <- sum(values < threshold, na.rm = TRUE)
  n_above <- sum(values > threshold, na.rm = TRUE)

  tibble::tibble(
    cen_peak = below_peak$peak_location,
    g1_peak = above_peak$peak_location,
    tumor_peak_1 = if (nrow(tumor_selected) >= 1) tumor_selected$peak_location[[1]] else NA_real_,
    tumor_peak_2 = if (nrow(tumor_selected) >= 2) tumor_selected$peak_location[[2]] else NA_real_,
    sum_below = sum_below,
    sum_above = sum_above,
    ratio_above_below = sum_above / pmax(sum_below, 1),
    n_below = n_below,
    n_above = n_above,
    count_ratio_above_below = n_above / pmax(n_below, 1),
    cen_peak_width_fwhm = below_peak$peak_width_fwhm,
    g1_peak_width_fwhm = above_peak$peak_width_fwhm,
    n_events = length(values)
  )
}

add_two_peak_model_columns <- function(input_tbl) {
  ratio_tbl <- input_tbl %>%
    dplyr::mutate(
      abs_log2_peak1_peak2_ratio = dplyr::if_else(
        is.finite(tumor_peak_1) & is.finite(tumor_peak_2) & tumor_peak_1 > 0 & tumor_peak_2 > 0,
        abs(log2(tumor_peak_1 / tumor_peak_2)),
        NA_real_
      )
    )

  max_abs_log2_ratio <- suppressWarnings(max(ratio_tbl$abs_log2_peak1_peak2_ratio, na.rm = TRUE))

  ratio_tbl %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      is_primary_outlier = is.finite(abs_log2_peak1_peak2_ratio) &&
        is.finite(max_abs_log2_ratio) &&
        abs_log2_peak1_peak2_ratio == max_abs_log2_ratio,
      dropped_outlier_peak = {
        candidate_peaks <- c(tumor_peak_1, tumor_peak_2)
        if (is_primary_outlier && sum(is.finite(candidate_peaks)) >= 2 && is.finite(g1_peak)) {
          finite_idx <- which(is.finite(candidate_peaks))
          candidate_peaks[[finite_idx[[which.max(abs(candidate_peaks[finite_idx] - g1_peak))]]]]
        } else {
          NA_real_
        }
      },
      modeled_peak_1 = {
        candidate_peaks <- c(tumor_peak_1, tumor_peak_2)
        if (is_primary_outlier && is.finite(dropped_outlier_peak)) {
          drop_idx <- which(is.finite(candidate_peaks) & candidate_peaks == dropped_outlier_peak)[[1]]
          candidate_peaks[[drop_idx]] <- NA_real_
        }
        modeled_peaks <- sort(candidate_peaks[is.finite(candidate_peaks)])
        if (length(modeled_peaks) >= 1) modeled_peaks[[1]] else NA_real_
      },
      modeled_peak_2 = {
        candidate_peaks <- c(tumor_peak_1, tumor_peak_2)
        if (is_primary_outlier && is.finite(dropped_outlier_peak)) {
          drop_idx <- which(is.finite(candidate_peaks) & candidate_peaks == dropped_outlier_peak)[[1]]
          candidate_peaks[[drop_idx]] <- NA_real_
        }
        modeled_peaks <- sort(candidate_peaks[is.finite(candidate_peaks)])
        if (length(modeled_peaks) >= 2) modeled_peaks[[2]] else NA_real_
      },
      modeled_peak_count = sum(is.finite(c(modeled_peak_1, modeled_peak_2)))
    ) %>%
    dplyr::ungroup()
}

detect_dataset_peaks <- function(
    dna_vectors,
    sample_metadata,
    dataset_id = NA_character_,
    dataset_name = NA_character_,
    threshold = 12000,
    tumor_peak_min_separation = 35000,
    tumor_peak_min_prominence = 0,
    density_n = 2048) {
  stopifnot(is.list(dna_vectors))

  sample_tbl <- tibble::tibble(
    sample_id = names(dna_vectors),
    raw_dna = unname(dna_vectors)
  )

  peak_tbl <- dplyr::bind_rows(lapply(seq_len(nrow(sample_tbl)), function(i) {
    detect_sample_peaks(
      values = sample_tbl$raw_dna[[i]],
      threshold = threshold,
      tumor_peak_min_separation = tumor_peak_min_separation,
      tumor_peak_min_prominence = tumor_peak_min_prominence,
      density_n = density_n
    ) %>%
      dplyr::mutate(
        sample_id = sample_tbl$sample_id[[i]],
        .before = 1
      )
  }))

  sample_metadata %>%
    dplyr::left_join(peak_tbl, by = "sample_id") %>%
    dplyr::mutate(
      dataset_id = dataset_id,
      dataset_name = dataset_name,
      .before = 1
    ) %>%
    add_two_peak_model_columns() %>%
    dplyr::arrange(sample_name)
}

build_peak_detection_overview <- function(peak_tbl) {
  peak_tbl %>%
    dplyr::summarise(
      dataset_id = dplyr::first(dataset_id),
      dataset_name = dplyr::first(dataset_name),
      n_samples = dplyr::n(),
      n_with_cen_peak = sum(is.finite(cen_peak)),
      n_with_g1_peak = sum(is.finite(g1_peak)),
      n_with_two_tumor_peaks = sum(is.finite(tumor_peak_1) & is.finite(tumor_peak_2)),
      n_primary_outliers = sum(is_primary_outlier, na.rm = TRUE),
      median_ratio_above_below = stats::median(ratio_above_below, na.rm = TRUE),
      median_count_ratio_above_below = stats::median(count_ratio_above_below, na.rm = TRUE)
    )
}

build_peak_plot_data <- function(peak_tbl, dna_vectors, max_signal = 200000, bins = 256) {
  hist_tbl <- dplyr::bind_rows(lapply(seq_along(dna_vectors), function(i) {
    sample_id <- names(dna_vectors)[[i]]
    values <- as.numeric(dna_vectors[[i]])
    values <- values[is.finite(values)]
    values <- values[values <= max_signal]

    if (!length(values)) {
      return(NULL)
    }

    hist_obj <- graphics::hist(values, breaks = bins, plot = FALSE)
    tibble::tibble(
      sample_id = sample_id,
      bin_center = hist_obj$mids,
      count = hist_obj$counts
    )
  }))

  peak_ann_tbl <- peak_tbl %>%
    dplyr::transmute(
      sample_id,
      sample_name,
      outlier_status = dplyr::if_else(is_primary_outlier, "Primary outlier", "Not primary outlier"),
      cen_peak,
      tumor_peak_1,
      tumor_peak_2
    )

  peak_line_tbl <- dplyr::bind_rows(
    peak_ann_tbl %>%
      dplyr::transmute(sample_id, sample_name, outlier_status, peak_type = "CEN peak", peak_location = cen_peak),
    peak_ann_tbl %>%
      dplyr::transmute(sample_id, sample_name, outlier_status, peak_type = "Tumor peak 1", peak_location = tumor_peak_1),
    peak_ann_tbl %>%
      dplyr::transmute(sample_id, sample_name, outlier_status, peak_type = "Tumor peak 2", peak_location = tumor_peak_2)
  ) %>%
    dplyr::filter(is.finite(peak_location), peak_location <= max_signal)

  list(
    hist_tbl = hist_tbl,
    peak_line_tbl = peak_line_tbl
  )
}

parse_peak_sample_metadata <- function(sample_names) {
  sample_id <- tools::file_path_sans_ext(basename(sample_names))
  cell_line <- rep(NA_character_, length(sample_id))
  condition <- rep(NA_character_, length(sample_id))

  cell_line[sample_id == "Sample_CEN"] <- "CEN"

  sum159_idx <- grepl("^Sample_SUM159_", sample_id)
  cell_line[sum159_idx] <- "SUM159"
  condition[sum159_idx] <- sub("^Sample_SUM159_([0-9]+_)?", "", sample_id[sum159_idx])

  nonsample_idx <- !sum159_idx & !grepl("^Sample_", sample_id)
  nonsample_parts <- strsplit(sample_id[nonsample_idx], "_", fixed = TRUE)
  cell_line[nonsample_idx] <- vapply(nonsample_parts, function(x) if (length(x)) x[[1]] else NA_character_, character(1))
  condition[nonsample_idx] <- vapply(
    nonsample_parts,
    function(x) if (length(x) >= 2) paste(x[-1], collapse = "_") else NA_character_,
    character(1)
  )

  tibble::tibble(
    sample_name = basename(sample_names),
    sample_id = sample_id,
    cell_line = cell_line,
    condition = condition
  )
}

infer_dna_channels <- function(flow_frame) {
  channel_names <- colnames(flow_frame)

  preferred_area <- c("450-Violet C-A", "450/50 Violet B-A")
  preferred_height <- c("450-Violet C-H", "450/50 Violet B-H")
  for (i in seq_along(preferred_area)) {
    if (preferred_area[[i]] %in% channel_names && preferred_height[[i]] %in% channel_names) {
      return(c(preferred_area[[i]], preferred_height[[i]]))
    }
  }

  area_candidates <- grep("450.*Violet.*-A$", channel_names, value = TRUE)
  if (!length(area_candidates)) {
    area_candidates <- grep("DNA.*-A$", channel_names, value = TRUE, ignore.case = TRUE)
  }

  for (area_name in area_candidates) {
    height_name <- sub("-A$", "-H", area_name)
    if (height_name %in% channel_names) {
      return(c(area_name, height_name))
    }
  }

  stop("Could not infer a DNA area/height channel pair from the FCS file.")
}

detect_dna_singlets <- function(values_area, values_height, nmads = 4) {
  xy <- cbind(values_area, values_height)
  ok <- is.finite(xy[, 1]) & is.finite(xy[, 2]) & xy[, 1] > 0 & xy[, 2] > 0

  if (sum(ok) < 100) {
    return(rep(FALSE, nrow(xy)))
  }

  fit <- MASS::rlm(xy[ok, 2] ~ xy[ok, 1])
  pred <- stats::coef(fit)[1] + stats::coef(fit)[2] * xy[ok, 1]
  res <- xy[ok, 2] - pred

  spread <- stats::mad(res, na.rm = TRUE)
  if (!is.finite(spread) || spread == 0) {
    spread <- stats::sd(res, na.rm = TRUE)
  }
  if (!is.finite(spread) || spread == 0) {
    spread <- 1e-8
  }

  keep_ok <- abs(res - stats::median(res, na.rm = TRUE)) <= nmads * spread
  keep <- rep(FALSE, nrow(xy))
  keep[which(ok)[keep_ok]] <- TRUE
  keep
}

build_peak_detection_cache <- function(dataset_dir, cache_dir, nmads = 4) {
  flow_files <- list.files(
    dataset_dir,
    pattern = "\\.fcs$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  flow_files <- flow_files[!startsWith(basename(flow_files), "._")]
  flow_files <- sort(normalizePath(flow_files, winslash = "/", mustWork = TRUE))

  if (!length(flow_files)) {
    stop(sprintf("No FCS files found under '%s'.", dataset_dir))
  }

  first_ff <- flowCore::read.FCS(flow_files[[1]], transformation = FALSE, truncate_max_range = FALSE)
  dna_channels <- infer_dna_channels(first_ff)

  raw_vectors <- setNames(lapply(flow_files, function(flow_path) {
    ff <- flowCore::read.FCS(flow_path, transformation = FALSE, truncate_max_range = FALSE)
    expr <- flowCore::exprs(ff)
    keep <- detect_dna_singlets(expr[, dna_channels[[1]]], expr[, dna_channels[[2]]], nmads = nmads)
    vals <- as.numeric(expr[keep, dna_channels[[1]]])
    vals[is.finite(vals)]
  }), tools::file_path_sans_ext(basename(flow_files)))

  metadata_tbl <- parse_peak_sample_metadata(flow_files)

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(list(raw = raw_vectors), file.path(cache_dir, "filtered_dna_area_vectors.rds"))
  write.csv(metadata_tbl, file.path(cache_dir, "sample_metadata.csv"), row.names = FALSE)

  list(
    dna_vectors = list(raw = raw_vectors),
    metadata_tbl = metadata_tbl,
    dna_channels = dna_channels
  )
}
