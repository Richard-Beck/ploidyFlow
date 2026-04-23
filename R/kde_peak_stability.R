`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

filesystem_safe <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  gsub("(^_+|_+$)", "", x)
}

normalize_filtered_dna_export <- function(payload, export_path = "<in-memory>") {
  if (!is.list(payload) || is.null(payload$raw) || !is.list(payload$raw)) {
    stop(sprintf("Expected '%s' to contain a list with a named '$raw' sample list.", export_path))
  }

  sample_names <- names(payload$raw)
  if (is.null(sample_names) || any(!nzchar(sample_names)) || anyDuplicated(sample_names)) {
    stop(sprintf("Invalid or duplicated sample names in '$raw' for '%s'.", export_path))
  }

  transformed <- payload$transformed
  if (is.null(transformed)) {
    transformed <- payload$raw
  }
  if (!is.list(transformed) || !identical(sample_names, names(transformed))) {
    stop(sprintf("Raw and transformed sample names do not align in '%s'.", export_path))
  }

  list(raw = payload$raw, transformed = transformed)
}

read_filtered_dna_export <- function(export_path) {
  export_path <- normalizePath(export_path, winslash = "/", mustWork = TRUE)
  normalize_filtered_dna_export(readRDS(export_path), export_path)
}

find_filtered_dna_exports <- function(processed_root = "processed_data") {
  files <- list.files(
    processed_root,
    pattern = "^filtered_dna_area_vectors\\.rds$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )
  sort(normalizePath(files, winslash = "/", mustWork = TRUE))
}

dataset_id_from_export <- function(export_path) {
  basename(dirname(normalizePath(export_path, winslash = "/", mustWork = TRUE)))
}

clean_positive_values <- function(values, min_n = 20L) {
  values <- as.numeric(values)
  values <- values[is.finite(values) & values > 0]

  if (length(values) < min_n) {
    stop(sprintf("Need at least %d finite positive values; found %d.", min_n, length(values)))
  }
  if (length(unique(values)) < 3L) {
    stop("Need at least three distinct positive values for KDE peak detection.")
  }

  values
}

compute_kde_peak_prominence <- function(y, peak_idx) {
  peak_y <- y[[peak_idx]]
  left_higher <- which(y[seq_len(peak_idx)] > peak_y)
  left_boundary <- if (length(left_higher)) max(left_higher) else 1L
  right_higher_rel <- which(y[peak_idx:length(y)] > peak_y)
  right_boundary <- if (length(right_higher_rel)) peak_idx + min(right_higher_rel) - 1L else length(y)

  left_base <- min(y[left_boundary:peak_idx], na.rm = TRUE)
  right_base <- min(y[peak_idx:right_boundary], na.rm = TRUE)
  peak_y - max(left_base, right_base)
}

detect_kde_peaks_for_bandwidth <- function(
    values,
    bandwidth_multiplier = 1,
    bandwidth_method = c("SJ", "nrd0"),
    base_bandwidth = NULL,
    transform = c("log10", "raw"),
    density_n = 2048L,
    min_relative_height = 0.01,
    min_relative_prominence = 0.008,
    min_n = 20L) {
  bandwidth_method <- match.arg(bandwidth_method)
  transform <- match.arg(transform)
  raw_values <- clean_positive_values(values, min_n = min_n)
  fit_values <- if (identical(transform, "log10")) log10(raw_values) else raw_values

  base_bw <- base_bandwidth
  if (is.null(base_bw)) {
    base_bw <- if (identical(bandwidth_method, "nrd0")) {
      stats::bw.nrd0(fit_values)
    } else {
      stats::bw.SJ(fit_values)
    }
  }
  if (!is.finite(base_bw) || base_bw <= 0) {
    base_bw <- stats::bw.SJ(fit_values)
  }
  if (!is.finite(base_bw) || base_bw <= 0) {
    stop("Could not estimate a finite KDE bandwidth.")
  }

  dens <- stats::density(
    fit_values,
    bw = base_bw * bandwidth_multiplier,
    n = density_n,
    from = min(fit_values),
    to = max(fit_values),
    na.rm = TRUE
  )

  x <- dens$x
  y <- dens$y
  local_max_idx <- which(y >= c(-Inf, head(y, -1L)) & y >= c(tail(y, -1L), -Inf))
  if (!length(local_max_idx)) {
    local_max_idx <- which.max(y)
  }

  max_y <- max(y, na.rm = TRUE)
  peaks <- lapply(local_max_idx, function(peak_idx) {
    prominence <- compute_kde_peak_prominence(y, peak_idx)
    peak_fit <- x[[peak_idx]]
    data.frame(
      bandwidth_multiplier = bandwidth_multiplier,
      bandwidth_method = bandwidth_method,
      base_bandwidth = base_bw,
      bandwidth = base_bw * bandwidth_multiplier,
      peak_fit = peak_fit,
      peak_raw = if (identical(transform, "log10")) 10^peak_fit else peak_fit,
      peak_height = y[[peak_idx]],
      relative_height = y[[peak_idx]] / max_y,
      prominence = prominence,
      relative_prominence = prominence / max_y,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, peaks)
  out <- out[
    is.finite(out$relative_height) &
      out$relative_height >= min_relative_height &
      is.finite(out$relative_prominence) &
      out$relative_prominence >= min_relative_prominence,
    ,
    drop = FALSE
  ]
  rownames(out) <- NULL
  out
}

cluster_kde_peak_calls <- function(
    peak_calls,
    cluster_tolerance = 0.035,
    min_bandwidth_support = 4L) {
  if (is.null(peak_calls) || !nrow(peak_calls)) {
    return(data.frame(
      cluster_id = integer(),
      stable = logical(),
      bandwidth_support = integer(),
      n_peak_calls = integer(),
      peak_fit_median = numeric(),
      peak_raw_median = numeric(),
      peak_raw_min = numeric(),
      peak_raw_max = numeric(),
      mean_relative_prominence = numeric(),
      stable_rank = integer()
    ))
  }

  peak_calls <- peak_calls[order(peak_calls$peak_fit), , drop = FALSE]
  clusters <- list()

  for (i in seq_len(nrow(peak_calls))) {
    assigned <- FALSE
    for (j in seq_along(clusters)) {
      center <- stats::median(peak_calls$peak_fit[clusters[[j]]], na.rm = TRUE)
      if (abs(peak_calls$peak_fit[[i]] - center) <= cluster_tolerance) {
        clusters[[j]] <- c(clusters[[j]], i)
        assigned <- TRUE
        break
      }
    }
    if (!assigned) {
      clusters[[length(clusters) + 1L]] <- i
    }
  }

  cluster_tbl <- do.call(rbind, lapply(seq_along(clusters), function(cluster_id) {
    idx <- clusters[[cluster_id]]
    cluster_peaks <- peak_calls[idx, , drop = FALSE]
    support <- length(unique(cluster_peaks$bandwidth_multiplier))

    data.frame(
      cluster_id = cluster_id,
      stable = support >= min_bandwidth_support,
      bandwidth_support = support,
      n_peak_calls = nrow(cluster_peaks),
      peak_fit_median = stats::median(cluster_peaks$peak_fit, na.rm = TRUE),
      peak_raw_median = stats::median(cluster_peaks$peak_raw, na.rm = TRUE),
      peak_raw_min = min(cluster_peaks$peak_raw, na.rm = TRUE),
      peak_raw_max = max(cluster_peaks$peak_raw, na.rm = TRUE),
      mean_relative_prominence = mean(cluster_peaks$relative_prominence, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))

  cluster_tbl <- cluster_tbl[order(cluster_tbl$peak_raw_median), , drop = FALSE]
  cluster_tbl$stable_rank <- ifelse(cluster_tbl$stable, cumsum(cluster_tbl$stable), NA_integer_)
  rownames(cluster_tbl) <- NULL
  cluster_tbl
}

detect_kde_stable_peaks <- function(
    values,
    bandwidth_multipliers = c(0.70, 0.85, 1.00, 1.15, 1.30),
    bandwidth_method = c("SJ", "nrd0"),
    base_bandwidth = NULL,
    transform = c("log10", "raw"),
    density_n = 2048L,
    min_relative_height = 0.01,
    min_relative_prominence = 0.008,
    cluster_tolerance = 0.035,
    min_bandwidth_support = 4L,
    min_n = 20L) {
  bandwidth_method <- match.arg(bandwidth_method)
  transform <- match.arg(transform)

  peak_calls <- do.call(rbind, lapply(
    bandwidth_multipliers,
    function(multiplier) {
      detect_kde_peaks_for_bandwidth(
        values = values,
        bandwidth_multiplier = multiplier,
        bandwidth_method = bandwidth_method,
        base_bandwidth = base_bandwidth,
        transform = transform,
        density_n = density_n,
        min_relative_height = min_relative_height,
        min_relative_prominence = min_relative_prominence,
        min_n = min_n
      )
    }
  ))

  if (is.null(peak_calls)) {
    peak_calls <- data.frame()
  }
  clusters <- cluster_kde_peak_calls(
    peak_calls,
    cluster_tolerance = cluster_tolerance,
    min_bandwidth_support = min_bandwidth_support
  )

  list(
    peak_calls = peak_calls,
    stable_peaks = clusters,
    settings = list(
      bandwidth_multipliers = bandwidth_multipliers,
      bandwidth_method = bandwidth_method,
      base_bandwidth = base_bandwidth,
      transform = transform,
      density_n = density_n,
      min_relative_height = min_relative_height,
      min_relative_prominence = min_relative_prominence,
      cluster_tolerance = cluster_tolerance,
      min_bandwidth_support = min_bandwidth_support,
      min_n = min_n
    )
  )
}

estimate_kde_base_bandwidth <- function(
    values,
    bandwidth_method = c("SJ", "nrd0"),
    transform = c("log10", "raw"),
    min_n = 20L) {
  bandwidth_method <- match.arg(bandwidth_method)
  transform <- match.arg(transform)
  raw_values <- clean_positive_values(values, min_n = min_n)
  fit_values <- if (identical(transform, "log10")) log10(raw_values) else raw_values

  bw <- if (identical(bandwidth_method, "nrd0")) {
    stats::bw.nrd0(fit_values)
  } else {
    stats::bw.SJ(fit_values)
  }
  if (!is.finite(bw) || bw <= 0) {
    stop("Could not estimate a finite positive KDE bandwidth.")
  }
  bw
}

estimate_global_kde_base_bandwidth <- function(
    value_list,
    bandwidth_method = c("SJ", "nrd0"),
    transform = c("log10", "raw"),
    mode = c("none", "median-sample", "pooled"),
    min_n = 20L) {
  bandwidth_method <- match.arg(bandwidth_method)
  transform <- match.arg(transform)
  mode <- match.arg(mode)

  if (identical(mode, "none")) {
    return(NULL)
  }

  if (identical(mode, "median-sample")) {
    sample_bandwidths <- vapply(
      value_list,
      estimate_kde_base_bandwidth,
      numeric(1),
      bandwidth_method = bandwidth_method,
      transform = transform,
      min_n = min_n
    )
    return(stats::median(sample_bandwidths, na.rm = TRUE))
  }

  pooled_values <- unlist(lapply(value_list, clean_positive_values, min_n = min_n), use.names = FALSE)
  estimate_kde_base_bandwidth(
    values = pooled_values,
    bandwidth_method = bandwidth_method,
    transform = transform,
    min_n = min_n
  )
}

plot_kde_stable_peak_histogram <- function(
    values,
    stable_peaks,
    output_path,
    sample_name = NA_character_,
    dataset_id = NA_character_,
    bins = 180L,
    upper_quantile = 0.995,
    width = 1400L,
    height = 900L,
    res = 150L) {
  plot_values <- clean_positive_values(values)
  upper <- as.numeric(stats::quantile(plot_values, upper_quantile, na.rm = TRUE))
  plot_values <- plot_values[plot_values <= upper]

  dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(output_path, width = width, height = height, res = res)
  on.exit(grDevices::dev.off(), add = TRUE)

  title <- sample_name
  if (!is.na(dataset_id) && nzchar(dataset_id)) {
    title <- sprintf("%s / %s", dataset_id, sample_name)
  }

  graphics::hist(
    plot_values,
    breaks = bins,
    col = "grey88",
    border = "grey55",
    main = title,
    xlab = "Raw DNA area",
    ylab = "Filtered singlet event count"
  )

  stable_peaks <- stable_peaks[
    stable_peaks$stable & stable_peaks$peak_raw_median <= upper,
    ,
    drop = FALSE
  ]
  if (nrow(stable_peaks)) {
    line_cols <- grDevices::hcl.colors(nrow(stable_peaks), palette = "Dark 3")
    graphics::abline(v = stable_peaks$peak_raw_median, col = line_cols, lwd = 2)
    graphics::legend(
      "topright",
      legend = sprintf("Stable peak %d: %.0f", stable_peaks$stable_rank, stable_peaks$peak_raw_median),
      col = line_cols,
      lwd = 2,
      bty = "n",
      cex = 0.85
    )
  }

  invisible(output_path)
}

write_kde_peak_outputs <- function(
    result,
    output_dir,
    values,
    sample_name = NA_character_,
    dataset_id = NA_character_,
    include_plot = TRUE,
    plot_bins = 180L) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  peak_calls <- result$peak_calls
  stable_peaks <- result$stable_peaks
  if (nrow(peak_calls)) {
    peak_calls$dataset_id <- dataset_id
    peak_calls$sample_name <- sample_name
    peak_calls <- peak_calls[, c("dataset_id", "sample_name", setdiff(names(peak_calls), c("dataset_id", "sample_name")))]
  }
  if (nrow(stable_peaks)) {
    stable_peaks$dataset_id <- dataset_id
    stable_peaks$sample_name <- sample_name
    stable_peaks <- stable_peaks[, c("dataset_id", "sample_name", setdiff(names(stable_peaks), c("dataset_id", "sample_name")))]
  }

  utils::write.csv(peak_calls, file.path(output_dir, "kde_peak_calls_by_bandwidth.csv"), row.names = FALSE)
  utils::write.csv(stable_peaks, file.path(output_dir, "kde_stable_peaks.csv"), row.names = FALSE)
  saveRDS(result, file.path(output_dir, "kde_peak_stability.rds"))

  plot_path <- NA_character_
  if (isTRUE(include_plot)) {
    plot_path <- file.path(output_dir, "raw_hist_stable_kde_peaks.png")
    plot_kde_stable_peak_histogram(
      values = values,
      stable_peaks = result$stable_peaks,
      output_path = plot_path,
      sample_name = sample_name,
      dataset_id = dataset_id,
      bins = plot_bins
    )
  }

  invisible(list(
    peak_calls = file.path(output_dir, "kde_peak_calls_by_bandwidth.csv"),
    stable_peaks = file.path(output_dir, "kde_stable_peaks.csv"),
    fit = file.path(output_dir, "kde_peak_stability.rds"),
    plot = plot_path
  ))
}
