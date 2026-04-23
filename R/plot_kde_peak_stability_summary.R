args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

parse_cli_args <- function(args) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected positional argument '%s'. Use --key value arguments.", key))
    }
    key <- sub("^--", "", key)
    if (i == length(args)) {
      stop(sprintf("Missing value for --%s.", key))
    }
    out[[key]] <- args[[i + 1L]]
    i <- i + 2L
  }
  out
}

require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required. Install it with install.packages('ggplot2').")
  }
}

filesystem_safe <- function(x) {
  x <- gsub("[^A-Za-z0-9_.-]+", "_", x)
  gsub("(^_+|_+$)", "", x)
}

build_pair_ratios <- function(peaks) {
  split_peaks <- split(peaks, peaks$sample_key)
  ratio_parts <- lapply(split_peaks, function(sample_peaks) {
    sample_peaks <- sample_peaks[order(sample_peaks$peak_raw_median), , drop = FALSE]
    if (nrow(sample_peaks) < 2L) {
      return(NULL)
    }

    pair_idx <- utils::combn(seq_len(nrow(sample_peaks)), 2L)
    data.frame(
      dataset_id = sample_peaks$dataset_id[[1L]],
      sample_name = sample_peaks$sample_name[[1L]],
      sample_key = sample_peaks$sample_key[[1L]],
      lower_peak_id = sample_peaks$peak_id[pair_idx[1L, ]],
      upper_peak_id = sample_peaks$peak_id[pair_idx[2L, ]],
      lower_stable_rank = sample_peaks$stable_rank[pair_idx[1L, ]],
      upper_stable_rank = sample_peaks$stable_rank[pair_idx[2L, ]],
      lower_peak_raw = sample_peaks$peak_raw_median[pair_idx[1L, ]],
      upper_peak_raw = sample_peaks$peak_raw_median[pair_idx[2L, ]],
      lower_mean_relative_prominence = sample_peaks$mean_relative_prominence[pair_idx[1L, ]],
      upper_mean_relative_prominence = sample_peaks$mean_relative_prominence[pair_idx[2L, ]],
      ratio = sample_peaks$peak_raw_median[pair_idx[2L, ]] / sample_peaks$peak_raw_median[pair_idx[1L, ]],
      adjacent_pair = pair_idx[2L, ] == pair_idx[1L, ] + 1L,
      stringsAsFactors = FALSE
    )
  })

  out <- do.call(rbind, ratio_parts[!vapply(ratio_parts, is.null, logical(1L))])
  rownames(out) <- NULL
  out$log2_ratio <- log2(out$ratio)
  out
}

estimate_g0g1_ratio_targets <- function(ratios, min_ratio = 1.55, max_ratio = 2.15, bin_width = 0.04) {
  candidates <- ratios[
    is.finite(ratios$ratio) &
      ratios$ratio >= min_ratio &
      ratios$ratio <= max_ratio,
    ,
    drop = FALSE
  ]
  if (!nrow(candidates)) {
    return(data.frame(
      dataset_id = character(),
      log2_ratio_bin = numeric(),
      target_ratio = numeric(),
      n_pairs = integer(),
      n_samples = integer(),
      lower_ratio = numeric(),
      upper_ratio = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  candidates$log2_ratio_bin <- round(candidates$log2_ratio / bin_width) * bin_width
  bins <- unique(candidates[, c("dataset_id", "log2_ratio_bin")])
  target_parts <- lapply(seq_len(nrow(bins)), function(i) {
    keep <- candidates$dataset_id == bins$dataset_id[[i]] &
      candidates$log2_ratio_bin == bins$log2_ratio_bin[[i]]
    bin_tbl <- candidates[keep, , drop = FALSE]
    data.frame(
      dataset_id = bins$dataset_id[[i]],
      log2_ratio_bin = bins$log2_ratio_bin[[i]],
      target_ratio = stats::median(bin_tbl$ratio, na.rm = TRUE),
      n_pairs = nrow(bin_tbl),
      n_samples = length(unique(bin_tbl$sample_key)),
      lower_ratio = min(bin_tbl$ratio, na.rm = TRUE),
      upper_ratio = max(bin_tbl$ratio, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })

  target_tbl <- do.call(rbind, target_parts)
  target_tbl <- target_tbl[order(
    target_tbl$dataset_id,
    -target_tbl$n_samples,
    -target_tbl$n_pairs,
    abs(target_tbl$log2_ratio_bin - 1)
  ), , drop = FALSE]
  target_tbl <- target_tbl[!duplicated(target_tbl$dataset_id), , drop = FALSE]
  target_tbl <- target_tbl[order(target_tbl$dataset_id), , drop = FALSE]
  rownames(target_tbl) <- NULL
  target_tbl
}

label_g0g1_peak_chains <- function(peaks, ratios, targets, tolerance_log2 = 0.12) {
  peaks$is_g0g1 <- FALSE
  peaks$g0g1_chain_id <- NA_character_
  peaks$g0g1_chain_size <- NA_integer_
  peaks$g0g1_role <- "Other stable peak"

  if (!nrow(targets)) {
    return(list(peaks = peaks, pairs = ratios[FALSE, , drop = FALSE]))
  }

  target_map <- stats::setNames(targets$target_ratio, targets$dataset_id)
  ratios$g0g1_target_ratio <- unname(target_map[ratios$dataset_id])
  g0g1_pairs <- ratios[
    is.finite(ratios$g0g1_target_ratio) &
      abs(log2(ratios$ratio / ratios$g0g1_target_ratio)) <= tolerance_log2,
    ,
    drop = FALSE
  ]
  g0g1_pairs$g0g1_ratio_error_log2 <- log2(g0g1_pairs$ratio / g0g1_pairs$g0g1_target_ratio)

  if (!nrow(g0g1_pairs)) {
    return(list(peaks = peaks, pairs = g0g1_pairs))
  }

  chain_parts <- list()
  split_pairs <- split(g0g1_pairs, g0g1_pairs$sample_key)
  chain_idx <- 0L
  for (sample_key in names(split_pairs)) {
    sample_pairs <- split_pairs[[sample_key]]
    sample_peak_ids <- unique(c(sample_pairs$lower_peak_id, sample_pairs$upper_peak_id))
    adjacency <- lapply(sample_peak_ids, function(id) {
      unique(c(
        sample_pairs$upper_peak_id[sample_pairs$lower_peak_id == id],
        sample_pairs$lower_peak_id[sample_pairs$upper_peak_id == id]
      ))
    })
    names(adjacency) <- as.character(sample_peak_ids)

    unvisited <- as.character(sample_peak_ids)
    while (length(unvisited)) {
      queue <- unvisited[[1L]]
      component <- character()
      while (length(queue)) {
        current <- queue[[1L]]
        queue <- queue[-1L]
        if (current %in% component) {
          next
        }
        component <- c(component, current)
        neighbors <- as.character(adjacency[[current]])
        queue <- unique(c(queue, setdiff(neighbors, component)))
      }
      unvisited <- setdiff(unvisited, component)
      component_ids <- as.integer(component)
      if (length(component_ids) >= 2L) {
        chain_idx <- chain_idx + 1L
        chain_parts[[length(chain_parts) + 1L]] <- data.frame(
          peak_id = component_ids,
          g0g1_chain_id = sprintf("g0g1_%03d", chain_idx),
          g0g1_chain_size = length(component_ids),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  chain_tbl <- do.call(rbind, chain_parts)
  match_idx <- match(peaks$peak_id, chain_tbl$peak_id)
  peaks$is_g0g1 <- !is.na(match_idx)
  peaks$g0g1_chain_id <- chain_tbl$g0g1_chain_id[match_idx]
  peaks$g0g1_chain_size <- chain_tbl$g0g1_chain_size[match_idx]
  peaks$g0g1_role <- ifelse(
    peaks$is_g0g1,
    ifelse(peaks$g0g1_chain_size >= 3L, "G0/G1 triplet", "G0/G1 pair"),
    "Other stable peak"
  )

  list(peaks = peaks, pairs = g0g1_pairs)
}

summarize_ratio_support <- function(ratios, bin_width = 0.05, max_ratio = 32) {
  in_range <- ratios[
    is.finite(ratios$ratio) &
      ratios$ratio >= 1 &
      ratios$ratio <= max_ratio,
    ,
    drop = FALSE
  ]
  in_range$log2_ratio_bin <- round(in_range$log2_ratio / bin_width) * bin_width
  in_range$ratio_bin_center <- 2^in_range$log2_ratio_bin
  in_range$ratio_label <- sprintf("%.2f", in_range$ratio_bin_center)

  keys <- unique(in_range[, c("dataset_id", "log2_ratio_bin", "ratio_bin_center", "ratio_label")])
  support <- lapply(seq_len(nrow(keys)), function(i) {
    keep <- in_range$dataset_id == keys$dataset_id[[i]] &
      in_range$log2_ratio_bin == keys$log2_ratio_bin[[i]]
    bin_tbl <- in_range[keep, , drop = FALSE]
    data.frame(
      dataset_id = keys$dataset_id[[i]],
      log2_ratio_bin = keys$log2_ratio_bin[[i]],
      ratio_bin_center = keys$ratio_bin_center[[i]],
      ratio_label = keys$ratio_label[[i]],
      n_pairs = nrow(bin_tbl),
      n_samples = length(unique(bin_tbl$sample_key)),
      n_adjacent_pairs = sum(bin_tbl$adjacent_pair),
      n_adjacent_samples = length(unique(bin_tbl$sample_key[bin_tbl$adjacent_pair])),
      stringsAsFactors = FALSE
    )
  })

  support <- do.call(rbind, support)
  support <- support[support$ratio_bin_center <= max_ratio, , drop = FALSE]
  support <- support[order(support$dataset_id, support$log2_ratio_bin), , drop = FALSE]
  rownames(support) <- NULL
  support
}

theme_peak <- function(base_size = 10) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid.minor = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold"),
      axis.text.y = ggplot2::element_text(size = 5.5),
      legend.position = "bottom"
    )
}

cli <- parse_cli_args(args)
input_path <- cli[["input"]] %||% file.path("processed_data", "kde_peak_stability", "combined_kde_stable_peaks.csv")
output_dir <- cli[["output-dir"]] %||% file.path("figure", "kde_peak_stability")
ratio_max <- as.numeric(cli[["ratio-max"]] %||% "32")
ratio_bin_width <- as.numeric(cli[["ratio-bin-width"]] %||% "0.05")
g0g1_min_ratio <- as.numeric(cli[["g0g1-min-ratio"]] %||% "1.55")
g0g1_max_ratio <- as.numeric(cli[["g0g1-max-ratio"]] %||% "2.15")
g0g1_bin_width <- as.numeric(cli[["g0g1-bin-width"]] %||% "0.04")
g0g1_tolerance_log2 <- as.numeric(cli[["g0g1-tolerance-log2"]] %||% "0.12")

if (!is.finite(ratio_max) || ratio_max <= 1) {
  stop("--ratio-max must be greater than 1.")
}
if (!is.finite(ratio_bin_width) || ratio_bin_width <= 0) {
  stop("--ratio-bin-width must be positive.")
}
if (!is.finite(g0g1_min_ratio) || !is.finite(g0g1_max_ratio) || g0g1_min_ratio <= 1 || g0g1_min_ratio >= g0g1_max_ratio) {
  stop("--g0g1-min-ratio and --g0g1-max-ratio must define a valid ratio interval above 1.")
}
if (!is.finite(g0g1_bin_width) || g0g1_bin_width <= 0) {
  stop("--g0g1-bin-width must be positive.")
}
if (!is.finite(g0g1_tolerance_log2) || g0g1_tolerance_log2 <= 0) {
  stop("--g0g1-tolerance-log2 must be positive.")
}

require_ggplot2()
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

peaks <- read.csv(input_path, stringsAsFactors = FALSE)
required <- c(
  "dataset_id", "sample_name", "stable", "bandwidth_support", "peak_fit_median",
  "peak_raw_median", "mean_relative_prominence", "stable_rank"
)
missing <- setdiff(required, names(peaks))
if (length(missing)) {
  stop(sprintf("Input is missing required columns: %s", paste(missing, collapse = ", ")))
}

peaks <- peaks[peaks$stable & is.finite(peaks$peak_raw_median) & peaks$peak_raw_median > 0, , drop = FALSE]
peaks$peak_id <- seq_len(nrow(peaks))
peaks$sample_key <- paste(peaks$dataset_id, peaks$sample_name, sep = " / ")
peaks$peak_log10 <- log10(peaks$peak_raw_median)
peaks$prominence_label <- cut(
  peaks$mean_relative_prominence,
  breaks = c(-Inf, 0.025, 0.10, 0.50, Inf),
  labels = c("<0.025", "0.025-0.10", "0.10-0.50", ">=0.50")
)

sample_order_tbl <- aggregate(
  peak_raw_median ~ dataset_id + sample_key,
  data = peaks,
  FUN = min
)
sample_order_tbl <- sample_order_tbl[order(sample_order_tbl$dataset_id, sample_order_tbl$peak_raw_median, sample_order_tbl$sample_key), ]
peaks$sample_key <- factor(peaks$sample_key, levels = rev(sample_order_tbl$sample_key))

location_distribution <- ggplot2::ggplot(peaks, ggplot2::aes(x = peak_log10)) +
  ggplot2::geom_histogram(ggplot2::aes(y = ggplot2::after_stat(density)), bins = 55, fill = "#8fb9aa", color = "white", linewidth = 0.15) +
  ggplot2::geom_density(color = "#253c4f", linewidth = 0.7, adjust = 0.8) +
  ggplot2::geom_rug(ggplot2::aes(color = prominence_label), alpha = 0.35, sides = "b") +
  ggplot2::facet_wrap(ggplot2::vars(dataset_id), ncol = 1L, scales = "free_y") +
  ggplot2::scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
  ggplot2::labs(
    title = "Distribution of KDE-stable peak locations",
    subtitle = "Rugs show individual stable peak calls; x-axis is log10 raw DNA-A peak location.",
    x = "log10(raw DNA-A peak location)",
    y = "Density",
    color = "Mean relative\nprominence"
  ) +
  theme_peak(base_size = 11)

location_raster <- ggplot2::ggplot(
  peaks,
  ggplot2::aes(x = peak_log10, y = sample_key)
) +
  ggplot2::geom_point(
    ggplot2::aes(size = bandwidth_support, color = prominence_label),
    alpha = 0.82
  ) +
  ggplot2::facet_grid(dataset_id ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_color_brewer(palette = "Dark2", na.translate = FALSE) +
  ggplot2::scale_y_discrete(labels = function(x) sub("^.* / ", "", x)) +
  ggplot2::scale_size_continuous(range = c(1.2, 3.2), breaks = sort(unique(peaks$bandwidth_support))) +
  ggplot2::labs(
    title = "KDE-stable peak locations by sample",
    subtitle = "Samples are ordered within dataset by their lowest stable raw DNA-A peak.",
    x = "log10(raw DNA-A peak location)",
    y = NULL,
    color = "Mean relative\nprominence",
    size = "Bandwidth\nsupport"
  ) +
  theme_peak(base_size = 10)

ratios <- build_pair_ratios(peaks)
ratio_support <- summarize_ratio_support(ratios, bin_width = ratio_bin_width, max_ratio = ratio_max)
write.csv(ratios, file.path(output_dir, "kde_stable_peak_pair_ratios.csv"), row.names = FALSE)
write.csv(ratio_support, file.path(output_dir, "kde_stable_peak_ratio_support_by_bin.csv"), row.names = FALSE)

g0g1_targets <- estimate_g0g1_ratio_targets(
  ratios,
  min_ratio = g0g1_min_ratio,
  max_ratio = g0g1_max_ratio,
  bin_width = g0g1_bin_width
)
g0g1_parts <- label_g0g1_peak_chains(
  peaks,
  ratios,
  g0g1_targets,
  tolerance_log2 = g0g1_tolerance_log2
)
peaks <- g0g1_parts$peaks
g0g1_pairs <- g0g1_parts$pairs
write.csv(g0g1_targets, file.path(output_dir, "kde_stable_peak_g0g1_ratio_targets.csv"), row.names = FALSE)
write.csv(g0g1_pairs, file.path(output_dir, "kde_stable_peak_g0g1_pairs.csv"), row.names = FALSE)
write.csv(peaks, file.path(output_dir, "kde_stable_peak_g0g1_assignments.csv"), row.names = FALSE)

g0g1_raster <- ggplot2::ggplot(peaks, ggplot2::aes(x = peak_log10, y = sample_key)) +
  ggplot2::geom_point(
    ggplot2::aes(size = bandwidth_support),
    color = "#b8b8b8",
    alpha = 0.55
  ) +
  ggplot2::geom_point(
    data = peaks[peaks$is_g0g1, , drop = FALSE],
    ggplot2::aes(size = bandwidth_support, shape = g0g1_role),
    color = "#d95f02",
    alpha = 0.95
  ) +
  ggplot2::facet_grid(dataset_id ~ ., scales = "free_y", space = "free_y") +
  ggplot2::scale_y_discrete(labels = function(x) sub("^.* / ", "", x)) +
  ggplot2::scale_size_continuous(range = c(1.2, 3.2), breaks = sort(unique(peaks$bandwidth_support))) +
  ggplot2::scale_shape_manual(values = c("G0/G1 pair" = 16, "G0/G1 triplet" = 17)) +
  ggplot2::labs(
    title = "KDE-stable peak locations with inferred G0/G1 peaks highlighted",
    subtitle = sprintf(
      "G0/G1 candidates are stable peaks connected by dataset-specific ratios from %.2f-%.2f, tolerance %.2f log2 units.",
      g0g1_min_ratio,
      g0g1_max_ratio,
      g0g1_tolerance_log2
    ),
    x = "log10(raw DNA-A peak location)",
    y = NULL,
    size = "Bandwidth\nsupport",
    shape = "G0/G1\nchain"
  ) +
  theme_peak(base_size = 10)

ratio_breaks <- c(1, 1.25, 1.5, 2, 3, 4, 6, 8, 16, 32)
ratio_breaks <- ratio_breaks[ratio_breaks <= ratio_max]
ratio_support_plot <- ggplot2::ggplot(
  ratio_support,
  ggplot2::aes(x = log2_ratio_bin, y = dataset_id, fill = n_samples)
) +
  ggplot2::geom_tile(height = 0.78) +
  ggplot2::geom_text(
    data = ratio_support[ratio_support$n_samples >= max(2L, stats::quantile(ratio_support$n_samples, 0.85, na.rm = TRUE)), , drop = FALSE],
    ggplot2::aes(label = n_samples),
    size = 2.4,
    color = "white"
  ) +
  ggplot2::scale_x_continuous(
    breaks = log2(ratio_breaks),
    labels = ratio_breaks,
    limits = c(0, log2(ratio_max)),
    expand = ggplot2::expansion(mult = c(0, 0))
  ) +
  ggplot2::scale_fill_gradient(low = "#d9e6e2", high = "#1f5b63") +
  ggplot2::labs(
    title = "Within-sample stable peak ratio support",
    subtitle = sprintf("All stable peak pairs per sample; bins are %.2f on the log2-ratio scale. Tile labels mark high-support bins.", ratio_bin_width),
    x = "Upper peak / lower peak ratio",
    y = NULL,
    fill = "Samples\nwith ratio"
  ) +
  theme_peak(base_size = 11) +
  ggplot2::theme(axis.text.y = ggplot2::element_text(size = 9))

adjacent_ratios <- ratios[ratios$adjacent_pair & ratios$ratio <= ratio_max, , drop = FALSE]
adjacent_ratio_plot <- ggplot2::ggplot(adjacent_ratios, ggplot2::aes(x = log2_ratio, fill = dataset_id)) +
  ggplot2::geom_histogram(bins = 70, alpha = 0.72, position = "identity") +
  ggplot2::facet_wrap(ggplot2::vars(dataset_id), ncol = 1L, scales = "free_y") +
  ggplot2::scale_x_continuous(
    breaks = log2(ratio_breaks),
    labels = ratio_breaks
  ) +
  ggplot2::coord_cartesian(xlim = c(0, log2(ratio_max))) +
  ggplot2::scale_fill_brewer(palette = "Set2", guide = "none") +
  ggplot2::labs(
    title = "Adjacent stable peak ratios",
    subtitle = "Uses only neighboring stable peaks after sorting each sample by raw DNA-A location.",
    x = "Adjacent upper peak / lower peak ratio",
    y = "Pair count"
  ) +
  theme_peak(base_size = 11)

ggplot2::ggsave(
  file.path(output_dir, "kde_stable_peak_location_distribution.png"),
  location_distribution,
  width = 9,
  height = 9,
  dpi = 220
)
ggplot2::ggsave(
  file.path(output_dir, "kde_stable_peak_location_by_sample.png"),
  location_raster,
  width = 11,
  height = max(8, length(unique(peaks$sample_key)) * 0.13),
  dpi = 220,
  limitsize = FALSE
)
ggplot2::ggsave(
  file.path(output_dir, "kde_stable_peak_pair_ratio_support.png"),
  ratio_support_plot,
  width = 10,
  height = 4.8,
  dpi = 220
)
ggplot2::ggsave(
  file.path(output_dir, "kde_stable_peak_g0g1_by_sample.png"),
  g0g1_raster,
  width = 11,
  height = max(8, length(unique(peaks$sample_key)) * 0.13),
  dpi = 220,
  limitsize = FALSE
)
ggplot2::ggsave(
  file.path(output_dir, "kde_stable_peak_adjacent_ratio_distribution.png"),
  adjacent_ratio_plot,
  width = 9,
  height = 9,
  dpi = 220
)

cat(sprintf("Read %d stable peaks from %d sample instances.\n", nrow(peaks), length(unique(peaks$sample_key))))
cat(sprintf("Flagged %d G0/G1 candidate peaks in %d pair/triplet chains.\n", sum(peaks$is_g0g1), length(unique(peaks$g0g1_chain_id[peaks$is_g0g1]))))
cat(sprintf("Wrote plots and ratio tables to %s\n", normalizePath(output_dir, winslash = "/", mustWork = TRUE)))
