args <- commandArgs(trailingOnly = TRUE)

project_root <- normalizePath(
  if (length(args) >= 1L) args[[1]] else ".",
  winslash = "/",
  mustWork = TRUE
)

example_dir <- file.path(project_root, "dev", "opencyto_minimal_example")
gated_root <- file.path(example_dir, "output")
figure_dir <- file.path(example_dir, "output", "peak_detection_figures")

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggplot2)
})

source(file.path(project_root, "R", "kde_peak_stability.R"))

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

discover_gated_dirs <- function(gated_root) {
  candidates <- list.dirs(gated_root, recursive = TRUE, full.names = TRUE)
  candidates[basename(candidates) == "gated_flow"]
}

resolve_population <- function(gs, population) {
  pop_paths <- flowWorkspace::gs_get_pop_paths(gs, path = "auto")
  if (population %in% pop_paths) {
    return(population)
  }

  matches <- pop_paths[basename(pop_paths) == population]
  if (length(matches) == 1L) {
    return(matches[[1]])
  }

  NA_character_
}

required_population <- function(gs, population) {
  resolved <- resolve_population(gs, population)
  if (is.na(resolved)) {
    stop("Could not find required population '", population, "'.")
  }
  resolved
}

optional_population <- function(gs, population) {
  resolved <- resolve_population(gs, population)
  if (is.na(resolved)) {
    return(NULL)
  }
  resolved
}

raw_dna_values <- function(gh, population, channel = "DNA_AREA") {
  events <- flowWorkspace::gh_pop_get_data(gh, population)
  if (!channel %in% colnames(events)) {
    stop("Missing channel '", channel, "' in population '", population, "'.")
  }

  values <- as.numeric(flowCore::exprs(events)[, channel])
  values[is.finite(values)]
}

compute_density_prominence <- function(y, peak_idx) {
  peak_y <- y[[peak_idx]]
  left_higher <- which(y[seq_len(peak_idx)] > peak_y)
  left_boundary <- if (length(left_higher)) max(left_higher) else 1L
  right_higher_rel <- which(y[peak_idx:length(y)] > peak_y)
  right_boundary <- if (length(right_higher_rel)) peak_idx + min(right_higher_rel) - 1L else length(y)

  left_base <- min(y[left_boundary:peak_idx], na.rm = TRUE)
  right_base <- min(y[peak_idx:right_boundary], na.rm = TRUE)
  peak_y - max(left_base, right_base)
}

find_density_peak_candidates <- function(values, density_n = 2048L, min_relative_height = 0.01) {
  values <- values[is.finite(values) & values > 0]
  if (length(values) < 50L || length(unique(values)) < 3L) {
    return(data.frame())
  }

  dens <- stats::density(values, n = density_n, from = min(values), to = max(values), na.rm = TRUE)
  x <- dens$x
  y <- dens$y
  local_max_idx <- which(y >= c(-Inf, head(y, -1L)) & y >= c(tail(y, -1L), -Inf))
  if (!length(local_max_idx)) {
    local_max_idx <- which.max(y)
  }

  max_y <- max(y, na.rm = TRUE)
  peaks <- do.call(rbind, lapply(local_max_idx, function(peak_idx) {
    prominence <- compute_density_prominence(y, peak_idx)
    data.frame(
      peak_x = x[[peak_idx]],
      peak_height = y[[peak_idx]],
      relative_height = y[[peak_idx]] / max_y,
      prominence = prominence,
      relative_prominence = prominence / max_y,
      stringsAsFactors = FALSE
    )
  }))

  peaks <- peaks[
    is.finite(peaks$relative_height) &
      peaks$relative_height >= min_relative_height,
    ,
    drop = FALSE
  ]
  peaks[order(peaks$prominence, peaks$peak_height, decreasing = TRUE), , drop = FALSE]
}

select_ratio_separated_peaks <- function(candidates, n_keep = 3L, min_ratio = 1.5) {
  if (!nrow(candidates)) {
    return(candidates)
  }

  keep_idx <- integer()
  for (i in seq_len(nrow(candidates))) {
    loc <- candidates$peak_x[[i]]
    separated <- !length(keep_idx) || all(
      pmax(loc, candidates$peak_x[keep_idx]) / pmin(loc, candidates$peak_x[keep_idx]) >= min_ratio
    )
    if (isTRUE(separated)) {
      keep_idx <- c(keep_idx, i)
    }
    if (length(keep_idx) >= n_keep) {
      break
    }
  }

  candidates[keep_idx, , drop = FALSE][order(candidates$peak_x[keep_idx]), , drop = FALSE]
}

detect_gate_split_peaks <- function(gh, gs) {
  cen_pop <- optional_population(gs, "cen_dna_support")
  tumor_pop <- optional_population(gs, "tumor_dna_support")
  peak_parts <- list()

  if (!is.null(cen_pop)) {
    cen_candidates <- find_density_peak_candidates(raw_dna_values(gh, cen_pop))
    if (nrow(cen_candidates)) {
      cen_peak <- cen_candidates[1L, , drop = FALSE]
      cen_peak$method <- "1. gate split KDE"
      cen_peak$peak_type <- "CEN"
      peak_parts[[length(peak_parts) + 1L]] <- cen_peak
    }
  }

  if (!is.null(tumor_pop)) {
    tumor_candidates <- find_density_peak_candidates(raw_dna_values(gh, tumor_pop))
    tumor_peaks <- select_ratio_separated_peaks(tumor_candidates, n_keep = 3L, min_ratio = 1.5)
    if (nrow(tumor_peaks)) {
      tumor_peaks$method <- "1. gate split KDE"
      tumor_peaks$peak_type <- paste0("Tumor ", seq_len(nrow(tumor_peaks)))
      peak_parts[[length(peak_parts) + 1L]] <- tumor_peaks
    }
  }

  if (!length(peak_parts)) {
    return(data.frame())
  }

  do.call(rbind, peak_parts)
}

detect_stability_peaks <- function(values) {
  result <- tryCatch(
    detect_kde_stable_peaks(
      values = values,
      bandwidth_multipliers = c(0.70, 0.85, 1.00, 1.15, 1.30),
      bandwidth_method = "SJ",
      transform = "log10",
      density_n = 2048L,
      min_relative_height = 0.01,
      min_relative_prominence = 0.008,
      cluster_tolerance = 0.035,
      min_bandwidth_support = 4L,
      min_n = 50L
    ),
    error = function(e) NULL
  )
  if (is.null(result) || !nrow(result$stable_peaks)) {
    return(data.frame())
  }

  stable_peaks <- result$stable_peaks[result$stable_peaks$stable, , drop = FALSE]
  if (!nrow(stable_peaks)) {
    return(data.frame())
  }

  data.frame(
    peak_x = stable_peaks$peak_raw_median,
    peak_height = NA_real_,
    relative_height = NA_real_,
    prominence = NA_real_,
    relative_prominence = stable_peaks$mean_relative_prominence,
    method = "2. log10 KDE stability",
    peak_type = paste0("Stable ", stable_peaks$stable_rank),
    stringsAsFactors = FALSE
  )
}

refine_consensus_peak <- function(values, x1, x2, pad_log10 = 0.025, density_n = 4096L) {
  values <- values[is.finite(values) & values > 0]
  if (length(values) < 50L || !is.finite(x1) || !is.finite(x2) || x1 <= 0 || x2 <= 0) {
    return(sqrt(x1 * x2))
  }

  fit_values <- log10(values)
  lower <- min(log10(c(x1, x2))) - pad_log10
  upper <- max(log10(c(x1, x2))) + pad_log10
  if (!is.finite(lower) || !is.finite(upper) || lower >= upper) {
    return(sqrt(x1 * x2))
  }

  dens <- stats::density(
    fit_values,
    bw = stats::bw.SJ(fit_values),
    n = density_n,
    from = min(fit_values),
    to = max(fit_values),
    na.rm = TRUE
  )
  in_window <- which(dens$x >= lower & dens$x <= upper)
  if (!length(in_window)) {
    return(sqrt(x1 * x2))
  }

  local_max_idx <- in_window[
    dens$y[in_window] >= c(-Inf, dens$y[head(in_window, -1L)]) &
      dens$y[in_window] >= c(dens$y[tail(in_window, -1L)], -Inf)
  ]
  if (!length(local_max_idx)) {
    local_max_idx <- in_window[[which.max(dens$y[in_window])]]
  } else {
    local_max_idx <- local_max_idx[[which.max(dens$y[local_max_idx])]]
  }

  10^dens$x[[local_max_idx]]
}

build_agreed_peak_set <- function(plot_values, peak_tbl, max_log10_distance = 0.08) {
  gate_peaks <- peak_tbl[peak_tbl$method == "1. gate split KDE" & is.finite(peak_tbl$peak_x), , drop = FALSE]
  stable_peaks <- peak_tbl[peak_tbl$method == "2. log10 KDE stability" & is.finite(peak_tbl$peak_x), , drop = FALSE]
  if (!nrow(gate_peaks) || !nrow(stable_peaks)) {
    return(data.frame())
  }

  distance_mat <- abs(outer(log10(gate_peaks$peak_x), log10(stable_peaks$peak_x), "-"))
  matches <- list()
  for (i in seq_len(nrow(gate_peaks))) {
    j <- which.min(distance_mat[i, ])
    stable_nearest_gate <- which.min(distance_mat[, j])
    if (identical(stable_nearest_gate, i) && distance_mat[i, j] <= max_log10_distance) {
      gate_x <- gate_peaks$peak_x[[i]]
      stable_x <- stable_peaks$peak_x[[j]]
      consensus_x <- refine_consensus_peak(plot_values, gate_x, stable_x)
      matches[[length(matches) + 1L]] <- data.frame(
        agreed_peak_id = length(matches) + 1L,
        peak_x = consensus_x,
        peak_type = gate_peaks$peak_type[[i]],
        gate_split_peak_x = gate_x,
        stability_peak_x = stable_x,
        log10_distance = distance_mat[i, j],
        ratio_distance = max(gate_x, stable_x) / min(gate_x, stable_x),
        stringsAsFactors = FALSE
      )
    }
  }

  if (!length(matches)) {
    return(data.frame())
  }

  out <- do.call(rbind, matches)
  out <- out[order(out$peak_x), , drop = FALSE]
  out$agreed_peak_id <- seq_len(nrow(out))
  out
}

plot_method_comparison <- function(plot_values, peak_tbl, agreed_tbl, dataset, sample_name, output_file) {
  plot_values <- plot_values[is.finite(plot_values)]
  methods <- c("1. gate split KDE", "2. log10 KDE stability")
  hist_df <- do.call(rbind, lapply(methods, function(method) {
    data.frame(value = plot_values, method = method, stringsAsFactors = FALSE)
  }))

  peak_tbl <- peak_tbl[is.finite(peak_tbl$peak_x), , drop = FALSE]
  hist_max <- max(graphics::hist(plot_values, breaks = 180, plot = FALSE)$counts, na.rm = TRUE)
  peak_tbl$label_y <- hist_max * 1.08
  peak_tbl$label <- sprintf("%s\n%.0f", peak_tbl$peak_type, peak_tbl$peak_x)
  peak_tbl$method <- factor(peak_tbl$method, levels = methods)
  hist_df$method <- factor(hist_df$method, levels = methods)
  agreed_plot_tbl <- NULL
  if (nrow(agreed_tbl)) {
    agreed_plot_tbl <- do.call(rbind, lapply(methods, function(method) {
      data.frame(
        agreed_tbl,
        method = method,
        label_y = hist_max * 1.19,
        label = sprintf("Agreed %d\n%.0f", agreed_tbl$agreed_peak_id, agreed_tbl$peak_x),
        stringsAsFactors = FALSE
      )
    }))
    agreed_plot_tbl$method <- factor(agreed_plot_tbl$method, levels = methods)
  }

  p <- ggplot(hist_df, aes(x = value)) +
    geom_histogram(bins = 180, fill = "grey82", color = "grey55", linewidth = 0.12) +
    facet_wrap(~method, ncol = 1, scales = "free_y") +
    scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.20))) +
    labs(
      title = paste(dataset, sample_name, sep = " / "),
      subtitle = paste0("cells_keep raw DNA_AREA | n = ", length(plot_values)),
      x = "Raw DNA-A",
      y = "Event count"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 9),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold")
    )

  if (nrow(peak_tbl)) {
    p <- p +
      geom_vline(
        data = peak_tbl,
        aes(xintercept = peak_x, color = peak_type),
        linewidth = 0.55,
        show.legend = FALSE
      ) +
      geom_label(
        data = peak_tbl,
        aes(x = peak_x, y = label_y, label = label, color = peak_type),
        inherit.aes = FALSE,
        vjust = 1.1,
        angle = 90,
        size = 2.5,
        label.size = 0.15,
        fill = "white",
        show.legend = FALSE
      )
  }

  if (!is.null(agreed_plot_tbl) && nrow(agreed_plot_tbl)) {
    p <- p +
      geom_vline(
        data = agreed_plot_tbl,
        aes(xintercept = peak_x),
        color = "black",
        linetype = "dashed",
        linewidth = 0.65
      ) +
      geom_label(
        data = agreed_plot_tbl,
        aes(x = peak_x, y = label_y, label = label),
        inherit.aes = FALSE,
        vjust = 1.1,
        angle = 90,
        size = 2.5,
        label.size = 0.2,
        fill = "white",
        color = "black"
      )
  }

  ggplot2::ggsave(output_file, p, width = 12, height = 8, dpi = 150, bg = "white")
}

dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(figure_dir, pattern = "\\.png$", full.names = TRUE), force = TRUE)
peak_call_file <- file.path(figure_dir, "peak_calls_by_method.csv")
agreed_peak_file <- file.path(figure_dir, "agreed_peaks.csv")

gated_dirs <- discover_gated_dirs(gated_root)
if (length(gated_dirs) == 0L) {
  stop("No gated_flow directories found under ", gated_root)
}

plot_count <- 0L
all_peak_calls <- list()
all_agreed_peaks <- list()
for (gated_dir in gated_dirs) {
  dataset <- basename(dirname(gated_dir))
  cat("Loading ", dataset, "\n", sep = "")
  gs <- flowWorkspace::load_gs(gated_dir)
  keep_pop <- required_population(gs, "cells_keep")

  for (sample_name in flowCore::sampleNames(gs)) {
    gh <- gs[[sample_name]]
    keep_values <- raw_dna_values(gh, keep_pop)

    gate_split_peaks <- detect_gate_split_peaks(gh, gs)
    stability_peaks <- detect_stability_peaks(keep_values)
    peak_tbl <- rbind(gate_split_peaks, stability_peaks)
    agreed_tbl <- build_agreed_peak_set(keep_values, peak_tbl)

    if (nrow(peak_tbl)) {
      peak_tbl$dataset <- dataset
      peak_tbl$sample_name <- sample_name
      peak_tbl$n_events_cells_keep <- length(keep_values)
      all_peak_calls[[length(all_peak_calls) + 1L]] <- peak_tbl[
        ,
        c("dataset", "sample_name", "n_events_cells_keep", setdiff(names(peak_tbl), c("dataset", "sample_name", "n_events_cells_keep"))),
        drop = FALSE
      ]
    }
    if (nrow(agreed_tbl)) {
      agreed_tbl$dataset <- dataset
      agreed_tbl$sample_name <- sample_name
      agreed_tbl$n_events_cells_keep <- length(keep_values)
      all_agreed_peaks[[length(all_agreed_peaks) + 1L]] <- agreed_tbl[
        ,
        c("dataset", "sample_name", "n_events_cells_keep", setdiff(names(agreed_tbl), c("dataset", "sample_name", "n_events_cells_keep"))),
        drop = FALSE
      ]
    }

    output_file <- file.path(
      figure_dir,
      paste0(
        sanitize_filename(dataset),
        "__",
        sanitize_filename(sample_name),
        ".png"
      )
    )

    plot_method_comparison(keep_values, peak_tbl, agreed_tbl, dataset, sample_name, output_file)
    plot_count <- plot_count + 1L
  }
}

peak_call_tbl <- if (length(all_peak_calls)) do.call(rbind, all_peak_calls) else data.frame()
agreed_peak_tbl <- if (length(all_agreed_peaks)) do.call(rbind, all_agreed_peaks) else data.frame()
utils::write.csv(peak_call_tbl, peak_call_file, row.names = FALSE)
utils::write.csv(agreed_peak_tbl, agreed_peak_file, row.names = FALSE)

cat("Wrote ", plot_count, " peak-detection histogram(s) to ", figure_dir, "\n", sep = "")
cat("Wrote method peak calls to ", peak_call_file, "\n", sep = "")
cat("Wrote agreed peak calls to ", agreed_peak_file, "\n", sep = "")
