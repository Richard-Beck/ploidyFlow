args <- commandArgs(trailingOnly = TRUE)

project_root <- normalizePath(
  if (length(args) >= 1L) args[[1]] else ".",
  winslash = "/",
  mustWork = TRUE
)

example_dir <- file.path(project_root, "dev", "opencyto_minimal_example")
output_dir <- file.path(example_dir, "output", "peak_consistency")
qc_plot_file <- file.path(output_dir, "qc_removed_by_dataset.png")
qc_table_file <- file.path(output_dir, "qc_removed_by_dataset.csv")
event_count_plot_file <- file.path(output_dir, "qc_event_counts_by_dataset.png")
event_count_table_file <- file.path(output_dir, "qc_event_counts_by_dataset.csv")
acquisition_time_plot_file <- file.path(output_dir, "qc_acquisition_time_by_dataset.png")
acquisition_time_table_file <- file.path(output_dir, "qc_acquisition_time_by_dataset.csv")
trend_count_plot_file <- file.path(output_dir, "qc_peacoqc_trend_channel_counts_by_dataset.png")
trend_count_table_file <- file.path(output_dir, "qc_peacoqc_trend_channel_counts_by_dataset.csv")
voltage_plot_file <- file.path(output_dir, "qc_voltage_by_dataset.png")
voltage_table_file <- file.path(output_dir, "qc_voltage_by_dataset.csv")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(flowCore)
  library(ggplot2)
})

datasets <- c(
  "Anoxia_FlowCytometry",
  "Hypoxia_SUM159",
  "Polyploidization_EthanolFixation"
)

step_levels <- c("Margin removal", "PeacoQC", "Gating")
event_metric_levels <- c("Raw total", "Gated CEN", "Gated tumor")
trend_metric_levels <- c("Increasing", "Decreasing")
cell_line_colors <- c(
  SUM159_2N = "#0072B2",
  SUM159_4N = "#D55E00",
  HCC1937 = "#009E73",
  Jurkat = "#CC79A7",
  KALS1 = "#E69F00",
  MDAMB231 = "#56B4E9",
  Unknown = "#666666"
)

infer_cell_line <- function(dataset, sample_name) {
  if (dataset == "Anoxia_FlowCytometry") {
    if (grepl("^Sample_CEN(?:\\.fcs|\\.fcs\\.h5)?$", sample_name)) {
      return("Unknown")
    }
    if (grepl("^Sample_SUM159_4N_Control(?:\\.fcs|\\.fcs\\.h5)?$", sample_name)) {
      return("SUM159_4N")
    }
    return("SUM159_2N")
  }

  if (dataset == "Hypoxia_SUM159") {
    if (grepl("_2N_", sample_name, fixed = TRUE)) {
      return("SUM159_2N")
    }
    if (grepl("_4N_", sample_name, fixed = TRUE)) {
      return("SUM159_4N")
    }
    return("Unknown")
  }

  if (dataset == "Polyploidization_EthanolFixation") {
    prefix <- strsplit(sample_name, "_", fixed = TRUE)[[1L]][[1L]]
    if (grepl("^DKMG", prefix)) {
      return("SUM159_2N")
    }
    return(prefix)
  }

  "Unknown"
}

dataset_channel_map <- function(dataset) {
  map_file <- file.path(example_dir, "config", dataset, "channel-map.csv")
  if (!file.exists(map_file)) {
    return(data.frame(raw_channel = character(), standard_channel = character()))
  }
  utils::read.csv(map_file, stringsAsFactors = FALSE, check.names = FALSE)
}

read_dataset_qc <- function(dataset) {
  debug_dir <- file.path(example_dir, "output", dataset, "debug")
  qc_file <- file.path(debug_dir, "qc_table.csv")
  stats_file <- file.path(debug_dir, "population_stats.csv")

  if (!file.exists(qc_file)) {
    stop("Missing QC table: ", qc_file)
  }
  if (!file.exists(stats_file)) {
    stop("Missing population stats table: ", stats_file)
  }

  qc <- utils::read.csv(qc_file, stringsAsFactors = FALSE, check.names = FALSE)
  stats <- utils::read.csv(stats_file, stringsAsFactors = FALSE, check.names = FALSE)
  required_qc_cols <- c("dataset", "sample", "margin_removed_percent", "peacoqc_removed_percent")
  missing_qc_cols <- setdiff(required_qc_cols, colnames(qc))
  if (length(missing_qc_cols) > 0L) {
    stop("QC table is missing columns: ", paste(missing_qc_cols, collapse = ", "))
  }

  root_counts <- stats[stats$pop == "root", c("sample", "count"), drop = FALSE]
  names(root_counts)[names(root_counts) == "count"] <- "root_count"
  keep_counts <- stats[stats$pop == "/PeacoQC_Clean/cells_keep", c("sample", "count"), drop = FALSE]
  names(keep_counts)[names(keep_counts) == "count"] <- "cells_keep_count"
  gating_counts <- merge(root_counts, keep_counts, by = "sample", all.x = TRUE, sort = FALSE)
  gating_counts$gating_removed_percent <- with(
    gating_counts,
    ifelse(root_count > 0, 100 * (root_count - cells_keep_count) / root_count, NA_real_)
  )

  qc <- merge(qc, gating_counts, by = "sample", all.x = TRUE, sort = FALSE)

  rows <- list(
    data.frame(
      dataset = qc$dataset,
      sample = qc$sample,
      cell_line = mapply(infer_cell_line, qc$dataset, qc$sample, USE.NAMES = FALSE),
      step = "Margin removal",
      removed_percent = qc$margin_removed_percent,
      stringsAsFactors = FALSE
    ),
    data.frame(
      dataset = qc$dataset,
      sample = qc$sample,
      cell_line = mapply(infer_cell_line, qc$dataset, qc$sample, USE.NAMES = FALSE),
      step = "PeacoQC",
      removed_percent = qc$peacoqc_removed_percent,
      stringsAsFactors = FALSE
    ),
    data.frame(
      dataset = qc$dataset,
      sample = qc$sample,
      cell_line = mapply(infer_cell_line, qc$dataset, qc$sample, USE.NAMES = FALSE),
      step = "Gating",
      removed_percent = qc$gating_removed_percent,
      stringsAsFactors = FALSE
    )
  )

  do.call(rbind, rows)
}

get_population_count <- function(stats, sample, candidates) {
  sample_stats <- stats[stats$sample == sample, , drop = FALSE]
  for (candidate in candidates) {
    idx <- which(sample_stats$pop == candidate)
    if (length(idx) > 0L) {
      return(sample_stats$count[[idx[[1L]]]])
    }
  }
  NA_real_
}

count_delimited_channels <- function(x) {
  if (is.na(x) || !nzchar(trimws(x))) {
    return(0L)
  }
  length(strsplit(x, ";", fixed = TRUE)[[1L]])
}

sample_fcs_path <- function(dataset, sample, debug_dir) {
  manifest_file <- file.path(debug_dir, "manifest.csv")
  if (file.exists(manifest_file)) {
    manifest <- utils::read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
    idx <- which(manifest$sample_name == sample)
    if (length(idx) > 0L && "fcs_file" %in% colnames(manifest)) {
      return(manifest$fcs_file[[idx[[1L]]]])
    }
  }

  file.path(project_root, "data", dataset, sample)
}

sample_acquisition_seconds <- function(fcs_path) {
  if (!file.exists(fcs_path)) {
    return(data.frame(
      acquisition_time_ticks = NA_real_,
      acquisition_time_step_seconds = NA_real_,
      acquisition_seconds = NA_real_
    ))
  }

  ff <- flowCore::read.FCS(
    fcs_path,
    transformation = FALSE,
    truncate_max_range = FALSE
  )
  expr <- flowCore::exprs(ff)
  if (!"Time" %in% colnames(expr)) {
    return(data.frame(
      acquisition_time_ticks = NA_real_,
      acquisition_time_step_seconds = NA_real_,
      acquisition_seconds = NA_real_
    ))
  }

  time_values <- as.numeric(expr[, "Time"])
  time_values <- time_values[is.finite(time_values)]
  if (length(time_values) < 2L) {
    return(data.frame(
      acquisition_time_ticks = NA_real_,
      acquisition_time_step_seconds = NA_real_,
      acquisition_seconds = NA_real_
    ))
  }

  time_step <- suppressWarnings(as.numeric(flowCore::keyword(ff)[["$TIMESTEP"]]))
  if (length(time_step) == 0L || is.na(time_step) || !is.finite(time_step)) {
    time_step <- 1
  }
  time_ticks <- max(time_values) - min(time_values)

  data.frame(
    acquisition_time_ticks = time_ticks,
    acquisition_time_step_seconds = time_step,
    acquisition_seconds = time_ticks * time_step
  )
}

sample_voltage_table <- function(dataset, sample, cell_line, fcs_path, channel_map) {
  if (!file.exists(fcs_path)) {
    return(data.frame())
  }

  ff <- flowCore::read.FCS(
    fcs_path,
    transformation = FALSE,
    truncate_max_range = FALSE
  )
  kw <- flowCore::keyword(ff)
  par_n <- suppressWarnings(as.integer(kw[["$PAR"]]))
  if (length(par_n) == 0L || is.na(par_n) || par_n <= 0L) {
    return(data.frame())
  }

  rows <- list()
  for (param_idx in seq_len(par_n)) {
    raw_channel <- kw[[paste0("$P", param_idx, "N")]]
    voltage <- suppressWarnings(as.numeric(kw[[paste0("$P", param_idx, "V")]]))
    if (length(raw_channel) == 0L || is.null(raw_channel) || is.na(raw_channel) || raw_channel == "Time") {
      next
    }
    if (length(voltage) == 0L || is.na(voltage) || !is.finite(voltage)) {
      next
    }

    standard_channel <- raw_channel
    map_idx <- match(raw_channel, channel_map$raw_channel)
    if (!is.na(map_idx)) {
      standard_channel <- channel_map$standard_channel[[map_idx]]
    }

    rows[[length(rows) + 1L]] <- data.frame(
      dataset = dataset,
      sample = sample,
      cell_line = cell_line,
      parameter_index = param_idx,
      raw_channel = raw_channel,
      standard_channel = standard_channel,
      voltage = voltage,
      stringsAsFactors = FALSE
    )
  }

  if (!length(rows)) {
    return(data.frame())
  }
  do.call(rbind, rows)
}

read_dataset_event_qc <- function(dataset) {
  debug_dir <- file.path(example_dir, "output", dataset, "debug")
  qc_file <- file.path(debug_dir, "qc_table.csv")
  stats_file <- file.path(debug_dir, "population_stats.csv")

  qc <- utils::read.csv(qc_file, stringsAsFactors = FALSE, check.names = FALSE)
  stats <- utils::read.csv(stats_file, stringsAsFactors = FALSE, check.names = FALSE)

  event_rows <- list()
  acquisition_rows <- list()
  trend_rows <- list()
  voltage_rows <- list()
  channel_map <- dataset_channel_map(dataset)

  for (row_idx in seq_len(nrow(qc))) {
    sample <- qc$sample[[row_idx]]
    cell_line <- infer_cell_line(dataset, sample)
    raw_total <- qc$margin_events_before[[row_idx]]
    cen_count <- get_population_count(
      stats,
      sample,
      c("/PeacoQC_Clean/cen/cen_dna_support", "/PeacoQC_Clean/cen")
    )
    tumor_count <- get_population_count(
      stats,
      sample,
      c(
        "/PeacoQC_Clean/tumor/tumor_fsc_singlets/tumor_singlets/tumor_dna_support",
        "/PeacoQC_Clean/tumor/tumor_dna_support",
        "/PeacoQC_Clean/tumor"
      )
    )

    event_rows[[length(event_rows) + 1L]] <- data.frame(
      dataset = dataset,
      sample = sample,
      cell_line = cell_line,
      metric = event_metric_levels,
      event_count = c(raw_total, cen_count, tumor_count),
      stringsAsFactors = FALSE
    )

    fcs_path <- sample_fcs_path(dataset, sample, debug_dir)
    acquisition_parts <- sample_acquisition_seconds(fcs_path)
    acquisition_rows[[length(acquisition_rows) + 1L]] <- cbind(data.frame(
      dataset = dataset,
      sample = sample,
      cell_line = cell_line,
      fcs_file = fcs_path,
      total_events = raw_total,
      stringsAsFactors = FALSE
    ), acquisition_parts)

    trend_rows[[length(trend_rows) + 1L]] <- data.frame(
      dataset = dataset,
      sample = sample,
      cell_line = cell_line,
      trend = trend_metric_levels,
      flagged_channel_count = c(
        count_delimited_channels(qc$peacoqc_increasing_channels[[row_idx]]),
        count_delimited_channels(qc$peacoqc_decreasing_channels[[row_idx]])
      ),
      stringsAsFactors = FALSE
    )

    voltage_rows[[length(voltage_rows) + 1L]] <- sample_voltage_table(
      dataset = dataset,
      sample = sample,
      cell_line = cell_line,
      fcs_path = fcs_path,
      channel_map = channel_map
    )
  }

  list(
    events = do.call(rbind, event_rows),
    acquisition = do.call(rbind, acquisition_rows),
    trends = do.call(rbind, trend_rows),
    voltages = do.call(rbind, voltage_rows)
  )
}

plot_tbl <- do.call(rbind, lapply(datasets, read_dataset_qc))
plot_tbl$step <- factor(plot_tbl$step, levels = step_levels)
plot_tbl$dataset <- factor(plot_tbl$dataset, levels = datasets)
plot_tbl$cell_line[is.na(plot_tbl$cell_line) | plot_tbl$cell_line == ""] <- "Unknown"

utils::write.csv(plot_tbl, qc_table_file, row.names = FALSE)

plot_height <- 5.2
qc_plot <- ggplot(plot_tbl, aes(x = step, y = removed_percent, group = sample)) +
  geom_line(color = "grey72", linewidth = 0.35, alpha = 0.75) +
  geom_boxplot(
    aes(group = step),
    width = 0.46,
    outlier.shape = NA,
    fill = "grey95",
    color = "grey45",
    linewidth = 0.35
  ) +
  geom_point(
    aes(color = cell_line),
    position = position_jitter(width = 0.12, height = 0, seed = 1),
    size = 2.2,
    alpha = 0.88
  ) +
  facet_grid(cols = vars(dataset), scales = "free_x", space = "free_x") +
  scale_color_manual(values = cell_line_colors, drop = FALSE) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    expand = ggplot2::expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Events removed",
    color = "Inferred cell line"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical"
  )

ggplot2::ggsave(qc_plot_file, qc_plot, width = 10.8, height = plot_height, dpi = 180, bg = "white", limitsize = FALSE)

event_parts <- lapply(datasets, read_dataset_event_qc)
event_tbl <- do.call(rbind, lapply(event_parts, `[[`, "events"))
acquisition_tbl <- do.call(rbind, lapply(event_parts, `[[`, "acquisition"))
trend_tbl <- do.call(rbind, lapply(event_parts, `[[`, "trends"))
voltage_tbl <- do.call(rbind, lapply(event_parts, `[[`, "voltages"))

event_tbl$dataset <- factor(event_tbl$dataset, levels = datasets)
event_tbl$metric <- factor(event_tbl$metric, levels = event_metric_levels)
event_tbl$cell_line[is.na(event_tbl$cell_line) | event_tbl$cell_line == ""] <- "Unknown"
utils::write.csv(event_tbl, event_count_table_file, row.names = FALSE)

event_plot_tbl <- event_tbl[is.finite(event_tbl$event_count), , drop = FALSE]
event_plot <- ggplot(event_plot_tbl, aes(x = metric, y = event_count, group = sample)) +
  geom_line(color = "grey72", linewidth = 0.35, alpha = 0.75, na.rm = TRUE) +
  geom_boxplot(
    aes(group = metric),
    width = 0.46,
    outlier.shape = NA,
    fill = "grey95",
    color = "grey45",
    linewidth = 0.35,
    na.rm = TRUE
  ) +
  geom_point(
    aes(color = cell_line),
    position = position_jitter(width = 0.12, height = 0, seed = 1),
    size = 2.2,
    alpha = 0.88,
    na.rm = TRUE
  ) +
  facet_grid(cols = vars(dataset), scales = "free_x", space = "free_x") +
  scale_color_manual(values = cell_line_colors, drop = FALSE) +
  scale_y_continuous(
    labels = scales::label_comma(),
    expand = ggplot2::expansion(mult = c(0.02, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Events",
    color = "Inferred cell line"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical"
  )
ggplot2::ggsave(event_count_plot_file, event_plot, width = 10.8, height = plot_height, dpi = 180, bg = "white", limitsize = FALSE)

acquisition_tbl$dataset <- factor(acquisition_tbl$dataset, levels = datasets)
acquisition_tbl$cell_line[is.na(acquisition_tbl$cell_line) | acquisition_tbl$cell_line == ""] <- "Unknown"
utils::write.csv(acquisition_tbl, acquisition_time_table_file, row.names = FALSE)

total_events_plot <- ggplot(acquisition_tbl, aes(x = dataset, y = total_events)) +
  geom_boxplot(
    aes(group = dataset),
    width = 0.46,
    outlier.shape = NA,
    fill = "grey95",
    color = "grey45",
    linewidth = 0.35
  ) +
  geom_jitter(aes(color = cell_line), width = 0.18, height = 0, size = 2.4, alpha = 0.88) +
  scale_color_manual(values = cell_line_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::label_comma(), expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
  labs(x = "Dataset", y = "Raw total events", color = "Inferred cell line") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical"
  )
ggplot2::ggsave(file.path(output_dir, "qc_total_events_by_dataset.png"), total_events_plot, width = 8.2, height = 5.4, dpi = 180, bg = "white", limitsize = FALSE)

acquisition_time_plot <- ggplot(acquisition_tbl, aes(x = dataset, y = acquisition_seconds)) +
  geom_boxplot(
    aes(group = dataset),
    width = 0.46,
    outlier.shape = NA,
    fill = "grey95",
    color = "grey45",
    linewidth = 0.35,
    na.rm = TRUE
  ) +
  geom_jitter(aes(color = cell_line), width = 0.18, height = 0, size = 2.4, alpha = 0.88, na.rm = TRUE) +
  scale_color_manual(values = cell_line_colors, drop = FALSE) +
  scale_y_continuous(labels = scales::label_comma(), expand = ggplot2::expansion(mult = c(0.02, 0.08))) +
  labs(x = "Dataset", y = "Acquisition time (seconds)", color = "Inferred cell line") +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical"
  )
ggplot2::ggsave(acquisition_time_plot_file, acquisition_time_plot, width = 8.2, height = 5.4, dpi = 180, bg = "white", limitsize = FALSE)

trend_tbl$dataset <- factor(trend_tbl$dataset, levels = datasets)
trend_tbl$trend <- factor(trend_tbl$trend, levels = trend_metric_levels)
trend_tbl$cell_line[is.na(trend_tbl$cell_line) | trend_tbl$cell_line == ""] <- "Unknown"
utils::write.csv(trend_tbl, trend_count_table_file, row.names = FALSE)

trend_plot <- ggplot(trend_tbl, aes(x = trend, y = flagged_channel_count, group = sample)) +
  geom_line(color = "grey72", linewidth = 0.35, alpha = 0.75) +
  geom_boxplot(
    aes(group = trend),
    width = 0.46,
    outlier.shape = NA,
    fill = "grey95",
    color = "grey45",
    linewidth = 0.35
  ) +
  geom_point(
    aes(color = cell_line),
    position = position_jitter(width = 0.12, height = 0, seed = 1),
    size = 2.2,
    alpha = 0.88
  ) +
  facet_grid(cols = vars(dataset), scales = "free_x", space = "free_x") +
  scale_color_manual(values = cell_line_colors, drop = FALSE) +
  scale_y_continuous(breaks = 0:6, expand = ggplot2::expansion(mult = c(0.02, 0.12))) +
  labs(
    x = NULL,
    y = "PeacoQC trend-flagged channels",
    color = "Inferred cell line"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 25, hjust = 1),
    legend.position = "bottom",
    legend.box = "vertical"
  )
ggplot2::ggsave(trend_count_plot_file, trend_plot, width = 10.8, height = plot_height, dpi = 180, bg = "white", limitsize = FALSE)

if (nrow(voltage_tbl) > 0L) {
  voltage_tbl$dataset <- factor(voltage_tbl$dataset, levels = datasets)
  voltage_tbl$cell_line[is.na(voltage_tbl$cell_line) | voltage_tbl$cell_line == ""] <- "Unknown"
  voltage_tbl$standard_channel <- factor(
    voltage_tbl$standard_channel,
    levels = unique(voltage_tbl$standard_channel[order(voltage_tbl$parameter_index)])
  )
  utils::write.csv(voltage_tbl, voltage_table_file, row.names = FALSE)

  voltage_plot <- ggplot(voltage_tbl, aes(x = standard_channel, y = voltage)) +
    geom_boxplot(
      aes(group = standard_channel),
      width = 0.46,
      outlier.shape = NA,
      fill = "grey95",
      color = "grey45",
      linewidth = 0.35
    ) +
    geom_jitter(aes(color = cell_line), width = 0.16, height = 0, size = 1.9, alpha = 0.82) +
    facet_grid(cols = vars(dataset), scales = "free_x", space = "free_x") +
    scale_color_manual(values = cell_line_colors, drop = FALSE) +
    labs(
      x = NULL,
      y = "FCS detector voltage ($PnV)",
      color = "Inferred cell line"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "bottom",
      legend.box = "vertical"
    )

  ggplot2::ggsave(voltage_plot_file, voltage_plot, width = 12.8, height = 5.6, dpi = 180, bg = "white", limitsize = FALSE)
} else {
  utils::write.csv(data.frame(), voltage_table_file, row.names = FALSE)
}

cat("Wrote QC removal table to ", qc_table_file, "\n", sep = "")
cat("Wrote QC removal plot to ", qc_plot_file, "\n", sep = "")
cat("Wrote QC event count table to ", event_count_table_file, "\n", sep = "")
cat("Wrote QC event count plot to ", event_count_plot_file, "\n", sep = "")
cat("Wrote QC acquisition time table to ", acquisition_time_table_file, "\n", sep = "")
cat("Wrote QC acquisition time plot to ", acquisition_time_plot_file, "\n", sep = "")
cat("Wrote QC PeacoQC trend count table to ", trend_count_table_file, "\n", sep = "")
cat("Wrote QC PeacoQC trend count plot to ", trend_count_plot_file, "\n", sep = "")
cat("Wrote QC voltage table to ", voltage_table_file, "\n", sep = "")
cat("Wrote QC voltage plot to ", voltage_plot_file, "\n", sep = "")
