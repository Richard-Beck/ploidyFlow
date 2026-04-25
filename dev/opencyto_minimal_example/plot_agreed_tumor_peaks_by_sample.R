args <- commandArgs(trailingOnly = TRUE)

project_root <- normalizePath(
  if (length(args) >= 1L) args[[1]] else ".",
  winslash = "/",
  mustWork = TRUE
)

example_dir <- file.path(project_root, "dev", "opencyto_minimal_example")
peak_dir <- file.path(example_dir, "output", "peak_detection_figures")
consistency_dir <- file.path(example_dir, "output", "peak_consistency")
agreed_peak_file <- file.path(peak_dir, "agreed_peaks.csv")
output_file <- file.path(peak_dir, "agreed_tumor_peaks_by_sample.png")
consistency_peak_file <- file.path(consistency_dir, "agreed_tumor_peaks_by_sample.png")
ratio_plot_file <- file.path(consistency_dir, "tumor_peak_ratio_by_sample.png")
annotated_peak_file <- file.path(peak_dir, "agreed_tumor_peaks_by_sample.csv")

dir.create(consistency_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggtext)
})

if (!file.exists(agreed_peak_file)) {
  stop("Missing agreed peak file: ", agreed_peak_file)
}

infer_cell_line <- function(dataset, sample_name) {
  if (dataset == "Anoxia_FlowCytometry") {
    if (grepl("^Sample_CEN(?:\\.fcs|\\.fcs\\.h5)?$", sample_name)) {
      return(NA_character_)
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
    return(NA_character_)
  }

  if (dataset == "Polyploidization_EthanolFixation") {
    prefix <- strsplit(sample_name, "_", fixed = TRUE)[[1L]][[1L]]
    if (grepl("^DKMG", prefix)) {
      return("SUM159_2N")
    }
    return(prefix)
  }

  NA_character_
}

escape_html <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

cell_line_colors <- c(
  SUM159_2N = "#0072B2",
  SUM159_4N = "#D55E00",
  HCC1937 = "#009E73",
  Jurkat = "#CC79A7",
  KALS1 = "#E69F00",
  MDAMB231 = "#56B4E9",
  Unknown = "#666666"
)

getcol <- function(df, nm) {
  if (nm %in% names(df)) {
    return(df[[nm]])
  }
  rep(NA_real_, nrow(df))
}

read_qc_flags <- function(example_dir) {
  datasets <- c("Anoxia_FlowCytometry", "Hypoxia_SUM159", "Polyploidization_EthanolFixation")
  qc_parts <- list()

  for (dataset in datasets) {
    debug_dir <- file.path(example_dir, "output", dataset, "debug")
    qc_file <- file.path(debug_dir, "qc_table.csv")
    pop_file <- file.path(debug_dir, "population_stats.csv")
    if (!file.exists(qc_file) || !file.exists(pop_file)) {
      next
    }

    qc <- utils::read.csv(qc_file, stringsAsFactors = FALSE)
    pop <- utils::read.csv(pop_file, stringsAsFactors = FALSE)
    pop_wide <- stats::reshape(pop, idvar = "sample", timevar = "pop", direction = "wide")
    qc <- merge(qc, pop_wide, by = "sample", all.x = TRUE)

    qc$cells_keep_n <- getcol(qc, "count./PeacoQC_Clean/cells_keep")
    qc$qc_flag <- qc$cells_keep_n < 7000 | qc$fsc_gate_qc != "ok"
    qc$qc_reason <- ifelse(
      qc$cells_keep_n < 7000,
      "low retained event count",
      ifelse(qc$fsc_gate_qc != "ok", "FSC gate QC", "")
    )

    qc_parts[[length(qc_parts) + 1L]] <- qc[, c("dataset", "sample", "qc_flag", "qc_reason"), drop = FALSE]
  }

  if (!length(qc_parts)) {
    return(data.frame(dataset = character(), sample_name = character(), qc_flag = logical(), qc_reason = character()))
  }

  out <- do.call(rbind, qc_parts)
  names(out)[names(out) == "sample"] <- "sample_name"
  out
}

peaks <- utils::read.csv(agreed_peak_file, stringsAsFactors = FALSE)
required_cols <- c("dataset", "sample_name", "peak_x", "peak_type")
missing_cols <- setdiff(required_cols, names(peaks))
if (length(missing_cols)) {
  stop("Missing required column(s): ", paste(missing_cols, collapse = ", "))
}

plot_tbl <- peaks[
  grepl("^Tumor", peaks$peak_type) &
    !(peaks$dataset == "Anoxia_FlowCytometry" & grepl("^Sample_CEN", peaks$sample_name)),
  ,
  drop = FALSE
]

plot_tbl$cell_line <- mapply(infer_cell_line, plot_tbl$dataset, plot_tbl$sample_name, USE.NAMES = FALSE)
plot_tbl$cell_line[is.na(plot_tbl$cell_line) | plot_tbl$cell_line == ""] <- "Unknown"

extra_cell_lines <- setdiff(sort(unique(plot_tbl$cell_line)), names(cell_line_colors))
if (length(extra_cell_lines)) {
  hue_cols <- grDevices::hcl.colors(length(extra_cell_lines), palette = "Dark 3")
  names(hue_cols) <- extra_cell_lines
  cell_line_colors <- c(cell_line_colors, hue_cols)
}

sample_tbl <- unique(plot_tbl[, c("dataset", "sample_name", "cell_line")])
sample_tbl <- sample_tbl[order(sample_tbl$dataset, sample_tbl$cell_line, sample_tbl$sample_name), , drop = FALSE]
sample_tbl$sample_y <- paste(sample_tbl$dataset, sample_tbl$sample_name, sep = "__")
sample_tbl$sample_y <- factor(sample_tbl$sample_y, levels = rev(sample_tbl$sample_y))
sample_tbl$axis_label <- sprintf(
  "<span style='color:%s'>%s</span>",
  cell_line_colors[sample_tbl$cell_line],
  escape_html(sample_tbl$sample_name)
)

plot_tbl <- merge(plot_tbl, sample_tbl, by = c("dataset", "sample_name", "cell_line"), all.x = TRUE, sort = FALSE)
plot_tbl$sample_y <- factor(plot_tbl$sample_y, levels = levels(sample_tbl$sample_y))

plot_tbl <- plot_tbl[order(plot_tbl$dataset, plot_tbl$sample_name, plot_tbl$peak_x), , drop = FALSE]
sample_key <- paste(plot_tbl$dataset, plot_tbl$sample_name, sep = "\r")
plot_tbl$plot_peak_type <- NA_character_
for (key in unique(sample_key)) {
  idx <- which(sample_key == key)
  idx <- idx[order(plot_tbl$peak_x[idx])]
  plot_tbl$plot_peak_type[idx] <- paste0("Tumor ", seq_along(idx))
}
plot_tbl$plot_peak_type <- factor(plot_tbl$plot_peak_type, levels = sort(unique(plot_tbl$plot_peak_type)))

qc_flags <- read_qc_flags(example_dir)
sample_tbl <- merge(sample_tbl, qc_flags, by = c("dataset", "sample_name"), all.x = TRUE, sort = FALSE)
sample_tbl$qc_flag[is.na(sample_tbl$qc_flag)] <- FALSE
sample_tbl$qc_reason[is.na(sample_tbl$qc_reason)] <- ""
sample_tbl$axis_label <- ifelse(
  sample_tbl$qc_flag,
  sprintf(
    "<span style='color:#B2182B; font-weight:700'>[QC] %s</span>",
    escape_html(sample_tbl$sample_name)
  ),
  sprintf(
    "<span style='color:%s'>%s</span>",
    cell_line_colors[sample_tbl$cell_line],
    escape_html(sample_tbl$sample_name)
  )
)

plot_tbl <- merge(
  plot_tbl,
  sample_tbl[, c("dataset", "sample_name", "qc_flag", "qc_reason"), drop = FALSE],
  by = c("dataset", "sample_name"),
  all.x = TRUE,
  sort = FALSE
)
plot_tbl$qc_flag[is.na(plot_tbl$qc_flag)] <- FALSE

axis_labels <- sample_tbl$axis_label
names(axis_labels) <- as.character(sample_tbl$sample_y)

p <- ggplot(plot_tbl, aes(x = peak_x, y = sample_y)) +
  geom_point(aes(color = cell_line, shape = plot_peak_type), size = 2.4, alpha = 0.9) +
  geom_point(
    data = plot_tbl[plot_tbl$qc_flag, , drop = FALSE],
    aes(x = peak_x, y = sample_y),
    inherit.aes = FALSE,
    shape = 21,
    size = 4.1,
    stroke = 1.2,
    color = "#B2182B",
    fill = NA
  ) +
  facet_grid(rows = vars(dataset), scales = "free_y", space = "free_y") +
  scale_x_continuous(
    labels = scales::label_comma(),
    expand = ggplot2::expansion(mult = c(0.02, 0.06))
  ) +
  scale_y_discrete(labels = axis_labels) +
  scale_color_manual(values = cell_line_colors, breaks = names(cell_line_colors)) +
  labs(
    x = "Raw DNA-A",
    y = "Sample",
    color = "Inferred cell line",
    shape = "Agreed tumor peak rank"
  ) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    strip.text.y = element_text(angle = 0),
    axis.text.y = ggtext::element_markdown(size = 6.7),
    legend.position = "bottom",
    legend.box = "vertical"
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))

plot_height <- max(7, 0.17 * nrow(sample_tbl) + 2.4)
ggplot2::ggsave(output_file, p, width = 13, height = plot_height, dpi = 180, bg = "white", limitsize = FALSE)
ggplot2::ggsave(consistency_peak_file, p, width = 13, height = plot_height, dpi = 180, bg = "white", limitsize = FALSE)

utils::write.csv(
  plot_tbl[, c("dataset", "sample_name", "cell_line", "qc_flag", "qc_reason", "peak_type", "plot_peak_type", "peak_x")],
  annotated_peak_file,
  row.names = FALSE
)

ratio_parts <- list()
for (sample_key_i in unique(plot_tbl$sample_y)) {
  z <- plot_tbl[plot_tbl$sample_y == sample_key_i, , drop = FALSE]
  z <- z[order(z$peak_x), , drop = FALSE]
  if (nrow(z) < 2L) {
    next
  }

  for (i in seq_len(nrow(z) - 1L)) {
    ratio_parts[[length(ratio_parts) + 1L]] <- data.frame(
      dataset = z$dataset[[1L]],
      sample_name = z$sample_name[[1L]],
      cell_line = z$cell_line[[1L]],
      sample_y = z$sample_y[[1L]],
      qc_flag = z$qc_flag[[1L]],
      qc_reason = z$qc_reason[[1L]],
      ratio_type = paste0(z$plot_peak_type[[i + 1L]], " / ", z$plot_peak_type[[i]]),
      ratio = z$peak_x[[i + 1L]] / z$peak_x[[i]],
      stringsAsFactors = FALSE
    )
  }
}

ratio_tbl <- if (length(ratio_parts)) do.call(rbind, ratio_parts) else data.frame()
if (nrow(ratio_tbl)) {
  ratio_tbl$sample_y <- factor(ratio_tbl$sample_y, levels = levels(sample_tbl$sample_y))
  ratio_tbl$ratio_type <- factor(ratio_tbl$ratio_type, levels = c("Tumor 2 / Tumor 1", "Tumor 3 / Tumor 2"))

  ratio_plot <- ggplot(ratio_tbl, aes(x = ratio, y = sample_y)) +
    geom_vline(xintercept = 2, linewidth = 0.45, linetype = "dashed", color = "grey30") +
    geom_point(aes(color = cell_line, shape = ratio_type), size = 2.3, alpha = 0.9) +
    geom_point(
      data = ratio_tbl[ratio_tbl$qc_flag, , drop = FALSE],
      aes(x = ratio, y = sample_y),
      inherit.aes = FALSE,
      shape = 21,
      size = 4.2,
      stroke = 1.2,
      color = "#B2182B",
      fill = NA
    ) +
    facet_grid(rows = vars(dataset), scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = axis_labels) +
    scale_color_manual(values = cell_line_colors, drop = FALSE) +
    labs(
      x = "Adjacent agreed peak DNA-A ratio",
      y = "Sample",
      color = "Inferred cell line",
      shape = "Ratio"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.major.y = element_line(color = "grey88", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      strip.text.y = element_text(angle = 0),
      axis.text.y = ggtext::element_markdown(size = 6.7),
      legend.position = "bottom",
      legend.box = "vertical"
    )

  ggplot2::ggsave(ratio_plot_file, ratio_plot, width = 11.5, height = plot_height, dpi = 180, bg = "white", limitsize = FALSE)
}

cat("Wrote agreed tumor peak plot to ", output_file, "\n", sep = "")
cat("Wrote QC-flagged agreed tumor peak plot to ", consistency_peak_file, "\n", sep = "")
cat("Wrote QC-flagged tumor peak ratio plot to ", ratio_plot_file, "\n", sep = "")
cat("Wrote annotated tumor peak table to ", annotated_peak_file, "\n", sep = "")
