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
tumor_2_1_ratio_plot_file <- file.path(consistency_dir, "tumor_2_1_peak_ratio_by_dataset.png")
tumor_1_cen_1_ratio_plot_file <- file.path(consistency_dir, "tumor_1_cen_1_peak_ratio_by_dataset.png")
tumor_1_cen_1_ratio_file <- file.path(consistency_dir, "tumor_1_cen_1_peak_ratio_by_dataset.csv")
lower_peak_ratio_plot_file <- file.path(consistency_dir, "lower_peak_location_vs_adjacent_ratio.png")
fwhm_cv_file <- file.path(consistency_dir, "agreed_tumor_peak_fwhm_cv.csv")
fwhm_cv_plot_file <- file.path(consistency_dir, "agreed_tumor_peak_fwhm_cv_by_sample.png")
tumor_1_fwhm_cv_file <- file.path(consistency_dir, "tumor_1_peak_fwhm_cv_by_dataset.csv")
tumor_1_fwhm_cv_plot_file <- file.path(consistency_dir, "tumor_1_peak_fwhm_cv_by_dataset.png")
cen_ratio_plot_file <- file.path(consistency_dir, "cen_peak_ratio_by_dataset.png")
cen_cv_file <- file.path(consistency_dir, "cen_largest_peak_fwhm_cv.csv")
cen_cv_plot_file <- file.path(consistency_dir, "cen_largest_peak_fwhm_cv_by_dataset.png")
annotated_peak_file <- file.path(peak_dir, "agreed_tumor_peaks_by_sample.csv")

dir.create(consistency_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
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
    return(matches[[1L]])
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
  values[is.finite(values) & values > 0]
}

interpolate_halfmax_crossing <- function(x1, y1, x2, y2, half_height) {
  if (!is.finite(y2 - y1) || abs(y2 - y1) < .Machine$double.eps) {
    return(mean(c(x1, x2)))
  }
  x1 + (half_height - y1) * (x2 - x1) / (y2 - y1)
}

estimate_peak_fwhm <- function(values, agreed_peak_x, density_n = 2048L) {
  values <- values[is.finite(values) & values > 0]
  if (length(values) < 50L || length(unique(values)) < 3L || !is.finite(agreed_peak_x) || agreed_peak_x <= 0) {
    return(data.frame(
      density_peak_x = NA_real_,
      peak_density = NA_real_,
      fwhm = NA_real_,
      cv_fwhm = NA_real_
    ))
  }

  dens <- stats::density(values, n = density_n, from = min(values), to = max(values), na.rm = TRUE)
  local_max_idx <- which(dens$y >= c(-Inf, head(dens$y, -1L)) & dens$y >= c(tail(dens$y, -1L), -Inf))
  if (!length(local_max_idx)) {
    local_max_idx <- which.max(dens$y)
  }

  peak_idx <- local_max_idx[which.min(abs(dens$x[local_max_idx] - agreed_peak_x))]
  half_height <- dens$y[[peak_idx]] / 2

  left_candidates <- which(dens$y[seq_len(peak_idx)] <= half_height)
  right_candidates <- which(dens$y[peak_idx:length(dens$y)] <= half_height)
  if (!length(left_candidates) || !length(right_candidates)) {
    return(data.frame(
      density_peak_x = dens$x[[peak_idx]],
      peak_density = dens$y[[peak_idx]],
      fwhm = NA_real_,
      cv_fwhm = NA_real_
    ))
  }

  left_idx <- max(left_candidates)
  right_idx <- peak_idx + min(right_candidates) - 1L
  if (left_idx >= peak_idx || right_idx <= peak_idx) {
    return(data.frame(
      density_peak_x = dens$x[[peak_idx]],
      peak_density = dens$y[[peak_idx]],
      fwhm = NA_real_,
      cv_fwhm = NA_real_
    ))
  }

  left_x <- interpolate_halfmax_crossing(
    dens$x[[left_idx]],
    dens$y[[left_idx]],
    dens$x[[left_idx + 1L]],
    dens$y[[left_idx + 1L]],
    half_height
  )
  right_x <- interpolate_halfmax_crossing(
    dens$x[[right_idx - 1L]],
    dens$y[[right_idx - 1L]],
    dens$x[[right_idx]],
    dens$y[[right_idx]],
    half_height
  )
  fwhm <- right_x - left_x

  data.frame(
    density_peak_x = dens$x[[peak_idx]],
    peak_density = dens$y[[peak_idx]],
    fwhm = fwhm,
    cv_fwhm = fwhm / (2.355 * agreed_peak_x)
  )
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

tumor_peak_idx <- if ("branch" %in% names(peaks)) {
  peaks$branch == "Tumor"
} else {
  grepl("^Tumor", peaks$peak_type)
}
cen_peak_idx <- if ("branch" %in% names(peaks)) {
  peaks$branch == "CEN"
} else {
  grepl("^CEN", peaks$peak_type)
}

plot_tbl <- peaks[
  tumor_peak_idx &
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
      lower_peak_x = z$peak_x[[i]],
      ratio = z$peak_x[[i + 1L]] / z$peak_x[[i]],
      regression_group = ifelse(
        z$dataset[[1L]] == "Polyploidization_EthanolFixation",
        "Ethanol fixation",
        "Hypoxia/anoxia"
      ),
      stringsAsFactors = FALSE
    )
  }
}

ratio_tbl <- if (length(ratio_parts)) do.call(rbind, ratio_parts) else data.frame()
if (nrow(ratio_tbl)) {
  ratio_tbl$sample_y <- factor(ratio_tbl$sample_y, levels = levels(sample_tbl$sample_y))
  ratio_tbl$ratio_type <- factor(ratio_tbl$ratio_type, levels = c("Tumor 2 / Tumor 1", "Tumor 3 / Tumor 2"))
  ratio_tbl$regression_group <- factor(ratio_tbl$regression_group, levels = c("Ethanol fixation", "Hypoxia/anoxia"))

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

  tumor_2_1_ratio_tbl <- ratio_tbl[ratio_tbl$ratio_type == "Tumor 2 / Tumor 1", , drop = FALSE]
  if (nrow(tumor_2_1_ratio_tbl)) {
    tumor_2_1_ratio_plot <- ggplot(tumor_2_1_ratio_tbl, aes(x = dataset, y = ratio)) +
      geom_hline(yintercept = 2, linewidth = 0.45, linetype = "dashed", color = "grey30") +
      geom_hline(yintercept = c(1.95, 2.05), linewidth = 0.35, linetype = "dashed", color = "grey55") +
      geom_jitter(aes(color = cell_line), width = 0.18, height = 0, size = 2.5, alpha = 0.85) +
      scale_color_manual(values = cell_line_colors, drop = FALSE) +
      labs(
        x = "Dataset",
        y = "Tumor 2 / Tumor 1 DNA-A ratio",
        color = "Inferred cell line"
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "bottom",
        legend.box = "vertical"
      )

    ggplot2::ggsave(tumor_2_1_ratio_plot_file, tumor_2_1_ratio_plot, width = 8.2, height = 5.4, dpi = 180, bg = "white", limitsize = FALSE)
  }

  tumor_1_tbl <- plot_tbl[plot_tbl$plot_peak_type == "Tumor 1", , drop = FALSE]
  cen_1_tbl <- peaks[cen_peak_idx & peaks$peak_type == "CEN 1", , drop = FALSE]
  if (nrow(tumor_1_tbl) && nrow(cen_1_tbl)) {
    cen_1_tbl <- cen_1_tbl[, c("dataset", "sample_name", "peak_x"), drop = FALSE]
    names(cen_1_tbl)[names(cen_1_tbl) == "peak_x"] <- "cen_1_peak_x"
    tumor_1_cen_1_ratio_tbl <- merge(
      tumor_1_tbl[, c("dataset", "sample_name", "cell_line", "qc_flag", "qc_reason", "peak_x"), drop = FALSE],
      cen_1_tbl,
      by = c("dataset", "sample_name"),
      all = FALSE,
      sort = FALSE
    )
    names(tumor_1_cen_1_ratio_tbl)[names(tumor_1_cen_1_ratio_tbl) == "peak_x"] <- "tumor_1_peak_x"
    tumor_1_cen_1_ratio_tbl$ratio <- tumor_1_cen_1_ratio_tbl$tumor_1_peak_x / tumor_1_cen_1_ratio_tbl$cen_1_peak_x

    if (nrow(tumor_1_cen_1_ratio_tbl)) {
      utils::write.csv(tumor_1_cen_1_ratio_tbl, tumor_1_cen_1_ratio_file, row.names = FALSE)

      tumor_1_cen_1_ratio_plot <- ggplot(tumor_1_cen_1_ratio_tbl, aes(x = dataset, y = ratio)) +
        geom_jitter(aes(color = cell_line), width = 0.18, height = 0, size = 2.5, alpha = 0.85) +
        scale_color_manual(values = cell_line_colors, drop = FALSE) +
        labs(
          x = "Dataset",
          y = "Tumor 1 / CEN 1 DNA-A ratio",
          color = "Inferred cell line"
        ) +
        theme_bw(base_size = 10) +
        theme(
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 25, hjust = 1),
          legend.position = "bottom",
          legend.box = "vertical"
        )

      ggplot2::ggsave(
        tumor_1_cen_1_ratio_plot_file,
        tumor_1_cen_1_ratio_plot,
        width = 8.2,
        height = 5.4,
        dpi = 180,
        bg = "white",
        limitsize = FALSE
      )
    }
  }

  lower_peak_ratio_plot <- ggplot(ratio_tbl, aes(x = lower_peak_x, y = ratio)) +
    geom_hline(yintercept = 2, linewidth = 0.45, linetype = "dashed", color = "grey30") +
    geom_smooth(
      data = ratio_tbl[!ratio_tbl$qc_flag, , drop = FALSE],
      aes(linetype = regression_group),
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      color = "black",
      linewidth = 0.7
    ) +
    geom_point(aes(color = cell_line, shape = dataset), size = 2.6, alpha = 0.9) +
    geom_point(
      data = ratio_tbl[ratio_tbl$qc_flag, , drop = FALSE],
      aes(x = lower_peak_x, y = ratio),
      inherit.aes = FALSE,
      shape = 21,
      size = 4.5,
      stroke = 1.2,
      color = "#B2182B",
      fill = NA
    ) +
    scale_x_continuous(
      labels = scales::label_comma(),
      expand = ggplot2::expansion(mult = c(0.04, 0.06))
    ) +
    scale_color_manual(values = cell_line_colors, drop = FALSE) +
    labs(
      x = "Lower agreed peak raw DNA-A",
      y = "Adjacent agreed peak DNA-A ratio",
      color = "Inferred cell line",
      shape = "Dataset",
      linetype = "Regression fit"
    ) +
    theme_bw(base_size = 10) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical"
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))

  ggplot2::ggsave(lower_peak_ratio_plot_file, lower_peak_ratio_plot, width = 8.2, height = 5.8, dpi = 180, bg = "white", limitsize = FALSE)
}

fwhm_parts <- list()
gated_dirs <- discover_gated_dirs(file.path(example_dir, "output"))
for (gated_dir in gated_dirs) {
  dataset <- basename(dirname(gated_dir))
  gs <- flowWorkspace::load_gs(gated_dir)
  keep_pop <- required_population(gs, "cells_keep")

  for (sample_name in flowCore::sampleNames(gs)) {
    peak_rows <- plot_tbl[plot_tbl$dataset == dataset & plot_tbl$sample_name == sample_name, , drop = FALSE]
    if (!nrow(peak_rows)) {
      next
    }

    keep_values <- raw_dna_values(gs[[sample_name]], keep_pop)
    for (i in seq_len(nrow(peak_rows))) {
      width_tbl <- estimate_peak_fwhm(keep_values, peak_rows$peak_x[[i]])
      fwhm_parts[[length(fwhm_parts) + 1L]] <- cbind(
        peak_rows[i, c("dataset", "sample_name", "cell_line", "sample_y", "qc_flag", "qc_reason", "peak_type", "plot_peak_type", "peak_x"), drop = FALSE],
        n_events_cells_keep = length(keep_values),
        width_tbl
      )
    }
  }
}

fwhm_tbl <- if (length(fwhm_parts)) do.call(rbind, fwhm_parts) else data.frame()
utils::write.csv(fwhm_tbl, fwhm_cv_file, row.names = FALSE)
if (nrow(fwhm_tbl)) {
  fwhm_tbl$sample_y <- factor(fwhm_tbl$sample_y, levels = levels(sample_tbl$sample_y))
  fwhm_tbl$plot_peak_type <- factor(fwhm_tbl$plot_peak_type, levels = levels(plot_tbl$plot_peak_type))

  fwhm_cv_plot <- ggplot(fwhm_tbl, aes(x = cv_fwhm, y = sample_y)) +
    geom_point(aes(color = cell_line, shape = plot_peak_type), size = 2.4, alpha = 0.9, na.rm = TRUE) +
    geom_point(
      data = fwhm_tbl[fwhm_tbl$qc_flag, , drop = FALSE],
      aes(x = cv_fwhm, y = sample_y),
      inherit.aes = FALSE,
      shape = 21,
      size = 4.1,
      stroke = 1.2,
      color = "#B2182B",
      fill = NA,
      na.rm = TRUE
    ) +
    facet_grid(rows = vars(dataset), scales = "free_y", space = "free_y") +
    scale_x_continuous(
      labels = scales::label_percent(accuracy = 0.1),
      expand = ggplot2::expansion(mult = c(0.02, 0.06))
    ) +
    scale_y_discrete(labels = axis_labels) +
    scale_color_manual(values = cell_line_colors, drop = FALSE) +
    labs(
      x = "FWHM-derived CV proxy",
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
    )

  ggplot2::ggsave(fwhm_cv_plot_file, fwhm_cv_plot, width = 11.5, height = plot_height, dpi = 180, bg = "white", limitsize = FALSE)

  tumor_1_fwhm_tbl <- fwhm_tbl[fwhm_tbl$plot_peak_type == "Tumor 1", , drop = FALSE]
  if (nrow(tumor_1_fwhm_tbl)) {
    utils::write.csv(tumor_1_fwhm_tbl, tumor_1_fwhm_cv_file, row.names = FALSE)

    tumor_1_fwhm_cv_plot <- ggplot(tumor_1_fwhm_tbl, aes(x = dataset, y = cv_fwhm)) +
      geom_boxplot(
        aes(group = dataset),
        width = 0.46,
        outlier.shape = NA,
        fill = "grey95",
        color = "grey45",
        linewidth = 0.35,
        na.rm = TRUE
      ) +
      geom_jitter(aes(color = cell_line), width = 0.18, height = 0, size = 2.5, alpha = 0.85, na.rm = TRUE) +
      scale_y_continuous(
        labels = scales::label_percent(accuracy = 0.1),
        expand = ggplot2::expansion(mult = c(0.02, 0.08))
      ) +
      scale_color_manual(values = cell_line_colors, drop = FALSE) +
      labs(
        x = "Dataset",
        y = "Tumor 1 FWHM-derived CV proxy",
        color = "Inferred cell line"
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "bottom",
        legend.box = "vertical"
      )

    ggplot2::ggsave(
      tumor_1_fwhm_cv_plot_file,
      tumor_1_fwhm_cv_plot,
      width = 8.2,
      height = 5.4,
      dpi = 180,
      bg = "white",
      limitsize = FALSE
    )
  }
} else {
  utils::write.csv(data.frame(), tumor_1_fwhm_cv_file, row.names = FALSE)
}

cen_tbl <- peaks[cen_peak_idx, , drop = FALSE]
if (nrow(cen_tbl)) {
  cen_tbl$cell_line <- mapply(infer_cell_line, cen_tbl$dataset, cen_tbl$sample_name, USE.NAMES = FALSE)
  cen_tbl$cell_line[is.na(cen_tbl$cell_line) | cen_tbl$cell_line == ""] <- "Unknown"
  cen_tbl <- merge(cen_tbl, qc_flags, by = c("dataset", "sample_name"), all.x = TRUE, sort = FALSE)
  cen_tbl$qc_flag[is.na(cen_tbl$qc_flag)] <- FALSE
  cen_tbl$qc_reason[is.na(cen_tbl$qc_reason)] <- ""

  cen_ratio_parts <- list()
  for (key in unique(paste(cen_tbl$dataset, cen_tbl$sample_name, sep = "\r"))) {
    z <- cen_tbl[paste(cen_tbl$dataset, cen_tbl$sample_name, sep = "\r") == key, , drop = FALSE]
    z <- z[order(z$peak_x), , drop = FALSE]
    if (nrow(z) < 2L) {
      next
    }
    cen_ratio_parts[[length(cen_ratio_parts) + 1L]] <- data.frame(
      dataset = z$dataset[[1L]],
      sample_name = z$sample_name[[1L]],
      cell_line = z$cell_line[[1L]],
      qc_flag = z$qc_flag[[1L]],
      qc_reason = z$qc_reason[[1L]],
      cen1_peak_x = z$peak_x[[1L]],
      cen2_peak_x = z$peak_x[[2L]],
      ratio = z$peak_x[[2L]] / z$peak_x[[1L]],
      stringsAsFactors = FALSE
    )
  }

  cen_ratio_tbl <- if (length(cen_ratio_parts)) do.call(rbind, cen_ratio_parts) else data.frame()
  if (nrow(cen_ratio_tbl)) {
    cen_ratio_plot <- ggplot(cen_ratio_tbl, aes(x = dataset, y = ratio)) +
      geom_hline(yintercept = 2, linewidth = 0.45, linetype = "dashed", color = "grey30") +
      geom_hline(yintercept = c(1.95, 2.05), linewidth = 0.35, linetype = "dashed", color = "grey55") +
      geom_jitter(aes(color = cell_line), width = 0.18, height = 0, size = 2.5, alpha = 0.85) +
      scale_color_manual(values = cell_line_colors, drop = FALSE) +
      labs(
        x = "Dataset",
        y = "CEN 2 / CEN 1 DNA-A ratio",
        color = "Inferred cell line"
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "bottom",
        legend.box = "vertical"
      )

    ggplot2::ggsave(cen_ratio_plot_file, cen_ratio_plot, width = 8.2, height = 5.4, dpi = 180, bg = "white", limitsize = FALSE)
  }

  cen_cv_parts <- list()
  for (gated_dir in gated_dirs) {
    dataset <- basename(dirname(gated_dir))
    gs <- flowWorkspace::load_gs(gated_dir)
    cen_pop <- optional_population(gs, "cen_dna_support")
    if (is.null(cen_pop)) {
      next
    }

    for (sample_name in flowCore::sampleNames(gs)) {
      z <- cen_tbl[cen_tbl$dataset == dataset & cen_tbl$sample_name == sample_name, , drop = FALSE]
      if (!nrow(z)) {
        next
      }
      z <- z[order(z$peak_x, decreasing = TRUE), , drop = FALSE]
      cen_values <- raw_dna_values(gs[[sample_name]], cen_pop)
      width_tbl <- estimate_peak_fwhm(cen_values, z$peak_x[[1L]])
      cen_cv_parts[[length(cen_cv_parts) + 1L]] <- cbind(
        z[1L, c("dataset", "sample_name", "cell_line", "qc_flag", "qc_reason", "peak_type", "peak_x"), drop = FALSE],
        n_events_cen = length(cen_values),
        width_tbl
      )
    }
  }

  cen_cv_tbl <- if (length(cen_cv_parts)) do.call(rbind, cen_cv_parts) else data.frame()
  utils::write.csv(cen_cv_tbl, cen_cv_file, row.names = FALSE)
  if (nrow(cen_cv_tbl)) {
    cen_cv_plot <- ggplot(cen_cv_tbl, aes(x = dataset, y = cv_fwhm)) +
      geom_hline(yintercept = 0.03, linewidth = 0.45, linetype = "dashed", color = "grey30") +
      geom_jitter(aes(color = cell_line), width = 0.18, height = 0, size = 2.5, alpha = 0.85, na.rm = TRUE) +
      scale_y_continuous(labels = scales::label_percent(accuracy = 0.1)) +
      scale_color_manual(values = cell_line_colors, drop = FALSE) +
      labs(
        x = "Dataset",
        y = "Largest CEN peak FWHM-derived CV proxy",
        color = "Inferred cell line"
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1),
        legend.position = "bottom",
        legend.box = "vertical"
      )

    ggplot2::ggsave(cen_cv_plot_file, cen_cv_plot, width = 8.2, height = 5.4, dpi = 180, bg = "white", limitsize = FALSE)
  }
} else {
  utils::write.csv(data.frame(), cen_cv_file, row.names = FALSE)
}

cat("Wrote agreed tumor peak plot to ", output_file, "\n", sep = "")
cat("Wrote QC-flagged agreed tumor peak plot to ", consistency_peak_file, "\n", sep = "")
cat("Wrote QC-flagged tumor peak ratio plot to ", ratio_plot_file, "\n", sep = "")
cat("Wrote Tumor 2 / Tumor 1 peak ratio plot to ", tumor_2_1_ratio_plot_file, "\n", sep = "")
cat("Wrote Tumor 1 / CEN 1 peak ratio table to ", tumor_1_cen_1_ratio_file, "\n", sep = "")
cat("Wrote Tumor 1 / CEN 1 peak ratio plot to ", tumor_1_cen_1_ratio_plot_file, "\n", sep = "")
cat("Wrote lower peak location vs adjacent ratio plot to ", lower_peak_ratio_plot_file, "\n", sep = "")
cat("Wrote FWHM-derived agreed tumor peak CV table to ", fwhm_cv_file, "\n", sep = "")
cat("Wrote FWHM-derived agreed tumor peak CV plot to ", fwhm_cv_plot_file, "\n", sep = "")
cat("Wrote Tumor 1 FWHM-derived CV table to ", tumor_1_fwhm_cv_file, "\n", sep = "")
cat("Wrote Tumor 1 FWHM-derived CV plot to ", tumor_1_fwhm_cv_plot_file, "\n", sep = "")
cat("Wrote CEN peak ratio plot to ", cen_ratio_plot_file, "\n", sep = "")
cat("Wrote largest CEN peak FWHM-derived CV table to ", cen_cv_file, "\n", sep = "")
cat("Wrote largest CEN peak FWHM-derived CV plot to ", cen_cv_plot_file, "\n", sep = "")
cat("Wrote annotated tumor peak table to ", annotated_peak_file, "\n", sep = "")
