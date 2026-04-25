args <- commandArgs(trailingOnly = TRUE)

project_root <- normalizePath(
  if (length(args) >= 1L) args[[1]] else ".",
  winslash = "/",
  mustWork = TRUE
)

example_dir <- file.path(project_root, "dev", "opencyto_minimal_example")
output_root <- file.path(example_dir, "output")
lab_dir <- file.path(output_root, "gate_lab")
figures_dir <- file.path(lab_dir, "figures")
debug_dir <- file.path(lab_dir, "debug")
max_events_per_sample_gate <- if (length(args) >= 2L) as.integer(args[[2]]) else 5000L

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggplot2)
})

source(file.path(example_dir, "R", "plotting.R"))

tumor_gate_specs <- data.frame(
  gate_spec = c(
    "current",
    "lower_dna_height",
    "wider_fsc",
    "wider_both"
  ),
  fsc_min = c(5.6, 5.6, 5.3, 5.3),
  fsc_max = c(8.2, 8.2, 8.5, 8.5),
  dna_height_min = c(5.3, 5.0, 5.3, 5.0),
  dna_height_max = c(8.2, 8.2, 8.2, 8.4),
  stringsAsFactors = FALSE
)

tumor_dims <- c("FSC-A_ASINH", "DNA_HEIGHT_ASINH")

load_gated_sets <- function(output_root) {
  gated_dirs <- Sys.glob(file.path(output_root, "*", "gated_flow"))
  if (length(gated_dirs) == 0L) {
    stop("No gated_flow directories found under ", output_root)
  }

  lapply(gated_dirs, function(gated_dir) {
    dataset <- basename(dirname(gated_dir))
    list(
      dataset = dataset,
      gated_dir = gated_dir,
      gs = flowWorkspace::load_gs(gated_dir)
    )
  })
}

sample_events <- function(events, n) {
  if (nrow(events) <= n) {
    return(events)
  }
  events[sample.int(nrow(events), n), , drop = FALSE]
}

apply_rect_spec <- function(expr, spec) {
  expr[, tumor_dims[[1]]] >= spec$fsc_min &
    expr[, tumor_dims[[1]]] <= spec$fsc_max &
    expr[, tumor_dims[[2]]] >= spec$dna_height_min &
    expr[, tumor_dims[[2]]] <= spec$dna_height_max
}

collect_gate_lab_data <- function(gated_sets, gate_specs, max_events) {
  plot_rows <- list()
  metric_rows <- list()

  for (set_info in gated_sets) {
    gs <- set_info$gs
    for (sample_name in flowCore::sampleNames(gs)) {
      gh <- gs[[sample_name]]
      parent <- flowWorkspace::gh_pop_get_data(gh, "PeacoQC_Clean")
      expr <- flowCore::exprs(parent)

      missing_dims <- setdiff(tumor_dims, colnames(expr))
      if (length(missing_dims) > 0L) {
        stop(
          "Sample ",
          sample_name,
          " is missing tumor gate dimensions: ",
          paste(missing_dims, collapse = ", ")
        )
      }

      parent_df <- data.frame(
        x = expr[, tumor_dims[[1]]],
        y = expr[, tumor_dims[[2]]]
      )

      for (i in seq_len(nrow(gate_specs))) {
        spec <- gate_specs[i, , drop = FALSE]
        in_gate <- apply_rect_spec(expr, spec)

        plot_df <- parent_df
        plot_df$in_gate <- in_gate
        plot_df <- sample_events(plot_df, max_events)
        plot_df$dataset <- set_info$dataset
        plot_df$sample <- sample_name
        plot_df$gate_spec <- spec$gate_spec
        plot_df$fsc_min <- spec$fsc_min
        plot_df$fsc_max <- spec$fsc_max
        plot_df$dna_height_min <- spec$dna_height_min
        plot_df$dna_height_max <- spec$dna_height_max
        plot_rows[[length(plot_rows) + 1L]] <- plot_df

        metric_rows[[length(metric_rows) + 1L]] <- data.frame(
          dataset = set_info$dataset,
          sample = sample_name,
          gate_spec = spec$gate_spec,
          parent_events = nrow(expr),
          tumor_events = sum(in_gate),
          tumor_fraction = sum(in_gate) / nrow(expr),
          fsc_min = spec$fsc_min,
          fsc_max = spec$fsc_max,
          dna_height_min = spec$dna_height_min,
          dna_height_max = spec$dna_height_max,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  list(
    plot_data = do.call(rbind, plot_rows),
    metrics = do.call(rbind, metric_rows)
  )
}

plot_dataset_gate_grid <- function(plot_data, gate_specs, dataset, output_file) {
  dataset_data <- plot_data[plot_data$dataset == dataset, , drop = FALSE]

  rects <- unique(dataset_data[, c(
    "gate_spec",
    "fsc_min",
    "fsc_max",
    "dna_height_min",
    "dna_height_max"
  )])

  p <- ggplot2::ggplot(
    dataset_data,
    ggplot2::aes(x = x, y = y)
  ) +
    ggplot2::geom_point(
      data = dataset_data[!dataset_data$in_gate, , drop = FALSE],
      color = "grey75",
      alpha = 0.25,
      size = 0.15
    ) +
    ggplot2::geom_point(
      data = dataset_data[dataset_data$in_gate, , drop = FALSE],
      color = "#1B7F79",
      alpha = 0.35,
      size = 0.18
    ) +
    ggplot2::geom_rect(
      data = rects,
      ggplot2::aes(
        xmin = fsc_min,
        xmax = fsc_max,
        ymin = dna_height_min,
        ymax = dna_height_max
      ),
      inherit.aes = FALSE,
      color = "#C74343",
      fill = NA,
      linewidth = 0.35
    ) +
    ggplot2::facet_grid(gate_spec ~ sample) +
    ggplot2::labs(
      title = paste0(dataset, ": tumor gate candidates"),
      x = tumor_dims[[1]],
      y = tumor_dims[[2]]
    ) +
    ggplot2::theme_bw(base_size = 9) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text.x = ggplot2::element_text(size = 7),
      strip.text.y = ggplot2::element_text(size = 8),
      axis.text = ggplot2::element_text(size = 6)
    )

  width <- max(8, length(unique(dataset_data$sample)) * 3.1)
  height <- max(6, nrow(gate_specs) * 2.1)
  ggplot2::ggsave(output_file, p, width = width, height = height, dpi = 150, bg = "white")
  output_file
}

set.seed(1)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(debug_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(figures_dir, pattern = "[.]png$", full.names = TRUE))

gated_sets <- load_gated_sets(output_root)
lab_data <- collect_gate_lab_data(
  gated_sets = gated_sets,
  gate_specs = tumor_gate_specs,
  max_events = max_events_per_sample_gate
)

metrics_file <- file.path(debug_dir, "tumor_gate_metrics.csv")
utils::write.csv(lab_data$metrics, metrics_file, row.names = FALSE)

for (dataset in unique(lab_data$plot_data$dataset)) {
  output_file <- file.path(figures_dir, paste0(sanitize_filename(dataset), "_tumor_gate_grid.png"))
  plot_dataset_gate_grid(lab_data$plot_data, tumor_gate_specs, dataset, output_file)
}

cat("Wrote tumor gate lab figures:\n")
cat(paste0("  ", list.files(figures_dir, pattern = "[.]png$", full.names = TRUE), collapse = "\n"), "\n")
cat("Wrote tumor gate metrics:\n")
cat("  ", metrics_file, "\n", sep = "")
