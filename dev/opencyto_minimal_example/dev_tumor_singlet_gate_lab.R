args <- commandArgs(trailingOnly = TRUE)

project_root <- normalizePath(
  if (length(args) >= 1L) args[[1]] else ".",
  winslash = "/",
  mustWork = TRUE
)

example_dir <- file.path(project_root, "dev", "opencyto_minimal_example")
output_root <- file.path(example_dir, "output")
lab_dir <- file.path(output_root, "gate_lab", "tumor_singlet")
figures_dir <- file.path(lab_dir, "figures")
debug_dir <- file.path(lab_dir, "debug")
max_events_per_sample_gate <- if (length(args) >= 2L) as.integer(args[[2]]) else 5000L
max_fit_events_per_sample_gate <- if (length(args) >= 3L) as.integer(args[[3]]) else 20000L

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggplot2)
})

source(file.path(example_dir, "R", "plotting.R"))
source(file.path(example_dir, "R", "custom_gates.R"))

singlet_gate_specs <- data.frame(
  gate_spec = c(
    "current",
    "narrow_90",
    "wide_98",
    "less_trim",
    "k2",
    "k4",
    "slope_06_13"
  ),
  slope_min = c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4, 0.6),
  slope_max = c(1.4, 1.4, 1.4, 1.4, 1.4, 1.4, 1.3),
  slope_steps = c(81L, 81L, 81L, 81L, 81L, 81L, 81L),
  x_trim = c(0.01, 0.01, 0.01, 0.005, 0.01, 0.01, 0.01),
  K = c(3L, 3L, 3L, 3L, 2L, 4L, 3L),
  central_prob = c(0.95, 0.90, 0.98, 0.95, 0.95, 0.95, 0.95),
  nstart = c(10L, 10L, 10L, 10L, 10L, 10L, 10L),
  stringsAsFactors = FALSE
)

singlet_dims <- c("DNA_AREA_ASINH", "DNA_HEIGHT_ASINH")

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

build_residual_band_gate <- function(fr, spec, max_fit_events) {
  expr <- flowCore::exprs(fr)
  fit_expr <- expr
  if (!is.na(max_fit_events) && nrow(fit_expr) > max_fit_events) {
    fit_expr <- fit_expr[sample.int(nrow(fit_expr), max_fit_events), , drop = FALSE]
  }

  fit <- residual_band_fit(
    x = fit_expr[, singlet_dims[[1]]],
    y = fit_expr[, singlet_dims[[2]]],
    slope_min = spec$slope_min,
    slope_max = spec$slope_max,
    slope_steps = spec$slope_steps,
    x_trim = spec$x_trim,
    K = spec$K,
    central_prob = spec$central_prob,
    nstart = spec$nstart
  )

  x_limits <- range(expr[, singlet_dims[[1]]], finite = TRUE)
  x_min <- x_limits[[1]]
  x_max <- x_limits[[2]]
  r_lo <- fit$residual_limits[[1]]
  r_hi <- fit$residual_limits[[2]]
  slope <- fit$slope

  boundaries <- matrix(
    c(
      x_min, slope * x_min + r_lo,
      x_max, slope * x_max + r_lo,
      x_max, slope * x_max + r_hi,
      x_min, slope * x_min + r_hi
    ),
    ncol = 2,
    byrow = TRUE
  )
  colnames(boundaries) <- singlet_dims

  list(
    gate = flowCore::polygonGate(.gate = boundaries, filterId = "residualBandGate"),
    fit = fit,
    fit_events = nrow(fit_expr)
  )
}

collect_singlet_lab_data <- function(gated_sets, gate_specs, max_events, max_fit_events) {
  plot_rows <- list()
  polygon_rows <- list()
  metric_rows <- list()

  for (set_info in gated_sets) {
    gs <- set_info$gs
    for (sample_name in flowCore::sampleNames(gs)) {
      gh <- gs[[sample_name]]
      parent <- flowWorkspace::gh_pop_get_data(gh, "tumor")
      expr <- flowCore::exprs(parent)

      missing_dims <- setdiff(singlet_dims, colnames(expr))
      if (length(missing_dims) > 0L) {
        stop(
          "Sample ",
          sample_name,
          " is missing tumor singlet dimensions: ",
          paste(missing_dims, collapse = ", ")
        )
      }

      parent_df <- data.frame(
        x = expr[, singlet_dims[[1]]],
        y = expr[, singlet_dims[[2]]]
      )

      for (i in seq_len(nrow(gate_specs))) {
        spec <- gate_specs[i, , drop = FALSE]
        gate_parts <- tryCatch(
          build_residual_band_gate(parent, spec, max_fit_events),
          error = function(e) e
        )

        if (inherits(gate_parts, "error")) {
          metric_rows[[length(metric_rows) + 1L]] <- data.frame(
            dataset = set_info$dataset,
            sample = sample_name,
            gate_spec = spec$gate_spec,
            status = "error",
            error = conditionMessage(gate_parts),
            parent_events = nrow(expr),
            fit_events = NA_integer_,
            singlet_events = NA_integer_,
            singlet_fraction = NA_real_,
            slope = NA_real_,
            score = NA_real_,
            fit_singlet_count = NA_integer_,
            spec[, setdiff(names(spec), "gate_spec"), drop = FALSE],
            stringsAsFactors = FALSE
          )
          next
        }

        in_gate <- flowCore::filter(parent, gate_parts$gate)@subSet

        plot_df <- parent_df
        plot_df$in_gate <- in_gate
        plot_df <- sample_events(plot_df, max_events)
        plot_df$dataset <- set_info$dataset
        plot_df$sample <- sample_name
        plot_df$gate_spec <- spec$gate_spec
        plot_rows[[length(plot_rows) + 1L]] <- plot_df

        polygon_df <- data.frame(
          x = gate_parts$gate@boundaries[, singlet_dims[[1]]],
          y = gate_parts$gate@boundaries[, singlet_dims[[2]]],
          dataset = set_info$dataset,
          sample = sample_name,
          gate_spec = spec$gate_spec,
          vertex = seq_len(nrow(gate_parts$gate@boundaries))
        )
        polygon_rows[[length(polygon_rows) + 1L]] <- polygon_df

        metric_rows[[length(metric_rows) + 1L]] <- data.frame(
          dataset = set_info$dataset,
          sample = sample_name,
          gate_spec = spec$gate_spec,
          status = "ok",
          error = "",
          parent_events = nrow(expr),
          fit_events = gate_parts$fit_events,
          singlet_events = sum(in_gate),
          singlet_fraction = sum(in_gate) / nrow(expr),
          slope = gate_parts$fit$slope,
          score = gate_parts$fit$score,
          fit_singlet_count = gate_parts$fit$singlet_count,
          spec[, setdiff(names(spec), "gate_spec"), drop = FALSE],
          stringsAsFactors = FALSE
        )
      }
    }
  }

  list(
    plot_data = do.call(rbind, plot_rows),
    polygon_data = do.call(rbind, polygon_rows),
    metrics = do.call(rbind, metric_rows)
  )
}

plot_dataset_singlet_grid <- function(plot_data, polygon_data, gate_specs, dataset, output_file) {
  dataset_data <- plot_data[plot_data$dataset == dataset, , drop = FALSE]
  dataset_polygons <- polygon_data[polygon_data$dataset == dataset, , drop = FALSE]

  p <- ggplot2::ggplot(dataset_data, ggplot2::aes(x = x, y = y)) +
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
    ggplot2::geom_polygon(
      data = dataset_polygons,
      ggplot2::aes(x = x, y = y, group = interaction(gate_spec, sample)),
      inherit.aes = FALSE,
      color = "#C74343",
      fill = NA,
      linewidth = 0.35
    ) +
    ggplot2::facet_grid(gate_spec ~ sample) +
    ggplot2::labs(
      title = paste0(dataset, ": tumor singlet gate candidates"),
      x = singlet_dims[[1]],
      y = singlet_dims[[2]]
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
lab_data <- collect_singlet_lab_data(
  gated_sets = gated_sets,
  gate_specs = singlet_gate_specs,
  max_events = max_events_per_sample_gate,
  max_fit_events = max_fit_events_per_sample_gate
)

metrics_file <- file.path(debug_dir, "tumor_singlet_gate_metrics.csv")
utils::write.csv(lab_data$metrics, metrics_file, row.names = FALSE)

for (dataset in unique(lab_data$plot_data$dataset)) {
  output_file <- file.path(figures_dir, paste0(sanitize_filename(dataset), "_tumor_singlet_gate_grid.png"))
  plot_dataset_singlet_grid(
    lab_data$plot_data,
    lab_data$polygon_data,
    singlet_gate_specs,
    dataset,
    output_file
  )
}

cat("Wrote tumor singlet gate lab figures:\n")
cat(paste0("  ", list.files(figures_dir, pattern = "[.]png$", full.names = TRUE), collapse = "\n"), "\n")
cat("Wrote tumor singlet gate metrics:\n")
cat("  ", metrics_file, "\n", sep = "")
