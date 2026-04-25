sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9_.-]+", "_", x)
}

gt_method_dims <- function(gt_method) {
  dims <- slot(gt_method, "dims")
  trimws(strsplit(dims, ",", fixed = TRUE)[[1]])
}

pop_alias <- function(pop_path) {
  if (pop_path == "root") {
    return("root")
  }
  sub("^.*/", "", pop_path)
}

resolve_population_path <- function(gt, population) {
  nodes <- graph::nodes(gt)
  if (population %in% nodes) {
    return(population)
  }

  matches <- nodes[vapply(nodes, pop_alias, character(1)) == population]
  if (length(matches) != 1L) {
    stop("Could not uniquely resolve population '", population, "'.")
  }

  matches[[1]]
}

edge_data_or_null <- function(gt, parent, child) {
  graph::edgeData(gt, from = parent, to = child)[[1]]
}

canonical_parent_edge <- function(gt, child) {
  template_edges <- graph::edges(gt)

  for (parent in names(template_edges)) {
    if (!child %in% template_edges[[parent]]) {
      next
    }

    edge_data <- edge_data_or_null(gt, parent, child)
    if (!isTRUE(edge_data$isReference)) {
      return(list(parent = parent, child = child, edge_data = edge_data))
    }
  }

  stop("Could not find a drawable parent edge for population '", child, "'.")
}

image_file_to_wrapped_plot <- function(image_file) {
  patchwork::wrap_elements(
    full = grid::rasterGrob(png::readPNG(image_file), interpolate = TRUE)
  )
}

plot_panels <- function(plots_config) {
  panels <- plots_config$panels
  if (is.data.frame(panels)) {
    return(lapply(seq_len(nrow(panels)), function(i) as.list(panels[i, , drop = FALSE])))
  }
  panels
}

plot_title <- function(panel, fallback) {
  if (!is.null(panel$title) && !is.na(panel$title) && nzchar(panel$title)) {
    return(panel$title)
  }
  fallback
}

render_gate_step_plot <- function(gs, gt, panel, sample_index) {
  gh <- gs[[sample_index]]
  sample_name <- flowCore::sampleNames(gs)[[sample_index]]
  child <- resolve_population_path(gt, panel$population)
  parent_edge <- canonical_parent_edge(gt, child)
  gt_method <- parent_edge$edge_data$gtMethod

  if (slot(gt_method, "name") == "boolGate") {
    stop("Population '", panel$population, "' is a boolGate and cannot be drawn as a gate_step panel.")
  }

  dims <- gt_method_dims(gt_method)
  parent_events <- flowWorkspace::gh_pop_get_data(gh, parent_edge$parent)
  gate <- flowWorkspace::gh_pop_get_gate(gh, child)

  if (length(dims) == 1L) {
    p <- ggcyto::autoplot(parent_events, x = dims[[1]]) +
      ggcyto::geom_gate(gate)
  } else if (length(dims) == 2L) {
    p <- ggcyto::autoplot(parent_events, x = dims[[1]], y = dims[[2]]) +
      ggcyto::geom_gate(gate)
  } else {
    stop(
      "Cannot draw population '",
      panel$population,
      "' as gate_step: expected one or two dimensions, got ",
      length(dims),
      "."
    )
  }

  p + ggplot2::ggtitle(plot_title(
    panel,
    paste0(
      sample_name,
      ": ",
      pop_alias(parent_edge$parent),
      " -> ",
      pop_alias(child),
      " (",
      slot(gt_method, "name"),
      ")"
    )
  ))
}

render_population_density_plot <- function(gs, gt, panel, sample_index) {
  gh <- gs[[sample_index]]
  population <- resolve_population_path(gt, panel$population)
  events <- flowWorkspace::gh_pop_get_data(gh, population)

  if (!panel$x %in% colnames(events)) {
    stop("Panel '", panel$id, "' references missing x channel '", panel$x, "'.")
  }

  ggcyto::autoplot(events, x = panel$x) +
    ggplot2::ggtitle(plot_title(panel, paste0(pop_alias(population), ": ", panel$x)))
}

render_population_scatter_plot <- function(gs, gt, panel, sample_index) {
  gh <- gs[[sample_index]]
  population <- resolve_population_path(gt, panel$population)
  events <- flowWorkspace::gh_pop_get_data(gh, population)

  missing_channels <- setdiff(c(panel$x, panel$y), colnames(events))
  if (length(missing_channels) > 0L) {
    stop("Panel '", panel$id, "' references missing channels: ", paste(missing_channels, collapse = ", "))
  }

  ggcyto::autoplot(events, x = panel$x, y = panel$y) +
    ggplot2::ggtitle(plot_title(panel, paste0(pop_alias(population), ": ", panel$x, " vs ", panel$y)))
}

render_base_graph_panel <- function(gs, gt, panel, output_file) {
  png(output_file, width = 1100, height = 800)
  on.exit(dev.off(), add = TRUE)

  if (panel$type == "gating_template_graph") {
    plot(gt)
  } else if (panel$type == "gating_set_graph") {
    plot(gs)
  } else {
    stop("Unsupported base graph panel type: ", panel$type)
  }

  dev.off()
  on.exit(NULL, add = FALSE)
  output_file
}

render_gg_panel <- function(plot, output_file, width = 7, height = 5, dpi = 150) {
  ggplot2::ggsave(output_file, plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  output_file
}

render_plot_panel <- function(gs, gt, panel, output_dir, sample_index, dpi) {
  output_file <- file.path(output_dir, paste0(sanitize_filename(panel$id), ".png"))

  if (panel$type %in% c("gating_template_graph", "gating_set_graph")) {
    return(render_base_graph_panel(gs, gt, panel, output_file))
  }

  p <- switch(
    panel$type,
    gate_step = render_gate_step_plot(gs, gt, panel, sample_index),
    population_density = render_population_density_plot(gs, gt, panel, sample_index),
    population_scatter = render_population_scatter_plot(gs, gt, panel, sample_index),
    stop("Unsupported plot panel type: ", panel$type)
  )

  render_gg_panel(p, output_file, dpi = dpi)
}

validate_plot_config <- function(plots_config) {
  required_layout <- c("design", "width", "height", "dpi")
  missing_layout <- setdiff(required_layout, names(plots_config$layout))
  if (length(missing_layout) > 0L) {
    stop("plots.layout is missing: ", paste(missing_layout, collapse = ", "))
  }

  required_panel <- c("id", "slot", "type")
  panels <- plot_panels(plots_config)
  for (panel in panels) {
    missing_panel <- setdiff(required_panel, names(panel))
    if (length(missing_panel) > 0L) {
      stop("Plot panel is missing: ", paste(missing_panel, collapse = ", "))
    }
    if (nchar(panel$slot) != 1L) {
      stop("Panel '", panel$id, "' slot must be a single character.")
    }
  }

  design_slots <- unique(strsplit(paste0(plots_config$layout$design, collapse = ""), "")[[1]])
  design_slots <- setdiff(design_slots, ".")
  panel_slots <- vapply(panels, function(panel) panel$slot, character(1))

  missing_slots <- setdiff(design_slots, panel_slots)
  if (length(missing_slots) > 0L) {
    stop("Layout uses slots without panels: ", paste(missing_slots, collapse = ", "))
  }

  unused_slots <- setdiff(panel_slots, design_slots)
  if (length(unused_slots) > 0L) {
    stop("Panels use slots absent from layout: ", paste(unused_slots, collapse = ", "))
  }

  invisible(TRUE)
}

render_configured_plots <- function(gs, gt, plots_config, sample_index = 1L) {
  validate_plot_config(plots_config)

  output_dir <- plots_config$output_dir
  summary_file <- file.path(output_dir, plots_config$summary_file)
  panel_dir <- file.path(
    if (!is.null(plots_config$debug_dir) && nzchar(plots_config$debug_dir)) {
      plots_config$debug_dir
    } else {
      output_dir
    },
    "panels",
    tools::file_path_sans_ext(basename(summary_file))
  )
  dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(list.files(panel_dir, pattern = "[.]png$", full.names = TRUE))

  panel_files <- list()
  wrapped_panels <- list()

  for (panel in plot_panels(plots_config)) {
    panel_file <- render_plot_panel(gs, gt, panel, panel_dir, sample_index, plots_config$layout$dpi)
    panel_files[[panel$id]] <- panel_file
    wrapped_panels[[panel$slot]] <- image_file_to_wrapped_plot(panel_file)
  }

  summary <- do.call(
    patchwork::wrap_plots,
    c(wrapped_panels, list(design = paste(plots_config$layout$design, collapse = "\n")))
  )

  ggplot2::ggsave(
    summary_file,
    plot = summary,
    width = plots_config$layout$width,
    height = plots_config$layout$height,
    dpi = plots_config$layout$dpi,
    bg = "white"
  )

  if (isTRUE(plots_config$cleanup_panels)) {
    unlink(unlist(panel_files, use.names = FALSE), force = TRUE)
    unlink(panel_dir, recursive = TRUE, force = TRUE)
    panel_parent <- dirname(panel_dir)
    if (dir.exists(panel_parent) && length(list.files(panel_parent, all.files = FALSE)) == 0L) {
      unlink(panel_parent, recursive = TRUE, force = TRUE)
    }
  }

  list(
    summary_file = summary_file,
    panel_files = unlist(panel_files, use.names = FALSE)
  )
}
