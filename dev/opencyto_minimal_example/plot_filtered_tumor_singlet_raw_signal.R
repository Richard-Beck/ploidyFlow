args <- commandArgs(trailingOnly = TRUE)

project_root <- normalizePath(
  if (length(args) >= 1L) args[[1]] else ".",
  winslash = "/",
  mustWork = TRUE
)

example_dir <- file.path(project_root, "dev", "opencyto_minimal_example")
consistency_dir <- file.path(example_dir, "output", "peak_consistency")
raw_signal_dir <- file.path(consistency_dir, "filtered_tumor_singlet_raw_signal")

dir.create(raw_signal_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggplot2)
})

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
    return(matches[[1L]])
  }

  NA_character_
}

raw_population_xy <- function(gh, population, x_channel = "DNA_AREA", y_channel = "DNA_HEIGHT") {
  events <- flowWorkspace::gh_pop_get_data(gh, population)
  missing_channels <- setdiff(c(x_channel, y_channel), colnames(events))
  if (length(missing_channels)) {
    stop("Missing channel(s) in population '", population, "': ", paste(missing_channels, collapse = ", "))
  }

  values <- flowCore::exprs(events)[, c(x_channel, y_channel), drop = FALSE]
  out <- data.frame(
    dna_area = as.numeric(values[, x_channel]),
    dna_height = as.numeric(values[, y_channel])
  )
  out[is.finite(out$dna_area) & is.finite(out$dna_height), , drop = FALSE]
}

unlink(list.files(raw_signal_dir, pattern = "[.]png$", full.names = TRUE), force = TRUE)

raw_signal_plot_count <- 0L
gated_dirs <- discover_gated_dirs(file.path(example_dir, "output"))
if (!length(gated_dirs)) {
  stop("No gated_flow directories found under ", file.path(example_dir, "output"))
}

for (gated_dir in gated_dirs) {
  dataset <- basename(dirname(gated_dir))
  gs <- flowWorkspace::load_gs(gated_dir)
  tumor_support_pop <- resolve_population(gs, "tumor_dna_support")
  if (is.na(tumor_support_pop)) {
    warning("Skipping ", dataset, ": missing tumor_dna_support population.")
    next
  }

  for (sample_name in flowCore::sampleNames(gs)) {
    event_tbl <- raw_population_xy(gs[[sample_name]], tumor_support_pop)
    if (!nrow(event_tbl)) {
      next
    }

    p_raw_signal <- ggplot(event_tbl, aes(x = dna_area, y = dna_height)) +
      geom_point(size = 0.18, alpha = 0.18, color = "black") +
      scale_x_continuous(labels = scales::label_comma()) +
      scale_y_continuous(labels = scales::label_comma()) +
      labs(
        title = paste0(dataset, ": ", sample_name),
        subtitle = paste0("tumor_dna_support raw signal | n = ", scales::comma(nrow(event_tbl))),
        x = "Raw DNA-A",
        y = "Raw DNA-H"
      ) +
      theme_bw(base_size = 10) +
      theme(
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 10),
        plot.subtitle = element_text(size = 8)
      )

    raw_signal_file <- file.path(
      raw_signal_dir,
      paste0(sanitize_filename(dataset), "__", sanitize_filename(sample_name), ".png")
    )
    ggplot2::ggsave(raw_signal_file, p_raw_signal, width = 6.6, height = 5.4, dpi = 150, bg = "white", limitsize = FALSE)
    raw_signal_plot_count <- raw_signal_plot_count + 1L
  }
}

cat("Wrote ", raw_signal_plot_count, " filtered tumor singlet raw signal plot(s) to ", raw_signal_dir, "\n", sep = "")
