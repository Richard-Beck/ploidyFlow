args <- commandArgs(trailingOnly = TRUE)

project_root <- normalizePath(
  if (length(args) >= 1L) args[[1]] else ".",
  winslash = "/",
  mustWork = TRUE
)

example_dir <- file.path(project_root, "dev", "opencyto_minimal_example")

suppressPackageStartupMessages({
  library(flowCore)
  library(flowWorkspace)
  library(ggcyto)
  library(ggplot2)
  library(graph)
  library(grid)
  library(openCyto)
  library(patchwork)
  library(png)
})

source(file.path(example_dir, "R", "custom_gates.R"))
source(file.path(example_dir, "R", "channel_mapping.R"))
source(file.path(example_dir, "R", "transformations.R"))
source(file.path(example_dir, "R", "template_validation.R"))
source(file.path(example_dir, "R", "plotting.R"))

resolve_config_path <- function(project_root, path) {
  normalizePath(file.path(project_root, path), winslash = "/", mustWork = FALSE)
}

config <- jsonlite::read_json(
  file.path(example_dir, "config", "config.json"),
  simplifyVector = TRUE
)
config$fcs_path <- resolve_config_path(project_root, config$fcs_path)
config$template_path <- resolve_config_path(project_root, config$template_path)
config$channel_map_path <- resolve_config_path(project_root, config$channel_map_path)
config$plots$output_dir <- resolve_config_path(project_root, config$plots$output_dir)

set.seed(config$random_seed)
register_custom_gates()

channel_map <- load_channel_map(config$channel_map_path)
template_table <- load_template_table(config$template_path)
preprocessing_channels <- template_preprocessing_channels(template_table, config)
margin_channels <- standard_to_raw_channels(preprocessing_channels, channel_map)

fs <- flowCore::read.flowSet(
  files = config$fcs_path,
  transformation = FALSE,
  truncate_max_range = FALSE
)
margin_parts <- remove_margins_flow_data(fs, margin_channels)
fs <- margin_parts$fs
fs <- standardize_channels(fs, channel_map)
config$preprocessing$channels <- preprocessing_channels
preprocessing_parts <- preprocess_flow_data(fs, config)
fs <- preprocessing_parts$fs

template_parts <- load_and_validate_template(config$template_path, fs, config)
fs <- add_template_transform_columns(fs, template_parts$table, config)

gs <- flowWorkspace::GatingSet(fs)
openCyto::gt_gating(template_parts$gt, gs)

stats <- flowWorkspace::gs_pop_get_stats(gs)
print(stats)

plot_files <- render_configured_plots(
  gs,
  template_parts$gt,
  config$plots
)
cat("Wrote configured plot panels:\n")
cat(paste0("  ", plot_files$panel_files, collapse = "\n"), "\n")
cat("Wrote configured plot summary:\n")
cat(paste0("  ", plot_files$summary_file), "\n")
