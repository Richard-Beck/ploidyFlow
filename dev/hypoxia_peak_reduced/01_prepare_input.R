source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

args <- commandArgs(trailingOnly = TRUE)
project_root <- normalizePath(if (length(args) >= 1) args[[1]] else ".", winslash = "/", mustWork = TRUE)
ratio_source <- if (length(args) >= 2) args[[2]] else "fluorescence_sum"
out_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_tbl <- prepare_peak_input(project_root = project_root, threshold = 12000, ratio_source = ratio_source)
date_levels <- attr(input_tbl, "date_levels")

stan_data <- list(
  N = nrow(input_tbl),
  N_date = length(date_levels),
  date_id = input_tbl$date_id,
  x_ratio = input_tbl$ratio_scaled,
  y_cen = input_tbl$cen_peak,
  y_g1 = input_tbl$g1_peak
)

saveRDS(
  list(
    input_tbl = input_tbl,
    stan_data = stan_data,
    threshold = 12000,
    date_levels = date_levels,
    ratio_source = attr(input_tbl, "ratio_source"),
    ratio_center = attr(input_tbl, "ratio_center"),
    ratio_scale = attr(input_tbl, "ratio_scale")
  ),
  file.path(out_dir, "prepared_input.rds")
)

write.csv(input_tbl, file.path(out_dir, "prepared_input_table.csv"), row.names = FALSE)
cat("Saved prepared input under:", out_dir, "\n")
cat("Modeling ratio source:", attr(input_tbl, "ratio_source"), "\n")
