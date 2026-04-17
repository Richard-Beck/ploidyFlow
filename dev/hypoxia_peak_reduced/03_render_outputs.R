source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

project_root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[[1]] else ".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output")
report_path <- file.path(project_root, "dev", "hypoxia_peak_reduced", "report_hypoxia_peak_reduced.Rmd")

summary_tbl <- render_fit_outputs(out_dir)

if (rmarkdown::pandoc_available("1.12.3")) {
  rmarkdown::render(
    report_path,
    output_dir = dirname(report_path),
    intermediates_dir = dirname(report_path),
    envir = new.env(parent = globalenv())
  )
} else {
  knitr::knit(
    report_path,
    output = file.path(dirname(report_path), "report_hypoxia_peak_reduced.md"),
    envir = new.env(parent = globalenv())
  )
  cat("Pandoc was not available; knitted markdown report instead of HTML.\n")
}

cat("Saved checkpoint outputs under:", out_dir, "\n")
print(summary_tbl)
