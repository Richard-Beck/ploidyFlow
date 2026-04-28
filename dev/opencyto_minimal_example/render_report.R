report_file <- file.path("dev", "opencyto_minimal_example", "report_opencyto_minimal_example.Rmd")

if (!file.exists(report_file)) {
  stop("Cannot find report source: ", report_file, call. = FALSE)
}

rmarkdown::render(
  report_file,
  output_dir = dirname(report_file),
  intermediates_dir = dirname(report_file),
  envir = new.env(parent = globalenv())
)
