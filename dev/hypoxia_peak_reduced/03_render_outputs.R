source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

project_root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[[1]] else ".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output")

fit_tbl <- read.csv(file.path(out_dir, "sample_predictions.csv"), stringsAsFactors = FALSE, check.names = FALSE)
summary_tbl <- build_sample_peak_summary(fit_tbl)

write.csv(summary_tbl, file.path(out_dir, "sample_peak_summary.csv"), row.names = FALSE)
ggsave(file.path(out_dir, "phi_counterfactual.png"), plot_phi_counterfactual(fit_tbl), width = 12, height = 8, dpi = 220)

cat("Saved checkpoint outputs under:", out_dir, "\n")
print(summary_tbl)
