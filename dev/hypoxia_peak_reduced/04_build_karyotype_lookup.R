source("dev/hypoxia_peak_reduced/hypoxia_peak_reduced_helpers.R")

project_root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE)) >= 1) commandArgs(trailingOnly = TRUE)[[1]] else ".", winslash = "/", mustWork = TRUE)
out_dir <- file.path(project_root, "dev", "hypoxia_peak_reduced", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

lookup_parts <- build_karyotype_flow_lookup(project_root)

write.csv(lookup_parts$lookup_tbl, file.path(out_dir, "karyotype_flow_lookup.csv"), row.names = FALSE)
write.csv(lookup_parts$ploidy_summary, file.path(out_dir, "karyotype_ploidy_summary.csv"), row.names = FALSE)

cat("Saved karyotype lookup outputs under:", out_dir, "\n")
print(lookup_parts$lookup_tbl)
