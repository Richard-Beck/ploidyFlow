# Fixed-init anoxia checkpoint archive

This archive captures the three anoxia-only fixed-initialization optimization
runs that are compared in `workflow/fixed_init_anoxia_checkpoint.Rmd` as of
April 13, 2026.

The archive is intended to be self-contained for report rendering. The R
Markdown checkpoint document reads the fit CSVs, manifests, Stan-data bundle,
and sample metadata from this directory rather than from `processed_data/` or
`figure/`.

This archive serves two purposes:

- it preserves the three model variants currently worth comparing
- it captures the current best explanation of the anoxia data

At present, the preferred model is `observed_cen_scaled`. The main reason is
that its mechanistic derivation based on the observed CEN peak better explains
the covariance between the CEN and tumor G0/G1 peaks. The `dye_ratio` model is
also a major improvement over the baseline `first_pass` model, but both newer
models still appear to use the inferred CEN/tumor ratio to absorb variability
that is not fully explained by the current mechanism.

The group-level effects in `observed_cen_scaled` are still more speculative
than the core observed-CEN-scaling idea, and a version of that model without
group effects remains an important missing comparison.

## Included fits

- `fits/first_pass/fixed_init_optimize-1.csv`
  First-pass ploidy model using `stan/ploidy_histogram_first_pass.stan`.
- `fits/dye_ratio/fixed_init_dye_ratio_optimize-1.csv`
  Dye-ratio model using `stan/ploidy_histogram_dye_ratio.stan`.
- `fits/observed_cen_scaled/fixed_init_dye_ratio_observed_cen_scaled_optimize-1.csv`
  Observed-CEN-scaled model using `stan/ploidy_histogram_dye_ratio_observed_cen_scaled.stan`.

Each fit subdirectory also contains the copied manifest from the original run.

## Experimental caveats

The 2N and 4N controls remain unresolved outliers here, as in the corresponding
hypoxia-style experiments. The main known experimental difference is that these
controls were formaldehyde-fixed separately, while permeabilization and
staining were carried out contemporaneously with the rest of the samples. That
separate fixation step is a plausible source of variability, but the mechanism
is still unresolved.

## Included inputs

- `inputs/stan_data.Rds`
- `inputs/sample_metadata.csv`

These are the inputs currently required by
`workflow/fixed_init_anoxia_checkpoint.Rmd` to reconstruct overlays and
diagnostic summaries from the archived optimizer outputs.

## Reproduction path

Run the canonical fit specs with:

```text
Rscript R/run_fit.R first_pass
Rscript R/run_fit.R dye_ratio
Rscript R/run_fit.R observed_cen_scaled
```

Then archive the resulting fit CSVs/manifests together with the current
`processed_data/stan_data.Rds` and
`processed_data/anoxia-flowcytometry/sample_metadata.csv`.

## Hash manifests

- `hashes/raw_fcs_sha256.csv`
  SHA-256 hashes for the raw anoxia `.fcs` files.
- `hashes/archive_input_and_result_sha256.csv`
  SHA-256 hashes for the archived inputs, fit CSVs, and copied manifests.
- `hashes/source_sha256.csv`
  SHA-256 hashes for the current fit runner/spec files, Stan data preparation
  script, checkpoint R Markdown, and the three Stan model files.
