# DNA Content from Flow Cytometry

## Objective

Infer latent DNA-content structure from flow cytometry histograms, with enough
model structure to explain both peak positions and population fractions across
related samples.

The current focus is the anoxia dataset, where a shared CEN spike-in provides a
useful internal reference across samples.

## Current Checkpoint

This repository now has a preserved anoxia-only checkpoint under
`results/2026-04-13-fixed-init-anoxia/`. That archive contains the three model
variants currently worth comparing:

- `first_pass`: baseline model with globally fixed peak locations
- `dye_ratio`: introduces a dye-ratio mechanism to explain coordinated peak
  shifts
- `observed_cen_scaled`: preferred current model; uses a more mechanistic peak
  derivation tied to the observed CEN location and includes small group-level
  effects

These three fits are compared in `workflow/fixed_init_anoxia_checkpoint.Rmd`.
That report is intended to render from the archived `results/` checkpoint
without depending on transient `processed_data/` or `figure/` outputs.

## Current Interpretation

The baseline `first_pass` model is useful as a reference, but it does not
explain the clear covariance between the CEN peak and the tumor G0/G1 peak.

The `dye_ratio` model is a major improvement because it explains much more of
the peak-position variation through a shared mechanism. However, it still uses
the inferred CEN/tumor ratio as its main degree of freedom for explaining peak
location shifts, so it can distort that ratio in conditions where additional
sources of variability are present.

The `observed_cen_scaled` model is the current preferred explanation of the
anoxia data. Its main strength is the mechanistic derivation based on the
observed CEN location, which better captures the covariance between the CEN and
tumor peaks. The group-level effects in that model are more speculative than
the core CEN-scaling idea, and a version without group effects is still a
useful missing comparison.

## Experimental Caveats

The 2N and 4N controls remain unresolved outliers in both the anoxia and
hypoxia-style datasets. The main known experimental difference is that these
controls were formaldehyde-fixed separately, whereas permeabilization and
staining were carried out contemporaneously with the other samples. That
separate fixation step is a plausible source of variability, but the mechanism
is not yet established.

More broadly, the current models still suggest that the observed peak shifts are
not explained entirely by the modeled dye-ratio / CEN-scaling mechanism. There
are likely additional sources of technical or biological variability that are
not yet represented.

## Main Entry Points

- `workflow/data_prep.Rmd`: gating, export, and exploratory data-preparation
  workflow
- `workflow/fixed_init_anoxia_checkpoint.Rmd`: archived three-model comparison
  for the anoxia checkpoint
- `R/run_fit.R`: canonical fit runner
- `R/fit_specs.R`: canonical preserved fit specifications

Run the current preserved fits with:

```text
Rscript R/run_fit.R first_pass
Rscript R/run_fit.R dye_ratio
Rscript R/run_fit.R observed_cen_scaled
```

## Repository Structure

- `data/`: raw flow cytometry inputs
- `processed_data/`: working derived data and non-archival intermediate outputs
- `figure/`: working figure outputs
- `results/`: archived checkpoint artifacts intended to be preserved
- `stan/`: Stan model definitions
- `R/`: reusable preparation, fitting, and plotting code
- `workflow/`: R Markdown analysis and reporting documents
