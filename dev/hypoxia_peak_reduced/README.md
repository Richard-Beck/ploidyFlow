# Hypoxia Peak Reduced Checkpoint

This folder is a checkpointed, refactored copy of the baseline hypoxia peak workflow.

The workflow is intentionally split into three steps:

1. `01_prepare_input.R`
   Builds the one-row-per-sample peak input table and the Stan data.
2. `02_run_model.R`
   Fits the original `hypoxia_peak_reduced.stan` model and writes the fitted per-sample predictions plus fitted parameters.
3. `03_render_outputs.R`
   Produces the checkpoint outputs to inspect first:
   - `sample_peak_summary.csv`
   - `phi_counterfactual.png`

## Current interpretation

The interpretation being preserved at this checkpoint is:

- the `Predicted if phiF = 1` peak locations carry real DNA-content signal across samples;
- this `phi = 1` counterfactual is the main object of interest in this workflow;
- `sample_peak_summary.csv` is the key table for reviewing those per-sample peak locations alongside the fitted parameters and outlier flag.

## Key files

- `hypoxia_peak_reduced.stan`: baseline Stan model
- `output/sample_peak_summary.csv`: key one-row-per-sample table
- `output/fit_parameters.csv`: fitted model parameters
- `output/date_phi.csv`: fitted date-level `phi`
- `output/phi_counterfactual.png`: main checkpoint figure

## How to rerun

From the repository root:

```powershell
Rscript dev\hypoxia_peak_reduced\01_prepare_input.R
Rscript dev\hypoxia_peak_reduced\02_run_model.R
Rscript dev\hypoxia_peak_reduced\03_render_outputs.R
```

All outputs are written under `dev/hypoxia_peak_reduced/output`.
