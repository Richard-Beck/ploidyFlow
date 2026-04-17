# Hypoxia Peak Reduced Checkpoint

This folder is a checkpointed, refactored copy of the baseline hypoxia peak
workflow.

The workflow is intentionally split into eight steps:

1. `01_prepare_input.R`
   Builds the one-row-per-sample peak input table and the Stan data.
   Optionally accepts a ratio source argument: `fluorescence_sum` (default) or
   `cell_count`.
2. `02_run_model.R`
   Fits the baseline `hypoxia_peak_reduced.stan` model and writes the fitted
   per-sample predictions plus fitted parameters.
3. `03_render_outputs.R`
   Produces the first-pass outputs to inspect:
   - `sample_peak_summary.csv`
   - `phi_counterfactual.png`
4. `04_build_karyotype_lookup.R`
   Builds the approved karyotyping-to-flow lookup table plus per-ID DNA ploidy
   summaries.
5. `05_run_ablations.R`
   Runs the baseline model plus the first ablation batch and writes a
   comparison table.
6. `06_run_extensions.R`
   Runs non-ablation model extensions and writes a comparison table.
7. `07_run_nuts_extension.R`
   Runs a NUTS fit for one selected extension, currently most useful for the
   `delta_dna` and `adjacent_8n_delta_dna` variants.
8. `08_render_nuts_counterfactual.R`
   Reconstructs posterior summaries for the delta-only counterfactual peaks
   from a completed NUTS fit.

## Current interpretation

The interpretation being preserved at this checkpoint is:

- date-level peak shifts are real and remain necessary in this workflow;
- the burden proxy is not currently the main explanatory result;
- a pure adjacent-state reinterpretation is too aggressive on its own, because
  it relabels some late 4N samples implausibly;
- the most interesting current hypothesis is the combined
  `adjacent_8n_delta_dna` model, where a smaller subset of ambiguous samples,
  especially in the O2 lineage, can be interpreted as an adjacent higher-ploidy
  state plus sample-specific DNA scaling;
- the karyotype-to-flow lookup is part of that interpretation, because the O2
  lineage looks more compatible with a hyperploid explanation than the cleaner
  control-like samples.

This is still a development checkpoint rather than a preserved final result.
The current goal is to record the best-supported working interpretation, not to
claim that the combined adjacent-state plus `delta_dna` story is already
settled.

## Key files

- `hypoxia_peak_reduced.stan`: baseline Stan model
- `output/sample_peak_summary.csv`: baseline one-row-per-sample summary
- `output/fit_parameters.csv`: baseline fitted model parameters
- `output/date_phi.csv`: fitted date-level `phi`
- `output/phi_counterfactual.png`: baseline checkpoint figure
- `output/karyotype_flow_lookup.csv`: karyotyping-to-flow lookup with DNA
  ploidy summary
- `output/karyotype_ploidy_summary.csv`: per-karyotyping-ID DNA ploidy summary
- `output/ablations/ablation_comparison.csv`: first-pass ablation comparison
  metrics
- `output/extensions/extension_comparison.csv`: comparison metrics for model
  extensions
- `output/extensions/*/sample_peak_summary.csv`: per-extension state calls and
  counterfactual peak summaries
- `output/extensions/*/nuts/delta_only_counterfactual_summary.csv`: posterior
  summary for the delta-only counterfactual peaks after NUTS

## How to rerun

From the repository root:

```powershell
Rscript dev\hypoxia_peak_reduced\01_prepare_input.R
Rscript dev\hypoxia_peak_reduced\02_run_model.R
Rscript dev\hypoxia_peak_reduced\03_render_outputs.R
Rscript dev\hypoxia_peak_reduced\04_build_karyotype_lookup.R
Rscript dev\hypoxia_peak_reduced\05_run_ablations.R
Rscript dev\hypoxia_peak_reduced\06_run_extensions.R
```

Optional posterior follow-up:

```powershell
Rscript dev\hypoxia_peak_reduced\07_run_nuts_extension.R . delta_dna 3 1000 1000 0.9 12
Rscript dev\hypoxia_peak_reduced\08_render_nuts_counterfactual.R . delta_dna

Rscript dev\hypoxia_peak_reduced\07_run_nuts_extension.R . adjacent_8n_delta_dna 3 1000 1000 0.9 12
Rscript dev\hypoxia_peak_reduced\08_render_nuts_counterfactual.R . adjacent_8n_delta_dna
```

You can optionally pass the optimizer iteration limit as the second argument:

```powershell
Rscript dev\hypoxia_peak_reduced\05_run_ablations.R . 8000
Rscript dev\hypoxia_peak_reduced\06_run_extensions.R . 8000
```

To fit using the gated-cell count ratio instead of the fluorescence-sum ratio:

```powershell
Rscript dev\hypoxia_peak_reduced\01_prepare_input.R . cell_count
Rscript dev\hypoxia_peak_reduced\02_run_model.R
Rscript dev\hypoxia_peak_reduced\03_render_outputs.R
Rscript dev\hypoxia_peak_reduced\04_build_karyotype_lookup.R
Rscript dev\hypoxia_peak_reduced\05_run_ablations.R
Rscript dev\hypoxia_peak_reduced\06_run_extensions.R
```

## Current ablations

The current ablation batch includes:

- `baseline`
- `no_date_effect`: forces `phi_date = 1`; this is the main negative control
  and performs much worse on the lower tumor peak
- `no_x_ratio_effect`: forces `beta = 0`; this currently changes the fit much
  less than removing the date effect
- `rho_fixed_one`: forces `rho = 1`; this is close to baseline and is not a
  major interpretive driver

## Current extensions

- `latent_u`: replaces deterministic `u_s = exp(-(alpha + beta x_s))` with a
  latent sample-level `log_u_s ~ N(alpha + beta x_s, sigma_u)`
- `delta_dna`: adds sample-specific DNA scaling while keeping the baseline
  state assignments
- `adjacent_8n`: fits the lower tumor peak with the original 2N/4N mixture,
  adds a global 8N/4N ratio, and constrains any retained second tumor peak to
  the adjacent state (4N if the lower peak is 2N, 8N if the lower peak is 4N)
- `adjacent_8n_delta_dna`: combines the adjacent-state interpretation with
  sample-specific DNA scaling and is the current most interesting working
  hypothesis

All outputs are written under `dev/hypoxia_peak_reduced/output`.
