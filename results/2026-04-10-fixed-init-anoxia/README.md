# Fixed-init anoxia optimization checkpoint

This archive captures the fixed-initialization CmdStan optimization result for
the anoxia-only subset from April 10, 2026.

The authoritative code state is the git commit that contains this directory.
The parent commit before the checkpoint was:

```text
744ed5857bfc9139d1e7cbdd017fcdcb319632ce
```

## Result

- `fixed_init_optimize-1.csv`: CmdStan optimize output for the anoxia-only fixed-init run.

The plots, RDS manifest, compiled CmdStan executable, gating cache, and full
`processed_data/` tree are intentionally excluded from this archive.

## Reproduction Path

1. Raw flow cytometry files under `data/Anoxia_FlowCytometry/`
2. Gated DNA-area export at `processed_data/anoxia-flowcytometry/filtered_dna_area_vectors.rds`
3. Stan data bundle at `processed_data/stan_data.Rds`
4. Stan model at `stan/ploidy_histogram_first_pass.stan`
5. Driver script at `R/fixed_init_optimize_and_plot.R`
6. Optimizer output at `processed_data/stan_optimize_fixed_init_anoxia/fixed_init_optimize-1.csv`

Run settings:

- `experiment_filter`: `anoxia`
- CmdStan method: `optimize`
- CmdStan algorithm: `lbfgs`
- CmdStan version: `2.37.0`
- model name in output header: `ploidy_histogram_first_pass_model`

## Hash Manifests

- `hashes/raw_fcs_sha256.csv`: SHA-256 hashes for the anoxia raw `.fcs` files.
- `hashes/derived_and_result_sha256.csv`: SHA-256 hashes for the gated DNA export, Stan data bundle, and optimizer CSV.
- `hashes/source_sha256.csv`: SHA-256 hashes for the R Markdown workflow, Stan data prep script, fixed-init optimization script, and Stan model.
