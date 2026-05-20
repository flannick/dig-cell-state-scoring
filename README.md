# `cell_state_de`

Reusable scripts for assigning cells to marker-defined states and running donor-aware differential expression on single-cell maps.

The toolkit is intentionally map-agnostic. PanKbase is the first target, but the interfaces are plain TSV plus JSON so the same scripts can be used after exporting selected genes and metadata from any Seurat or AnnData-style map.

## Selected-Gene Workflow

Run these commands from the analysis project root, not from inside this directory.

1. Export selected genes and metadata from a single-cell object.

```bash
R_LIBS_USER=../.Rlib /opt/homebrew/bin/Rscript --vanilla cell_state_de/scripts/extract_selected_expression_from_seurat.R \
  --rds data/external/pankbase/060425_scRNA_v3.3.rds \
  --genes cell_state_de/configs/example_genes.txt \
  --metadata-out cell_state_de/results/metadata.tsv.gz \
  --expression-out cell_state_de/results/expression_long.tsv.gz \
  --layer data
```

2. Assign marker-defined states.

```bash
../.venv/bin/python cell_state_de/scripts/assign_cell_states.py \
  --metadata cell_state_de/results/metadata.tsv.gz \
  --expression cell_state_de/results/expression_long.tsv.gz \
  --state-spec cell_state_de/configs/example_state_spec.json \
  --out cell_state_de/results/cell_state_membership.tsv.gz
```

3. Run donor-aware differential expression.

```bash
../.venv/bin/python cell_state_de/scripts/donor_pseudobulk_de.py \
  --metadata cell_state_de/results/metadata.tsv.gz \
  --expression cell_state_de/results/expression_long.tsv.gz \
  --states cell_state_de/results/cell_state_membership.tsv.gz \
  --genes cell_state_de/configs/example_genes.txt \
  --state mature_UCN3_MAFA_high \
  --donor-col donor_id \
  --group-col disease_group \
  --case T2D \
  --control ND \
  --out cell_state_de/results/mature_UCN3_MAFA_high_T2D_vs_ND.tsv
```

## Interfaces

- Metadata TSV: one row per cell. Must include a cell ID column and donor/group columns needed for DE.
- Expression long TSV: one row per observed cell-gene pair with `cell_id`, `gene`, and `expression`; omitted pairs are interpreted as zero by the Python scripts.
- State spec JSON: state names and rules for marker-module quantiles or metadata predicates.
- State membership TSV: one row per cell-state pair with boolean membership and the score used for assignment when applicable.
- DE TSV: one row per gene with donor counts, mean expression, log2 fold change, p-value, and BH q-value.

See [docs/interfaces.md](docs/interfaces.md) for the exact schema.

## GMT State Workflow

Use this workflow when state definitions come from a GMT marker file such as `../cell_states/out/pancreas/pancreas_cell_state_markers.gmt`.

For PanKbase, first stage the Seurat RDS from the tarball and remove it after derived outputs are complete.

```bash
mkdir -p results/cell_state_de/pankbase_beta/tmp_rds
tar -xzf data/external/pankbase/pankbase-scrna-umap-v3.3.tar.gz \
  -C results/cell_state_de/pankbase_beta/tmp_rds \
  --strip-components 1 \
  pankbase-scrna-umap-v3.3/060425_scRNA_v3.3.rds
```

1. Score cells for marker-defined states. The production default is local UCell-style
   rank scoring. The script can also write calibrated matched-random null
   thresholds for state calling.

```bash
R_LIBS_USER=../.Rlib /opt/homebrew/bin/Rscript --vanilla cell_state_de/scripts/score_gmt_states_from_seurat.R \
  --rds results/cell_state_de/pankbase_beta/tmp_rds/060425_scRNA_v3.3.rds \
  --gmt ../cell_states/out/pancreas/pancreas_cell_state_markers.gmt \
  --state-regex '^pancreas_beta_cell_' \
  --cell-filter-col Cell_Type \
  --cell-filter-values Beta \
  --metadata-cols 'Cell_Type,center_donor_id,description_of_diabetes_status,treatments' \
  --score-method ucell \
  --thresholds-out results/cell_state_de/pankbase_beta_state_thresholds.tsv.gz \
  --null-n 500 \
  --null-percentile 0.99 \
  --null-max-cells 20000 \
  --scores-out results/cell_state_de/pankbase_beta_state_scores.tsv.gz \
  --wide-out results/cell_state_de/pankbase_beta_state_scores_wide.tsv.gz \
  --metadata-out results/cell_state_de/pankbase_beta_state_metadata.tsv.gz
```

2. Call multi-label states from scores and calibrated thresholds.

```bash
../.venv/bin/python cell_state_de/scripts/call_states_from_scores.py \
  --scores results/cell_state_de/pankbase_beta_state_scores.tsv.gz \
  --thresholds results/cell_state_de/pankbase_beta_state_thresholds.tsv.gz \
  --metadata results/cell_state_de/pankbase_beta_state_metadata.tsv.gz \
  --parent-cell-type-col Cell_Type \
  --rules cell_state_de/configs/example_state_call_rules.json \
  --out results/cell_state_de/pankbase_beta_state_calls.tsv.gz \
  --annotation-out results/cell_state_de/pankbase_beta_cell_annotations.tsv.gz
```

Legacy exploratory quantile assignment remains available, but should not be used
as the production state-active definition because it forces the same approximate
fraction of cells into every state.

```bash
../.venv/bin/python cell_state_de/scripts/assign_states_from_scores.py \
  --scores results/cell_state_de/pankbase_beta_state_scores.tsv.gz \
  --metadata results/cell_state_de/pankbase_beta_state_metadata.tsv.gz \
  --cell-type-col Cell_Type \
  --method quantile \
  --within cell_type \
  --quantile 0.75 \
  --out results/cell_state_de/pankbase_beta_state_membership.tsv.gz
```

3. Run donor-pseudobulk differential expression from Seurat counts.

```bash
R_LIBS_USER=../.Rlib /opt/homebrew/bin/Rscript --vanilla cell_state_de/scripts/pseudobulk_de_from_seurat.R \
  --rds results/cell_state_de/pankbase_beta/tmp_rds/060425_scRNA_v3.3.rds \
  --membership results/cell_state_de/pankbase_beta_state_calls.tsv.gz \
  --analysis-types cell_type,state,state_association \
  --donor-col center_donor_id \
  --group-col description_of_diabetes_status \
  --cell-type-col Cell_Type \
  --cell-filter-col Cell_Type \
  --cell-filter-values Beta \
  --treatment-col treatments \
  --treatment-values no_treatment \
  --case-values 'type 2 diabetes' \
  --control-values non-diabetic \
  --out results/cell_state_de/pankbase_beta_state_de.tsv.gz
```

4. Assign genes to states from curated markers plus state association DE.

```bash
../.venv/bin/python cell_state_de/scripts/assign_genes_to_states.py \
  --gmt ../cell_states/out/pancreas/pancreas_cell_state_markers.gmt \
  --state-regex '^pancreas_beta_cell_' \
  --state-association-de results/cell_state_de/pankbase_beta_state_de.tsv.gz \
  --out results/cell_state_de/pankbase_beta_gene_state_assignments.tsv.gz
```

## Dependencies

Python scripts use `pandas`, `numpy`, and `scipy`. The R exporter requires `Seurat`.
GMT scoring requires `Seurat` and `Matrix`; pseudobulk DE requires `Seurat`, `Matrix`, and `edgeR`.
UCell-style scoring is implemented locally and does not require the R `UCell` package.

## General CMDKP State-Scoring Workflow

Use the general runner when applying CMDKP signatures to any exported map,
tissue, or annotated cell type. The input expression table must be a long TSV
with `cell_id`, `gene`, and `expression`. For Seurat maps, first export the
desired raw-count or log1p-normalized expression layer using the exporter above
or an equivalent map-specific export step. Do not use UMAP coordinates,
integrated embeddings, PCA coordinates, or other latent spaces as scoring input.

```bash
../.venv/bin/python cell_state_de/scripts/run_cmdkp_state_scoring.py \
  --expression-matrix results/cell_state_de/example_expression_long.tsv.gz \
  --cell-metadata results/cell_state_de/example_metadata.tsv.gz \
  --states-gmt out/pancreas/pancreas_cell_state_markers.gmt \
  --qc-states-gmt out/qc/cmdkp_all_tissues_minimal_bad_cell_qc_signatures.gmt \
  --out-dir results/cell_state_de/cmdkp_state_scoring \
  --expression-kind log1p_normalized \
  --map-id-col map_id \
  --tissue-col tissue \
  --cell-type-col cell_type \
  --donor-col donor_id \
  --sample-col sample_id \
  --query-genes results/cell_state_de/query_genes.txt \
  --mode both
```

Optional Ensembl-to-HGNC mapping is supported with:

```bash
--gene-map data/gene_id_to_hgnc.tsv
```

where the mapping table contains `gene_id` and `gene_symbol` columns.

The runner writes:

- `ucell_scores.tsv.gz`
- `cell_state_activity.tsv.gz`
- `cell_state_hard_assignments.tsv.gz`
- `qc_exclusions.tsv.gz`
- `expression_expected_assignments.tsv.gz`
- `expression_hard_assignments.tsv.gz`
- `de_expected_assignments.tsv.gz`
- `de_hard_assignments.tsv.gz`
- `state_summary.tsv.gz`
- `run_summary.json`
- `cell_state_scores.tsv.gz`
- `cell_state_thresholds.tsv.gz`
- `cell_state_calls.tsv.gz`
- `qc_signature_scores.tsv.gz`
- `bad_cell_qc_flags.tsv.gz`
- `cell_multilabel_state_summary.tsv.gz`
- `state_call_summary.tsv.gz`
- `state_scoring_method.md`

### Interpreting scores, activity, and hard calls

The primary score is a local UCell-style rank score:

```text
u_is = UCell(x_i, G_s)
```

where `x_i` is the expression vector for cell `i` and `G_s` is the marker set
for state `s`. The raw UCell score is the state activity weight used everywhere
downstream:

```text
a_is = u_is
```

There is no probability calibration, matched random gene-set null, local-FDR
step, or requirement that activities sum to one. Interpret `a_is` as a
continuous signature activity score, not as a posterior probability.

Expression and DE summaries use all expression-matrix genes by default. To
restrict those summaries to a query set, pass a newline-delimited list with
`--query-genes` or a comma-separated list with `--query-gene`. UCell scoring
still uses the full expression matrix and full state/QC GMTs; the query gene
options only restrict `expression_*` and `de_*` outputs. The
`--mode expected|hard|both` option controls whether expected/activity-weighted
summaries, hard-assignment summaries, or both are computed. Unrequested summary
files are still written as header-only tables for interface stability.

Hard calls are optional thresholded derivatives:

```text
I_is = 1[a_is >= tau_s]
```

Defaults are `tau_s = 0.80` for biological states and `0.95` for QC states,
with per-state overrides from YAML. Marker coverage must pass
`n_markers_present >= 5` and `marker_coverage_fraction >= 0.50` for hard calls.

QC exclusion is never applied silently. Cells are excluded only if
`--exclude-qc-above` is supplied, and `qc_exclusions.tsv.gz` records the
triggering QC states and reasons.

Expected expression uses raw UCell activity weights:

```text
E[g | s] = sum_i a_is x_ig / sum_i a_is
E_d[g | s] = sum_{i in donor d} a_is x_ig / sum_{i in donor d} a_is
```

When a query gene is part of a state marker set, leave-one-gene-out scoring is
used by default for expression summaries so the gene's own expression does not
directly define its state weight.

Phenotype summaries use simple donor-level linear models:

```text
E_d[g | s] ~ phenotype_d + covariates_d
```

If no phenotype table is supplied, DE output files are still written with
headers and DE is recorded as skipped in `run_summary.json`.
