# `cell_state_de`

Reusable scripts for assigning cells to marker-defined states and running donor-aware differential expression on single-cell maps.

The toolkit is intentionally map-agnostic. The interfaces are plain TSV plus JSON so the same scripts can be used after exporting selected genes and metadata from any Seurat or AnnData-style map.

## Selected-Gene Workflow

Run these commands from the analysis project root, not from inside this directory.

1. Export selected genes and metadata from a single-cell object.

```bash
R_LIBS_USER=../.Rlib /opt/homebrew/bin/Rscript --vanilla cell_state_de/scripts/extract_selected_expression_from_seurat.R \
  --rds data/external/example_map/example.rds \
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

Use this workflow when state definitions come from a GMT marker file.

1. Score cells for marker-defined states. The production default is local UCell-style
   rank scoring. The script can also write calibrated matched-random null
   thresholds for state calling.

```bash
R_LIBS_USER=../.Rlib /opt/homebrew/bin/Rscript --vanilla cell_state_de/scripts/score_gmt_states_from_seurat.R \
  --rds data/external/example_map/example.rds \
  --gmt results/cell_state_de/example_cell_state_markers.gmt \
  --state-regex '^tissue_a_cell_type_a_' \
  --cell-filter-col cell_type \
  --cell-filter-values 'Type A' \
  --metadata-cols 'cell_type,donor_id,condition,treatment' \
  --score-method ucell \
  --thresholds-out results/cell_state_de/example_state_thresholds.tsv.gz \
  --null-n 500 \
  --null-percentile 0.99 \
  --null-max-cells 20000 \
  --scores-out results/cell_state_de/example_state_scores.tsv.gz \
  --wide-out results/cell_state_de/example_state_scores_wide.tsv.gz \
  --metadata-out results/cell_state_de/example_state_metadata.tsv.gz
```

2. Call multi-label states from scores and calibrated thresholds.

```bash
../.venv/bin/python cell_state_de/scripts/call_states_from_scores.py \
  --scores results/cell_state_de/example_state_scores.tsv.gz \
  --thresholds results/cell_state_de/example_state_thresholds.tsv.gz \
  --metadata results/cell_state_de/example_state_metadata.tsv.gz \
  --parent-cell-type-col cell_type \
  --rules cell_state_de/configs/example_state_call_rules.json \
  --out results/cell_state_de/example_state_calls.tsv.gz \
  --annotation-out results/cell_state_de/example_cell_annotations.tsv.gz
```

Legacy exploratory quantile assignment remains available, but should not be used
as the production state-active definition because it forces the same approximate
fraction of cells into every state.

```bash
../.venv/bin/python cell_state_de/scripts/assign_states_from_scores.py \
  --scores results/cell_state_de/example_state_scores.tsv.gz \
  --metadata results/cell_state_de/example_state_metadata.tsv.gz \
  --cell-type-col cell_type \
  --state-cell-type-map results/cell_state_de/example_state_cell_type_map.tsv \
  --method quantile \
  --within cell_type \
  --quantile 0.75 \
  --out results/cell_state_de/example_state_membership.tsv.gz
```

3. Run donor-pseudobulk differential expression from Seurat counts.

```bash
R_LIBS_USER=../.Rlib /opt/homebrew/bin/Rscript --vanilla cell_state_de/scripts/pseudobulk_de_from_seurat.R \
  --rds data/external/example_map/example.rds \
  --membership results/cell_state_de/example_state_calls.tsv.gz \
  --analysis-types cell_type,state,state_association \
  --donor-col donor_id \
  --group-col condition \
  --cell-type-col cell_type \
  --cell-filter-col cell_type \
  --cell-filter-values 'Type A' \
  --case-values case \
  --control-values control \
  --out results/cell_state_de/example_state_de.tsv.gz
```

4. Assign genes to states from curated markers plus state association DE.

```bash
../.venv/bin/python cell_state_de/scripts/assign_genes_to_states.py \
  --gmt results/cell_state_de/example_cell_state_markers.gmt \
  --state-regex '^tissue_a_cell_type_a_' \
  --state-association-de results/cell_state_de/example_state_de.tsv.gz \
  --out results/cell_state_de/example_gene_state_assignments.tsv.gz
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
  --rank-10x-dir results/cell_state_de/example_rank_universe_10x \
  --expression-matrix results/cell_state_de/query_gene_expression_long.tsv.gz \
  --cell-metadata results/cell_state_de/example_metadata.tsv.gz \
  --states-gmt results/cell_state_de/example_cell_state_markers.gmt \
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
- `aucell_state_activity.tsv.gz`
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

Use separate inputs for scoring and expression summaries:

- AUCell/UCell state scoring should use a full sparse rank universe supplied
  with `--rank-10x-dir` or `--rank-matrix-mtx --rank-genes --rank-cells`.
- Expression and DE summaries should use a small query-gene expression matrix
  supplied with `--expression-matrix`.
- Do not export the full cell x gene matrix as long-form text for scoring, and
  do not use the small query-gene matrix as the AUCell rank universe for
  production runs.
- Sparse rank-universe scoring is performed directly from the Matrix Market/10x
  matrix and scores all signatures in a GMT pass from one per-cell sparse rank
  ordering. The broad matrix is not dense-pivoted in pandas.

For Seurat objects, export the rank universe with:

```bash
R_LIBS_USER=../.Rlib /opt/homebrew/bin/Rscript --vanilla cell_state_de/scripts/export_rank_universe_10x_from_seurat.R \
  --rds data/external/example_map/example.rds \
  --out-dir results/cell_state_de/example_rank_universe_10x \
  --counts-layer counts \
  --cell-filter-col cell_type \
  --cell-filter-values 'Type A' \
  --min-detection 0.01
```

The runner computes two local rank-based scores:

```text
u_is = UCell(x_i, G_s)
a_is = AUCell(x_i, G_s)
```

where `x_i` is the expression vector for cell `i` and `G_s` is the marker set
for state `s`. UCell is retained as a scalable signature score. AUCell is used
for hard state calls and soft state-excess weights.

The minimal AUCell state activity interface is
`aucell_state_activity.tsv.gz`, with `cell_id`, `state_name`, `aucell_score`,
`threshold_status`, `hard_call`, and `state_activity_weight`.

For each biological state, the runner explores the AUCell score distribution
within each `map_id + tissue + annotated_cell_type + state_name` group. A hard
threshold is accepted only when a minimum-density valley between lower- and
higher-AUCell regions supports a state-active population, or when an explicit
YAML threshold is supplied. Otherwise the state is marked `continuous_only` and
no hard-positive cells are forced.

Soft weights are threshold-centered when a hard threshold exists:

```text
w_is = clip((a_is - T_s) / (Q99_s - T_s), 0, 1)
```

For continuous-only states, soft weights use the high tail:

```text
w_is = clip((a_is - Q75_s) / (Q99_s - Q75_s), 0, 1)
```

`Q90_s` and `Q95_s` are reported only as score-distribution diagnostics. They
are not used as default hard-call thresholds because that would force fixed
fractions of cells into each state.

Weights are normalized within state across cells, not across states within a
cell. They do not sum to one and are not probabilities.

Expression and DE summaries use all expression-matrix genes by default. To
restrict those summaries to a query set, pass a newline-delimited list with
`--query-genes` or a comma-separated list with `--query-gene`. State scoring
uses the full sparse rank universe and full state/QC GMTs; the query gene
options only restrict `expression_*` and `de_*` outputs. The
`--mode expected|hard|both` option controls whether expected/activity-weighted
summaries, hard-assignment summaries, or both are computed. Unrequested summary
files are still written as header-only tables for interface stability.

Hard calls are optional thresholded derivatives:

```text
I_is = 1[a_is >= T_s]
```

Marker coverage must pass `n_markers_present >= 5` and
`marker_coverage_fraction >= 0.50` for hard calls. QC signatures use the QC
activity threshold, default `0.95`, unless overridden by YAML.

QC exclusion is never applied silently. Cells are excluded only if
`--exclude-qc-above` is supplied, and `qc_exclusions.tsv.gz` records the
triggering QC states and reasons.

Expected expression uses state-excess activity weights:

```text
E[g | s] = sum_i w_is x_ig / sum_i w_is
E_d[g | s] = sum_{i in donor d} w_is x_ig / sum_{i in donor d} w_is
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

## All-Gene State Expression

After scoring states, summarize every gene against the continuous state
activity weights using sparse raw counts:

```bash
../.venv/bin/python cell_state_de/scripts/summarize_state_expression.py \
  --raw-10x-dir results/cell_state_de/example_rank_universe_10x \
  --metadata results/cell_state_de/example_metadata.tsv.gz \
  --cell-state-activity results/cell_state_de/cmdkp_state_scoring/cell_state_activity.tsv.gz \
  --states-gmt results/cell_state_de/example_cell_state_markers.gmt \
  --out-dir results/cell_state_de/state_expression
```

The script computes CP10K from total raw counts per cell, then writes:

- `all_gene_all_parent_cp10k.tsv.gz`
- `all_gene_state_expression_cp10k.tsv.gz`
- `all_gene_state_specificity_cp10k.tsv.gz`
- `all_gene_state_expression_specificity_cp10k.tsv.gz`
- `state_expression_summary.json`

The primary expression unit is CP10K plus `log1p(mean CP10K)`. Specificity is
reported as a Spearman association between `log1p(CP10K)` and the state weight,
with BH FDR columns. Rows from low-signal or composite-required states are kept
and flagged rather than silently removed.

Two non-exclusive state weights are used:

- `gradient_percentile_squared`: within-state AUCell percentile squared.
- `high_tail_percentile_90_100`: the top AUCell percentile tail, scaled from
  percentile 0.90 to 1.00.

## State-Derived GMTs

Build separate GMTs for curated markers and state-derived expression signatures:

```bash
../.venv/bin/python cell_state_de/scripts/make_state_expression_gmts.py \
  --state-expression-specificity results/cell_state_de/state_expression/all_gene_state_expression_specificity_cp10k.tsv.gz \
  --original-state-gmt results/cell_state_de/example_cell_state_markers.gmt \
  --out-dir results/cell_state_de/state_gmts \
  --top-n 250
```

Outputs:

- `gmt/original_markers.gmt`
- `gmt/top_absolute_expression.gmt`
- `gmt/top_specific_fc.gmt`
- `gmt/top_specific_logp.gmt`
- `gmt_membership.tsv.gz`
- `gmt_build_summary.tsv`

The signature methods are deliberately separated so downstream PIGEAN runs do
not mix curated markers with expanded expression-derived signatures.

## PIGEAN Multi-Y

Stage or run PIGEAN separately for each GMT method:

```bash
../.venv/bin/python cell_state_de/scripts/run_pigean_multi_y_for_state_gmts.py \
  --pigean-command pigean \
  --multi-y-input results/example_multi_y.tsv.gz \
  --gmt-dir results/cell_state_de/state_gmts/gmt \
  --out-dir results/cell_state_de/pigean_state_gmts
```

For each method, the wrapper writes `run_command.txt`, `run_summary.json`, and
either PIGEAN outputs or a runnable `run_pigean.sh` if PIGEAN is unavailable or
`--dry-run` is used.

## Donor Phenotype Regression

Run donor-level models from sparse counts and cell-state activity:

```bash
../.venv/bin/python cell_state_de/scripts/run_state_phenotype_regression.py \
  --raw-10x-dir results/cell_state_de/example_rank_universe_10x \
  --metadata results/cell_state_de/example_metadata.tsv.gz \
  --cell-state-activity results/cell_state_de/cmdkp_state_scoring/cell_state_activity.tsv.gz \
  --phenotypes results/cell_state_de/example_phenotypes.tsv \
  --out-dir results/cell_state_de/state_de
```

Outputs:

- `de_whole_parent.tsv.gz`
- `de_state_weighted_gradient.tsv.gz`
- `de_state_weighted_hightail.tsv.gz`
- `de_state_specific_contrast.tsv.gz`
- `de_run_summary.json`

The model response is `log1p` donor-level mean CP10K. Numeric phenotypes are
standardized; categorical phenotypes use the first sorted level as the reference
group. Output tables include coefficient units and global, trait-level, and
gene-level FDR columns.

## Methods

CP10K:

```text
CP10K_ig = 10000 * c_ig / sum_g c_ig
```

State-weighted mean:

```text
mean_gs = sum_i w_is CP10K_ig / sum_i w_is
```

Detection:

```text
pct_gs = sum_i w_is I(c_ig > 0) / sum_i w_is
```

Gradient weight:

```text
r_is = percentile_rank(AUCell_is within calibration group)
w_gradient_is = r_is^2
```

High-tail weight:

```text
w_tail_is = clip((r_is - 0.90) / 0.10, 0, 1)
```

Donor state-weighted expression:

```text
y_dgs = sum_{i in donor d} w_is CP10K_ig / sum_{i in donor d} w_is
```

Regression:

```text
log1p(y_dgs) ~ phenotype_d + covariates_d
```
