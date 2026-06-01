# Interfaces

## Metadata TSV

Required columns:

- `cell_id`: stable cell identifier.
- `map_id` or equivalent configured map/dataset identifier.
- `tissue` or equivalent configured tissue column.
- `annotated_cell_type` or equivalent configured parent cell-type column.
- `donor_id` or equivalent configured donor column.

Conditionally required columns:

- `sample_id` or equivalent library/sample column: if absent in the general runner, donor ID is used as the sample grouping.
- Disease or case/control group column: required for differential expression.

Common optional columns:

- `cell_type`
- `disease_group`
- `sex`
- `age`
- `source`
- `chemistry`
- `treatment`

Column names are configurable in the command-line arguments, so maps do not need to use these exact names internally.

## Expression Long TSV

Required columns:

- `cell_id`
- `gene`
- `expression`

Expression should be non-negative normalized expression. Sparse long tables are allowed: missing cell-gene pairs are interpreted as zero by state assignment and differential expression scripts. For Seurat objects, the exporter uses the default assay `data` slot by default; pass `--layer` when extracting from a Seurat v5 layer.

For the general CMDKP runner, expression should be raw counts or log1p-normalized
expression. UMAP coordinates, PCA coordinates, integrated embeddings, and
batch-corrected latent spaces are invalid scoring inputs.

## State Spec JSON

Top-level format:

```json
{
  "cell_id_col": "cell_id",
  "donor_col": "donor_id",
  "states": [
    {
      "name": "mature_UCN3_MAFA_high",
      "type": "module_quantile",
      "genes": ["UCN3", "MAFA"],
      "direction": "high",
      "quantile": 0.75,
      "within": "donor"
    }
  ]
}
```

Supported state types:

- `module_quantile`: score each cell by mean expression of marker genes and threshold within all cells or within donor.
- `metadata_equals`: assign cells where a metadata column equals one value.
- `metadata_in`: assign cells where a metadata column matches any value from a list.

For `module_quantile`, `direction` must be `high` or `low`. If `within` is `donor`, the threshold is computed separately per donor and the metadata must include the configured donor column. If `within` is `all`, no donor column is required.

## State Membership TSV

Columns:

- `cell_id`
- `state`
- `in_state`
- `score`
- `threshold`
- `rule`

The table is long by state. A cell may belong to multiple states.

## Cell-State Score TSV

Long score table columns:

- `cell_id`
- `state`
- `score`
- `ucell_score`
- `score_method`
- `score_percentile_within_scope`
- `n_markers`
- `marker_genes_total`
- `n_markers_found`
- `marker_genes_present`
- `marker_coverage_fraction`
- `cell_type_scope`
- `library`

The wide score matrix has one row per cell, `cell_id` in the first column, and one score column per state.

For production GMT scoring, `score_method=ucell` is the default. The UCell-style
score is computed locally from per-cell gene ranks using a Mann-Whitney/rank
statistic and is normalized to the interval `[0, 1]`. Legacy
`score_method=mean_expression` is retained for comparability.

## State Threshold TSV

Columns:

- `tissue`
- `cell_type`
- `state`
- `threshold_method`
- `threshold_value`
- `null_percentile`
- `mixture_boundary`
- `n_cells_used`
- `marker_coverage_fraction`
- `null_n`
- `score_method`
- `max_rank`

The production default is `matched_random_gene_set_null99`: random gene sets of
the same size are matched by expression and detection bins, scored with the same
method, and the state threshold is the configured null percentile.

## State Calls TSV

Columns include:

- `cell_id`
- `state`
- `state_kind`
- `score_for_call`
- `threshold_value`
- `call`
- `confidence`
- `reason`
- `in_state`
- `call_score_column`

Supported `call` values are `active`, `inactive`, `ambiguous`,
`insufficient_coverage`, and `qc_flagged`. The compatibility column `in_state`
is true only when `call == active`.

## Multi-Label Cell Annotation TSV

Columns:

- `cell_id`
- `parent_cell_type`
- `active_biological_states`
- `active_process_flags`
- `qc_flags`
- `primary_interpretation`

`qc_*` signatures should be represented as flags rather than portal-facing
biological state labels.

## General CMDKP Runner Outputs

`run_cmdkp_state_scoring.py` writes the standard output bundle below.

The simplified continuous workflow additionally writes AUCell/UCell scores,
AUCell-based hard calls, state-excess soft weights, expression, DE, and
run-summary outputs. These are the preferred downstream interfaces.

The Python runner separates the broad scoring input from the expression-summary
input. Use `--rank-10x-dir` or `--rank-matrix-mtx --rank-genes --rank-cells` for
the full sparse rank universe used by AUCell/UCell scoring. Use
`--expression-matrix` for the small query-gene expression table used by expression
and DE summaries. The full rank universe should not be exported as long-form
text or dense-pivoted in pandas. Sparse rank-universe scoring is done directly
from the Matrix Market/10x matrix and scores all signatures in a GMT pass from
one per-cell sparse rank ordering. For Seurat inputs,
`export_rank_universe_10x_from_seurat.R` writes the required
`matrix.mtx.gz`, `features.tsv.gz`, and `barcodes.tsv.gz` files.

### `ucell_scores.tsv.gz`

Columns:

- `cell_id`
- `state_type`
- `state_name`
- `ucell_score`
- `aucell_score`
- `n_markers_present`
- `n_markers_total`
- `marker_coverage_fraction`
- `markers_present`
- `markers_missing`

### `aucell_state_activity.tsv.gz`

Minimal AUCell state activity output. It contains only
biological states:

- `cell_id`
- `state_name`
- `aucell_score`
- `threshold_status`
- `hard_call`
- `state_activity_weight`

### `cell_state_activity.tsv.gz`

Columns:

- `cell_id`
- `state_type`
- `state_name`
- `ucell_score`
- `aucell_score`
- `state_activity_weight`
- `soft_weight_method`
- `threshold_value`
- `threshold_status`
- `threshold_method`
- `q75_score`
- `q90_score_diagnostic`
- `q95_score_diagnostic`
- `q99_score`
- `marker_coverage_fraction`

Activity details:

- `ucell_score` is retained as a scalable signature score.
- `aucell_score` is used for biological hard-call thresholding and biological
  soft weights.
- `state_activity_weight` is threshold-to-q99 scaled for hard-callable states
  and q75-to-q99 scaled for continuous-only states.
- `q90_score_diagnostic` and `q95_score_diagnostic` are reported only as
  distribution diagnostics and are not used as default hard-call thresholds.
- Weights are independent across states and are not probabilities.

Expression and DE summary interfaces:

- `--query-genes` accepts a newline-delimited gene list and `--query-gene`
  accepts comma-separated genes. These options restrict
  `expression_expected_assignments.tsv.gz`,
  `expression_hard_assignments.tsv.gz`, `de_expected_assignments.tsv.gz`, and
  `de_hard_assignments.tsv.gz`.
- Query-gene restriction does not change state scoring or hard state calls.
- `--mode expected|hard|both` controls which expression and DE summary families
  are computed. Unrequested summary outputs are written as header-only files.

### `cell_state_hard_assignments.tsv.gz`

Columns:

- `cell_id`
- `state_type`
- `state_name`
- `aucell_score`
- `ucell_score`
- `state_activity_weight`
- `threshold`
- `hard_call`
- `threshold_source`
- `marker_coverage_pass`
- `reason`

### `qc_exclusions.tsv.gz`

Columns:

- `cell_id`
- `excluded`
- `exclusion_reason`
- `triggering_qc_states`
- `max_qc_activity`

### `expression_expected_assignments.tsv.gz`

Columns:

- `gene`
- `state_name`
- `state_type`
- `weighted_mean_expression`
- `weighted_detection_fraction`
- `whole_cell_type_mean_expression`
- `log2_weighted_vs_all`
- `n_cells`
- `sum_state_weight`
- `leave_one_gene_out_used`

### `expression_hard_assignments.tsv.gz`

Columns:

- `gene`
- `state_name`
- `state_type`
- `mean_expression_state_positive`
- `mean_expression_state_negative`
- `detection_fraction_positive`
- `detection_fraction_negative`
- `log2_fc_positive_vs_negative`
- `p_value`
- `q_value`
- `n_positive_cells`
- `n_negative_cells`

### `de_expected_assignments.tsv.gz` and `de_hard_assignments.tsv.gz`

These files are donor-level phenotype association summaries. They are written
with headers even when no phenotype table is supplied.

### `state_summary.tsv.gz` and `run_summary.json`

`state_summary.tsv.gz` provides one row per state with mean score/activity
and hard-call counts. `run_summary.json` records inputs, parameters, dimensions,
software versions, timestamp, and the number of QC-excluded cells.

### `cell_state_scores.tsv.gz`

One row per cell per relevant biological state:

- `cell_id`
- `map_id`
- `tissue`
- `annotated_cell_type`
- `state_name`
- `ucell_score`
- `aucell_score`
- `score_percentile_within_calibration_group`
- `n_markers_requested`
- `n_markers_present`
- `marker_coverage_fraction`
- `markers_present`
- `markers_missing`

### `cell_state_thresholds.tsv.gz`

One row per map/tissue/cell-type/state threshold eligibility group:

- `map_id`
- `tissue`
- `annotated_cell_type`
- `state_name`
- `threshold_method`
- `threshold_value`
- `mixture_threshold`
- `n_cells_in_calibration_group`
- `score_iqr`
- `q50_score`
- `q75_score`
- `q90_score_diagnostic`
- `q95_score_diagnostic`
- `q99_score`
- `active_fraction_at_threshold`
- `threshold_reason`
- `n_markers_requested`
- `n_markers_present`
- `marker_coverage_fraction`
- `threshold_status`

### `cell_state_calls.tsv.gz`

One row per cell per biological state:

- `cell_id`
- `map_id`
- `tissue`
- `annotated_cell_type`
- `state_name`
- `ucell_score`
- `threshold_value`
- `call`
- `confidence`
- `reason`
- `hard_exclusion_flag`
- `review_flag`
- `composite_rule_used`
- `requires_composite_validation`

### `qc_signature_scores.tsv.gz`

One row per cell per auxiliary QC signature:

- `cell_id`
- `map_id`
- `tissue`
- `annotated_cell_type`
- `qc_signature_name`
- `ucell_score`
- `aucell_score`
- `score_percentile_within_sample`
- `n_markers_present`
- `marker_coverage_fraction`

### `bad_cell_qc_flags.tsv.gz`

One row per cell with technical metrics and bad-cell flags:

- `cell_id`
- `map_id`
- `tissue`
- `annotated_cell_type`
- `total_counts`
- `n_genes_detected`
- `percent_mitochondrial`
- `percent_ribosomal`
- `percent_hemoglobin`
- `percent_malat1`
- `doublet_score`
- `hard_exclusion_flag`
- `review_flag`
- `bad_cell_reason`
- `top_offtarget_identity_signature`
- `top_ambient_signature`
- `parent_identity_score`

### `cell_multilabel_state_summary.tsv.gz`

One row per cell with active states, flags, recommendations, and a compact label.

### `state_call_summary.tsv.gz`

One row per map/tissue/cell-type/state with scored/active cell counts, donor
summaries, threshold details, and confidence summary.

### `state_scoring_method.md`

Narrative method record including matrix type, local UCell-style parameters, GMT
paths, marker coverage rules, thresholding, QC flag logic, exclusion policy, and
multi-label caveats.

## Differential Expression TSV

Columns:

- `analysis_type`
- `cell_type`
- `state`
- `gene`
- `log_fc`
- `pvalue`
- `qvalue`
- `case_donors`
- `control_donors`
- `case_cells`
- `control_cells`
- `case`
- `control`
- `test`

DE is donor-pseudobulk by default. The Seurat pseudobulk workflow sums counts by donor and group, then uses `edgeR` quasi-likelihood testing.

Supported `analysis_type` values:

- `cell_type`: disease contrast within each cell type.
- `state`: disease contrast within each assigned state.
- `state_association`: state-positive versus state-negative cells within the matching cell type.

## Gene-State Assignment TSV

Columns:

- `gene`
- `state`
- `library`
- `assignment_source`
- `is_marker`
- `log_fc`
- `pvalue`
- `qvalue`
- `assignment_pass`

`assignment_source` is `gmt_marker` for curated marker membership or `state_association_de` for data-driven state-positive associations.

## Matrix Value Types

Sparse rank and expression inputs can declare their value type explicitly. Supported values are `raw_counts`, `linear_cp10k`, `log1p_cp10k`, `linear_normalized`, `log1p_normalized`, and `auto`.

Use `--rank-value-type` with `run_cmdkp_state_scoring.py`. Rank-based AUCell/UCell scoring accepts non-negative raw, normalized, or log-normalized values because scoring depends on within-cell gene ranks. Negative scaled values are rejected by default.

Use `--expression-value-type` with `summarize_state_expression.py`. Raw counts are converted to CP10K using cell totals. Linear CP10K-like values are averaged directly. Log-normalized CP10K-like values are first transformed with `expm1(value)` and then averaged. Generic normalized values are handled the same way but are labeled as normalized expression rather than strict CP10K.

When `auto` is used, the code first looks for `matrix_value_type.json` in the 10x-like directory and then falls back to distribution-based inference. Low-confidence inference fails and requires an explicit value type.

For normalized inputs, direct QC metrics should come from metadata fields such as `QC:nCount_RNA`, `QC:nFeature_RNA`, and `QC:percent.mt`. QC signature scores can still be computed from the normalized rank matrix.
