# Program-to-State Matching

`scripts/match_programs_to_cell_states.py` matches data-driven gene programs to
curated biological states and optional QC/artifact signatures for one
tissue/cell type at a time.

The workflow combines:

- Gene-level enrichment of curated state markers in program loadings.
- Optional cell-level correlation between program usage and cell-state activity.
- Optional gene-level comparison between program genes and state-specific
  all-gene expression scores from `cell_state_de`.

It uses only Python standard library, `pandas`, `numpy`, `scipy`, and `PyYAML`.

## Minimal Marker-Only Run

```bash
python scripts/match_programs_to_cell_states.py \
  --program-loadings results/program_loadings.tsv.gz \
  --state-gmt dat/pancreas/pancreas_beta_cell_states.gmt \
  --out-dir results/program_state_matching/beta
```

This writes marker-enrichment outputs and header-only optional outputs for
cell-correlation and state-expression matching.

## Full Run

```bash
python scripts/match_programs_to_cell_states.py \
  --program-loadings results/program_loadings.tsv.gz \
  --state-gmt dat/pancreas/pancreas_beta_cell_states.gmt \
  --qc-gmt dat/qc/cmdkp_all_tissues_minimal_bad_cell_qc_signatures.gmt \
  --program-cell-activity results/program_cell_activity.tsv.gz \
  --cell-state-activity results/cell_state_activity.tsv.gz \
  --state-expression results/all_gene_state_expression_specificity_cp10k.tsv.gz \
  --metadata results/metadata.tsv.gz \
  --donor-col donor_id \
  --tissue pancreas \
  --cell-type beta_cell \
  --out-dir results/program_state_matching/beta
```

QC signatures are matched with the same gene-level and cell-level methods as
biological states, but they are reported as QC/artifact matches and are not
interpreted as biological labels.

## Input Formats

Program loadings may be long:

```text
program_id    gene    loading
program_1     INS     1.2
program_1     IAPP    1.1
```

Aliases accepted for `program_id` are `program`, `factor`, `module`, and
`component`. Aliases accepted for `gene` are `gene_symbol` and `gene_name`.
Aliases accepted for `loading` are `weight`, `value`, and `loading_score`.

Program loadings may also be wide, where the first column is gene and all
remaining columns are program IDs.

Program cell activity may be long:

```text
cell_id    program_id    program_activity
cell_1     program_1     0.9
cell_2     program_1     0.1
```

Aliases accepted for `program_activity` are `usage`, `score`, `weight`, and
`activity`. Wide format is also accepted, where the first column is `cell_id`
and remaining columns are program IDs.

State activity should be the `cell_state_activity.tsv.gz` output from
`cell_state_de`, with `cell_id`, `state_name`, `state_type`,
`state_activity_weight_gradient`, `state_activity_weight_hightail`,
`aucell_score`, and `ucell_score`.

State expression should be
`all_gene_state_expression_specificity_cp10k.tsv.gz`.

## Marker Universe

Program loadings are often available only for selected variable or relevant
genes. The matcher therefore uses the program-specific loading universe for
gene-level tests.

Absent state marker genes are reported in `missing_state_markers`; they are not
treated as zero-loading genes. If fewer than `--min-marker-overlap` markers are
available or marker coverage is below `--min-marker-coverage`, the result is
flagged as `insufficient_marker_coverage`.

## Outputs

- `program_state_marker_enrichment.tsv.gz`: GSEA-like enrichment, AUROC/MWU,
  top-N overlaps, marker coverage, missing markers, and q-values.
- `program_state_expression_score_match.tsv.gz`: optional state-expression
  score matching for gradient and high-tail state weights.
- `program_state_cell_correlation.tsv.gz`: optional cell-level and donor-level
  program/state activity correlations.
- `program_state_match_summary.tsv.gz`: one row per program/state with combined
  evidence, match class, interpretation, and QC caveat.
- `program_state_heatmap_matrix.tsv.gz`: wide program x state combined-score
  matrix.
- `program_state_heatmap_long.tsv.gz`: long-form heatmap values and match
  classes.
- `program_qc_match_summary.tsv.gz`: best QC/artifact match per program when
  QC signatures are supplied.
- `program_label_suggestions.tsv.gz`: suggested program label and quality class.
- `run_summary.json`: input files, counts, parameters, versions, warnings, and
  timestamp.

## Interpretation

The combined score is intentionally transparent:

```text
gene_score = max(0, GSEA_NES) * -log10(GSEA_q)
cell_score = max(0, cell_spearman_r_gradient) * -log10(cell_spearman_q_gradient)
combined_match_score = gene_score + cell_score
```

The `-log10` terms are capped at 50. Scores are for ranking and review; the
`match_class` and evidence columns should be inspected before assigning a final
program label.
