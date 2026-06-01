# Portal API Output Requirements

This document lists the minimum output files and columns needed to populate the current portal API examples. The API wrapper fields such as `profile`, `index`, `q`, `count`, `page`, and `nonce` are response-envelope fields and do not need to be produced by the cell-state pipeline.

The examples inspected were:

- `cell-state-metadata`
- `cell-state-metadata-extended`
- `qc-metadata`
- `qc-metadata-extended`
- `pankbase-scb-expression`
- `pankbase-scb-gene-set-factor-cs`
- `pankbase-scb-gene-factor`
- `cell-state-expression`
- `pankbase-scb-factor-cs`
- `pankbase-scb-factor`
- `cell-state-pigean`
- `cell-state-heat-map`

## Required API-Ready Files

### 1. Curated Cell-State Metadata

Source output:

- `dat/api/cell_state_index.json`
- `dat/api/cell_state_details_by_id.json`
- `dat/api/curated_cell_state_manifest.tsv`
- `dat/api/curated_cell_state_markers.tsv`
- `dat/api/curated_cell_state_citations.tsv`

Used by:

- `cell-state-metadata`
- `cell-state-metadata-extended`

Minimal metadata columns for the compact endpoint:

- `state_id`
- `display_name`
- `tissue_id` or `tissue`
- `tissue_label`
- `cell_type_id` or `cell_type`
- `cell_type_label`
- `state_label`
- `state_class`
- `release_class`
- `interpretation_status`
- `manual_review_status`
- `n_markers`

Extended endpoint fields, from details JSON:

- `state_id`
- `display_name`
- `tissue`
- `tissue_label`
- `cell_type`
- `cell_type_label`
- `state`
- `summary`
- `marker_set`
- `state_level_citations`
- `curation`
- `scoring`
- `human_genetics`
- `quality`
- `related`

Important nested fields:

- `state.label`, `state.class`, `state.release_class`, `state.interpretation_status`, `state.is_composite_required`, `state.is_qc`, `state.allow_hard_call`, `state.score_scope`
- `summary.short_description`, `summary.curation_notes`, `summary.recommended_display`
- `marker_set.source_gmt`, `marker_set.source_workbook`, `marker_set.gene_set_description`, `marker_set.n_markers`, `marker_set.markers`
- marker rows need `gene`, `role`, `evidence_level`, `marker_notes`, `citations`, `source_type`, `from_excel`, `from_gmt`
- `curation.curation_version`, `curation.manual_review_status`, `curation.curated_by`, `curation.last_reviewed`, `curation.provenance_files`, `curation.provenance_warnings`

### 2. QC Metadata

Source output:

- `dat/api/qc_state_index.json`
- `dat/api/qc_state_details_by_id.json`

Used by:

- `qc-metadata`
- `qc-metadata-extended`

Compact endpoint columns:

- `qc_signature_id`
- `display_name`
- `category`
- `tier`
- `recommended_use`
- `dummy`

Extended endpoint columns:

- `qc_signature_id`
- `display_name`
- `category`
- `tier`
- `recommended_use`
- `source`
- `source_gmt`
- `exclude_when`
- `markers`
- `dummy`

`dummy` is a constant query helper column with value `1`.

### 3. Curated Cell-State Expression

Source output:

- `expression/curated_state_expression_specificity_cp10k.tsv.gz`

Fallback source if not split:

- `expression/all_gene_state_expression_specificity_cp10k.tsv.gz`, filtered to curated biological states

Used by:

- `cell-state-expression`

Required columns:

- `gene`
- `tissue`
- `cell_type`
- `state_name`
- `log10_cpk`
- `log2fc_weighted_vs_all_parent`
- `p_value`

Pipeline columns to retain or transform:

- `gene`
- `tissue`
- `annotated_cell_type` as `cell_type`
- `state_name`
- `weighted_mean_cp10k`
- `weighted_mean_expression`
- `log2fc_weighted_vs_all_parent`
- `p_value`
- `q_value`
- `state_weight_type`
- `expression_unit`

Recommended portal transform:

- Keep `state_weight_type == gradient_percentile_squared` as the default display weight.
- Compute `log10_cpk = log10(weighted_mean_cp10k + 1)` when `weighted_mean_cp10k` is available.
- If the run used normalized/log-normalized input, compute the same display field from `weighted_mean_expression` and keep `expression_unit` for provenance.

### 4. Program/Factor Expression

Source output:

- `expression/program_expression_specificity_cp10k.tsv.gz`

Used by:

- `pankbase-scb-expression`

Required columns:

- `dataset`
- `cell_type`
- `model`
- `factor`
- `gene`
- `log10_cpk`
- `log2fc_weighted_vs_all_parent`
- `p_value`

Pipeline columns to retain or transform:

- `gene`
- `annotated_cell_type` as `cell_type`
- `state_name` or program signature name as `factor`
- `weighted_mean_cp10k`
- `weighted_mean_expression`
- `log2fc_weighted_vs_all_parent`
- `p_value`
- `q_value`
- `state_weight_type`
- `expression_unit`
- constant or config fields `dataset` and `model`

Recommended portal transform:

- Convert program state names such as `pancreas_beta_liger_factor12` to display `Factor12`.
- Use `state_weight_type == gradient_percentile_squared` by default.
- Compute `log10_cpk = log10(weighted_mean_cp10k + 1)` or the expression-unit equivalent.

### 5. Program/Factor Metadata

Required output:

- `program_factor_metadata.tsv.gz`

Used by:

- `pankbase-scb-factor`
- `pankbase-scb-factor-cs`

Required columns for `pankbase-scb-factor`:

- `dataset`
- `cell_type`
- `model`
- `factor`
- `importance`
- `label`
- `top_cells`
- `top_gene_sets`
- `top_genes`
- `top_traits`

Additional columns for `pankbase-scb-factor-cs`:

- `significant_cell_states`
- `qc_cell_states`

Likely sources:

- LIGER or program loadings for `top_genes`
- Program cell activity/usage for `top_cells`
- PIGEAN factor results for `top_traits`
- Program-to-state matching summary for `significant_cell_states` and `qc_cell_states`
- Existing model metadata or downstream curation for `label` and `importance`

### 6. Program Gene Loadings

Required output:

- `program_gene_loadings.tsv.gz`

Used by:

- `pankbase-scb-gene-factor`

Required columns:

- `dataset`
- `cell_type`
- `model`
- `factor`
- `gene`
- `value`

Likely source:

- Per-cell-type program loading files such as `inputs/programs/<cell_type>.program_loadings.tsv.gz`.

Portal transform:

- Convert wide loadings to long format.
- Rename loading column to `value`.
- Convert program IDs to display factors such as `Factor1`.

### 7. Program/Factor PIGEAN Trait Associations

Required output:

- `program_pigean_trait_results.tsv.gz`

Used by:

- `pankbase-scb-gene-set-factor-cs`

Required columns:

- `dataset`
- `cell_type`
- `model`
- `factor`
- `trait`
- `beta`
- `beta_uncorrected`

Likely source:

- PIGEAN multi-y results run on program-derived GMTs, one run per cell type/signature method as appropriate.

### 8. Curated Cell-State PIGEAN Trait Associations

Required output:

- `cell_state_pigean_trait_results.tsv.gz`

Used by:

- `cell-state-pigean`

Required columns:

- `tissue`
- `cell_type`
- `state_name`
- `trait`
- `beta`
- `beta_uncorrected`

Likely source:

- PIGEAN multi-y results run on curated cell-state GMTs.

### 9. Program-to-State Heat Map

Required output:

- `program_state_heatmap.tsv.gz`

Used by:

- `cell-state-heat-map`

Required columns:

- `tissue`
- `cell_type`
- `state_name`
- `program_id`
- `correlation`
- `gsea_p`
- `gsea_q`

Likely source:

- `program_state_matches/<cell_type>/program_state_match_summary.tsv.gz`

Pipeline columns to retain or transform:

- `state_id` or `state_name`
- `program_id`
- `cell_spearman_r_gradient` as `correlation`
- `gsea_p`
- `gsea_q`
- `state_type`
- `match_class`
- `combined_match_score`
- configured `tissue` and `cell_type`

## Intermediate Files Still Needed to Build API Tables

These do not need to be served directly by the API but are needed to compute the API-ready tables above:

- `scoring/cell_state_activity.tsv.gz`
- `scoring/curated_state_activity.tsv.gz`, if split
- `scoring/program_activity_from_signature_scoring.tsv.gz`, if using program signatures as state-like activity
- `inputs/programs/<cell_type>.program_loadings.tsv.gz`
- `inputs/programs/<cell_type>.program_cell_activity.tsv.gz`
- `program_state_matches/<cell_type>/program_state_match_summary.tsv.gz`
- `program_state_matches/<cell_type>/program_qc_match_summary.tsv.gz`
- `pigean/<method>/gene_set_stats.out.gz` or combined PIGEAN result tables

## Minimum Pipeline Deliverable Set

For the listed APIs, the minimum final deliverable should be these API-ready files:

1. `dat/api/cell_state_index.json`
2. `dat/api/cell_state_details_by_id.json`
3. `dat/api/qc_state_index.json`
4. `dat/api/qc_state_details_by_id.json`
5. `portal/cell_state_expression.tsv.gz`
6. `portal/program_expression.tsv.gz`
7. `portal/program_factor_metadata.tsv.gz`
8. `portal/program_gene_loadings.tsv.gz`
9. `portal/program_pigean_trait_results.tsv.gz`
10. `portal/cell_state_pigean_trait_results.tsv.gz`
11. `portal/program_state_heatmap.tsv.gz`

The existing pipeline already produces most source material for files 1-4 and 5-6. Files 7-11 require portal-specific post-processing and/or PIGEAN runs from the existing program-state matching, program loading, and PIGEAN outputs.
