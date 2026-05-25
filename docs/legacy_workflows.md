# Legacy Workflows

These commands are retained for compatibility and small exploratory analyses.
They are not the recommended production workflow for portal-scale state scoring.
Use the production workflow in the main README for multi-tissue or multi-cell-type runs.

## Selected-Gene Workflow

Run these commands from the analysis project root, not from inside this directory.

```bash
R_LIBS_USER=../.Rlib /opt/homebrew/bin/Rscript --vanilla cell_state_de/scripts/extract_selected_expression_from_seurat.R \
  --rds data/external/example_map/example.rds \
  --genes cell_state_de/configs/examples/example_genes.txt \
  --metadata-out cell_state_de/results/metadata.tsv.gz \
  --expression-out cell_state_de/results/expression_long.tsv.gz \
  --layer data

../.venv/bin/python cell_state_de/scripts/assign_cell_states.py \
  --metadata cell_state_de/results/metadata.tsv.gz \
  --expression cell_state_de/results/expression_long.tsv.gz \
  --state-spec cell_state_de/configs/examples/example_state_spec.json \
  --out cell_state_de/results/cell_state_membership.tsv.gz

../.venv/bin/python cell_state_de/scripts/donor_pseudobulk_de.py \
  --metadata cell_state_de/results/metadata.tsv.gz \
  --expression cell_state_de/results/expression_long.tsv.gz \
  --states cell_state_de/results/cell_state_membership.tsv.gz \
  --genes cell_state_de/configs/examples/example_genes.txt \
  --state example_marker_high \
  --donor-col donor_id \
  --group-col disease_group \
  --case case \
  --control control \
  --out cell_state_de/results/example_marker_high_case_vs_control.tsv
```

## Older GMT Score/Call Workflow

This older workflow scores GMTs from a Seurat object and then calls states from
precomputed score thresholds. The current production workflow scores sparse
rank-universe inputs directly with `run_cmdkp_state_scoring.py`.

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

../.venv/bin/python cell_state_de/scripts/call_states_from_scores.py \
  --scores results/cell_state_de/example_state_scores.tsv.gz \
  --thresholds results/cell_state_de/example_state_thresholds.tsv.gz \
  --metadata results/cell_state_de/example_state_metadata.tsv.gz \
  --parent-cell-type-col cell_type \
  --rules cell_state_de/configs/examples/example_state_call_rules.json \
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

## Legacy Pseudobulk and Gene Assignment

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

../.venv/bin/python cell_state_de/scripts/assign_genes_to_states.py \
  --gmt results/cell_state_de/example_cell_state_markers.gmt \
  --state-regex '^tissue_a_cell_type_a_' \
  --state-association-de results/cell_state_de/example_state_de.tsv.gz \
  --out results/cell_state_de/example_gene_state_assignments.tsv.gz
```
