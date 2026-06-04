# PIGEAN Full-Universe Default Run With Trait Blacklist

This is the PIGEAN command style to use for the production cell-state pipeline. The default blacklist excludes traits beginning with `HP_` or `exomes_`, and traits containing `gcat_` or `Orphanet`. It mirrors the validated run script:

`../blanc_screen/results/pigean_all_cell_states_full_universe_default_no_hpo_exomes_gcat_orphanet/run_pigean_all_cell_states_full_universe_default_no_hpo_exomes_gcat_orphanet.sh`

Run PIGEAN once per signature method within each cell type/group. Do not
concatenate GMTs across cell types before running PIGEAN.

## Inputs

```bash
PYTHON_BIN="../.venv/bin/python"
PIGEAN_SRC="/Users/flannick/codex-workspace/analysis/pigean_optimize/pigean/src"
MULTI_Y_IN="../resources/pigean/data/large/all.gene_stats.large.gt1.out.gz"
GENE_UNIVERSE="../resources/pigean/data/reference/NCBI37.3.plink.gene.loc"
TRAIT_BLACKLIST="results/pigean_all_cell_states_full_universe_default_no_hpo_exomes_gcat_orphanet/trait_blacklist_hp_exomes_gcat_orphanet.txt"
GMT_ROOT="results/cell_state_de/production_run/groups/<group_label>/state_expression_gmts/gmt"
OUT_ROOT="results/cell_state_de/production_run/groups/<group_label>/pigean"
```

`GMT_ROOT` should contain one GMT per signature method:

- `original_markers.gmt`
- `top_absolute_expression.gmt`
- `top_specific_fc.gmt`
- `top_specific_logp.gmt`

## Command

```bash
run_pigean() {
  local method="$1"
  local gmt="${GMT_ROOT}/${method}.gmt"
  local out_dir="${OUT_ROOT}/${method}"

  mkdir -p "${out_dir}"

  env PYTHONPATH="${PIGEAN_SRC}" \
  "${PYTHON_BIN}" -m pigean betas \
    --X-in "${gmt}" \
    --multi-y-in "${MULTI_Y_IN}" \
    --multi-y-id-col Gene \
    --multi-y-pheno-col Trait_Internal \
    --multi-y-log-bf-col Direct \
    --multi-y-combined-col Combined \
    --multi-y-prior-col Indirect \
    --multi-y-trait-blacklist-in "${TRAIT_BLACKLIST}" \
    --gene-universe-in "${GENE_UNIVERSE}" \
    --gene-universe-id-col 6 \
    --gene-universe-no-header \
    --gene-set-stats-out "${out_dir}/gene_set_stats.debug.out.gz" \
    --params-out "${out_dir}/params.out.gz" \
    --log-file "${out_dir}/run.log" \
    --warnings-file "${out_dir}/warnings.log" \
    --output-detail debug \
    --deterministic \
    --hide-progress \
    --min-gene-set-size 1 \
    --filter-gene-set-p 1 \
    --max-gene-set-read-p 1 \
    --no-filter-negative \
    --prune-gene-sets 1.1 \
    --weighted-prune-gene-sets 1.1
}

run_pigean original_markers
run_pigean top_absolute_expression
run_pigean top_specific_fc
run_pigean top_specific_logp
```

## Batch Workflow Config

The same command can be used through `run_all_cell_state_workflow.py` with a
PIGEAN template. The batch runner supplies `{gmt}`, `{multi_y}`, `{out_dir}`,
`{out}`, `{method}`, and `{extra_args}` and runs the template separately for
each signature method in each group.

```yaml
pigean:
  pigean_command: ../.venv/bin/python
  pigean_command_template: >-
    env PYTHONPATH=/Users/flannick/codex-workspace/analysis/pigean_optimize/pigean/src
    {pigean} -m pigean betas
    --X-in {gmt}
    --multi-y-in {multi_y}
    --multi-y-id-col Gene
    --multi-y-pheno-col Trait_Internal
    --multi-y-log-bf-col Direct
    --multi-y-combined-col Combined
    --multi-y-prior-col Indirect
    --multi-y-trait-blacklist-in results/pigean_all_cell_states_full_universe_default_no_hpo_exomes_gcat_orphanet/trait_blacklist_hp_exomes_gcat_orphanet.txt
    --gene-universe-in ../resources/pigean/data/reference/NCBI37.3.plink.gene.loc
    --gene-universe-id-col 6
    --gene-universe-no-header
    --gene-set-stats-out {out}
    --params-out {out_dir}/params.out.gz
    --log-file {out_dir}/run.log
    --warnings-file {out_dir}/warnings.log
    --output-detail debug
    --deterministic
    --hide-progress
    --min-gene-set-size 1
    --filter-gene-set-p 1
    --max-gene-set-read-p 1
    --no-filter-negative
    --prune-gene-sets 1.1
    --weighted-prune-gene-sets 1.1
    {extra_args}
  multi_y_input: ../resources/pigean/data/large/all.gene_stats.large.gt1.out.gz
  methods: [original_markers, top_absolute_expression, top_specific_fc, top_specific_logp]
  extra_args: ""
  dry_run: false
```
