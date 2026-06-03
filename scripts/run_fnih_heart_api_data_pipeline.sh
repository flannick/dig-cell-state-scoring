#!/usr/bin/env bash
set -euo pipefail

# Example API-minimal heart run using local FNIH heart expression/metadata
# and LIGER program files.
#
# Input files used:
#   data/external/fnih/heart/norm_counts.tsv.gz
#   data/external/fnih/heart/sample_metadata.tsv.gz
#   cell_state_de/dat/heart/heart_cell_state_markers.gmt
#   cell_state_de/dat/api/curated_cell_state_manifest.tsv
#   data/external/liger/heart/<cell type>/gene_loadings.tsv
#   data/external/liger/heart/<cell type>/cell_scores.tsv, when present
#
# Final API-ready files written to results/fnih/heart/out:
#   cell_state_expression.tsv.gz
#   cell_type_expression.tsv.gz
#   program_expression.tsv.gz
#   program_factor_metadata.tsv.gz
#   program_gene_loadings.tsv.gz
#   program_pigean_trait_results.tsv.gz
#   cell_state_pigean_trait_results.tsv.gz
#   program_state_heatmap.tsv.gz
#
# Intermediates are written only under results/fnih/heart/tmp:
#   input_metadata.tsv.gz and program_cell_type_map.tsv
#   rank_10x/ sparse matrix converted from norm_counts.tsv.gz
#   combined signatures/manifests, scoring, expression, matching, and logs.
#
# Default behavior uses CELL_SAMPLE_FRACTION=0.10 for a first-pass heart run.
# Heart has roughly 1.38M cells, so use CELL_SAMPLE_FRACTION=1.0 only for a full run.
# Override examples:
#   TISSUE_ROOT=/path/to/heart bash cell_state_de/scripts/run_fnih_heart_api_data_pipeline.sh
#   TOP_PROGRAM_GENES=200 bash cell_state_de/scripts/run_fnih_heart_api_data_pipeline.sh
#   CELL_SAMPLE_FRACTION=1.0 bash cell_state_de/scripts/run_fnih_heart_api_data_pipeline.sh   # full run
#   CELL_SAMPLE_FRACTION=0.10 CELL_SAMPLE_SEED=1 bash cell_state_de/scripts/run_fnih_heart_api_data_pipeline.sh

PYTHON_CMD="${PYTHON_CMD:-../.venv/bin/python}"
TISSUE_ROOT="${TISSUE_ROOT:-results/fnih/heart}"
CELL_STATE_DE_DIR="${CELL_STATE_DE_DIR:-cell_state_de}"
CELL_SAMPLE_FRACTION="${CELL_SAMPLE_FRACTION:-0.10}"
CELL_SAMPLE_SEED="${CELL_SAMPLE_SEED:-1}"
CELL_INCLUDE="${CELL_INCLUDE:-}"
PIGEAN_ENABLE="${PIGEAN_ENABLE:-1}"
PIGEAN_PYTHONPATH="${PIGEAN_PYTHONPATH:-/Users/flannick/codex-workspace/analysis/pigean_optimize/pigean/src}"
PIGEAN_MULTI_Y_IN="${PIGEAN_MULTI_Y_IN:-../resources/pigean/data/large/all.gene_stats.large.gt1.out.gz}"
PIGEAN_TRAIT_BLACKLIST="${PIGEAN_TRAIT_BLACKLIST:-auto}"
PIGEAN_GENE_UNIVERSE="${PIGEAN_GENE_UNIVERSE:-../resources/pigean/data/reference/NCBI37.3.plink.gene.loc}"
DEFAULT_STATE_WEIGHT_TYPE="${DEFAULT_STATE_WEIGHT_TYPE:-gradient_percentile_squared}"
EXPRESSION_VALUE_TYPE="${EXPRESSION_VALUE_TYPE:-auto}"
FORCE="${FORCE:-0}"
mkdir -p "${TISSUE_ROOT}/tmp"

TISSUE_ROOT="${TISSUE_ROOT}" "${PYTHON_CMD}" - <<'PY'
from pathlib import Path
import os
import pandas as pd

root = Path(os.environ.get('TISSUE_ROOT', 'results/fnih/heart'))
tmp = root / 'tmp'
tmp.mkdir(parents=True, exist_ok=True)
metadata = pd.read_csv('data/external/fnih/heart/sample_metadata.tsv.gz', sep='	', compression='infer', low_memory=False)

cell_type_map = {
    'ventricular cardiomyocyte': 'ventricular_cardiomyocyte',
    'atrial cardiomyocyte': 'atrial_cardiomyocyte',
    'fibroblast': 'fibroblast',
    'endothelial': 'endothelial',
    'pericyte': 'pericyte',
    'macrophage': 'macrophage',
    'T cell': 't_cell',
    'dendritic cell': 'dendritic_cell',
    'endocardial cell': 'endocardial_cell',
    'SMC': 'smc',
    'adipocyte': 'adipocyte',
    'neuron': 'neuron',
    'NK cell': 'nk_cell',
    'LEC': 'lec',
    'mast cell': 'mast_cell',
    'epicardial cell': 'epicardial_cell',
    'B cell': 'b_cell',
}
cell_type_source = metadata['cell_type__kp'].astype(str)
cell_type = cell_type_source.map(cell_type_map).fillna(
    cell_type_source.str.lower().str.replace(r'[^a-z0-9]+', '_', regex=True).str.strip('_')
)
out = pd.DataFrame({
    'cell_id': metadata['ID'].astype(str),
    'map_id': 'fnih_heart',
    'tissue': 'heart',
    'cell_type': cell_type,
    'donor_id': metadata['donor_id'].astype(str),
    'sample_id': metadata['biosample_id'].astype(str),
})
out.to_csv(tmp / 'input_metadata.tsv.gz', sep='	', index=False, compression='gzip')
program_map = pd.DataFrame([
    ('B cell', 'b_cell'),
    ('LEC', 'lec'),
    ('endocardial cell', 'endocardial_cell'),
    ('endothelial', 'endothelial'),
    ('pericyte', 'pericyte'),
], columns=['program_dir', 'cell_type'])
program_map.to_csv(tmp / 'program_cell_type_map.tsv', sep='	', index=False)
PY

TISSUE_ROOT="${TISSUE_ROOT}" \
TISSUE_ID="heart" \
EXPRESSION_TSV="data/external/fnih/heart/norm_counts.tsv.gz" \
EXPRESSION_VALUE_TYPE="${EXPRESSION_VALUE_TYPE}" \
METADATA="${TISSUE_ROOT}/tmp/input_metadata.tsv.gz" \
STATES_GMT="${CELL_STATE_DE_DIR}/dat/heart/heart_cell_state_markers.gmt" \
STATE_MANIFEST="${CELL_STATE_DE_DIR}/dat/api/curated_cell_state_manifest.tsv" \
CELL_ID_COL="cell_id" \
TISSUE_COL="tissue" \
CELL_TYPE_COL="cell_type" \
DONOR_COL="donor_id" \
SAMPLE_COL="sample_id" \
MAP_ID_COL="map_id" \
DATASET="scRNA" \
MODEL="mouse_msigdb" \
PROGRAM_ROOT="data/external/liger/heart" \
PROGRAM_CELL_TYPE_MAP="${TISSUE_ROOT}/tmp/program_cell_type_map.tsv" \
TOP_PROGRAM_GENES="${TOP_PROGRAM_GENES:-100}" \
CELL_SAMPLE_FRACTION="${CELL_SAMPLE_FRACTION}" \
CELL_SAMPLE_SEED="${CELL_SAMPLE_SEED}" \
CELL_INCLUDE="${CELL_INCLUDE}" \
PIGEAN_ENABLE="${PIGEAN_ENABLE}" \
PIGEAN_PYTHONPATH="${PIGEAN_PYTHONPATH}" \
PIGEAN_MULTI_Y_IN="${PIGEAN_MULTI_Y_IN}" \
PIGEAN_TRAIT_BLACKLIST="${PIGEAN_TRAIT_BLACKLIST}" \
PIGEAN_GENE_UNIVERSE="${PIGEAN_GENE_UNIVERSE}" \
DEFAULT_STATE_WEIGHT_TYPE="${DEFAULT_STATE_WEIGHT_TYPE}" \
FORCE="${FORCE}" \
PYTHON_CMD="${PYTHON_CMD}" \
CELL_STATE_DE_DIR="${CELL_STATE_DE_DIR}" \
"${CELL_STATE_DE_DIR}/scripts/run_tissue_api_data_pipeline.sh"
