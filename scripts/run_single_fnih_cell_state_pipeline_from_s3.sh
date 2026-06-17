#!/usr/bin/env bash
set -euo pipefail

# Download FNIH single-cell normalized-expression bundles from S3, convert to
# sparse 10x-like inputs, run cell_state_de scoring/expression, and remove large
# staged matrix files after successful completion.
#
# Examples:
#   scripts/run_fnih_cell_state_pipeline_from_s3.sh
#   DATASET=heart scripts/run_fnih_cell_state_pipeline_from_s3.sh
#   CLEANUP_LARGE_FILES=0 DATASET=heart scripts/run_fnih_cell_state_pipeline_from_s3.sh

PYTHON_CMD="${PYTHON_CMD:-python}"
S3_PREFIX="${S3_PREFIX:-s3://dig-analysis-data/single_cell}"
SCORING_ROOT="${SCORING_ROOT:-dig-cell-state-scoring}"
DATA_ROOT="${DATA_ROOT:-inputs}"
RESULTS_ROOT="${RESULTS_ROOT:-results/fnih_cell_state_scoring}"
S3_DATASET_NAME="${S3_DATASET_NAME:-}"
TISSUE="${TISSUE:-}"
RUN_EXPRESSION="${RUN_EXPRESSION:-1}"
CLEANUP_LARGE_FILES="${CLEANUP_LARGE_FILES:-1}"
VALUE_TYPE="${VALUE_TYPE:-auto}"
RANK_VALUE_TYPE="${RANK_VALUE_TYPE:-auto}"
EXPRESSION_VALUE_TYPE="${EXPRESSION_VALUE_TYPE:-auto}"
PROGRESS_EVERY_CELLS="${PROGRESS_EVERY_CELLS:-10000}"
ALLOW_SMALL_RANK_UNIVERSE="${ALLOW_SMALL_RANK_UNIVERSE:-0}"
AWS_CMD="${AWS_CMD:-aws}"
QC_GMT="${QC_GMT:-${SCORING_ROOT}/dat/qc/cmdkp_all_tissues_minimal_bad_cell_qc_signatures.gmt}"
STATE_MANIFEST="${STATE_MANIFEST:-${SCORING_ROOT}/dat/api/curated_cell_state_manifest_v2.tsv}"

if [[ ! -n "${S3_DATASET_NAME}" ]]; then
  echo "S3 Dataset Name not found or empty: ${S3_DATASET_NAME}" >&2
  exit 1
fi

if [[ ! -n "${TISSUE}" ]]; then
  echo "Tissue not found or empty: ${TISSUE}" >&2
  exit 1
fi

mkdir -p "${RESULTS_ROOT}"

log() {
  printf '[%s] %s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "$*"
}

safe_tissue_for_dat() {
  case "$1" in
    var) printf 'vat' ;;
    *) printf '%s' "$1" ;;
  esac
}

maybe_allow_small_rank_arg=()
if [[ "${ALLOW_SMALL_RANK_UNIVERSE}" == "1" ]]; then
  maybe_allow_small_rank_arg=(--allow-small-rank-universe)
fi

dat_tissue="$(safe_tissue_for_dat "${TISSUE}")"
tissue_dir="${DATA_ROOT}/${TISSUE}"
out_dir="${RESULTS_ROOT}/${TISSUE}"
input_dir="${out_dir}/inputs"
rank_dir="${out_dir}/rank_10x"
scoring_dir="${out_dir}/scoring"
expression_dir="${out_dir}/expression"
logs_dir="${out_dir}/logs"
state_gmt="${SCORING_ROOT}/dat/${dat_tissue}/${dat_tissue}_cell_state_markers.gmt"
scoring_done="${scoring_dir}/.complete"
expression_done="${expression_dir}/.complete"

mkdir -p "${tissue_dir}" "${input_dir}" "${rank_dir}" "${scoring_dir}" "${expression_dir}" "${logs_dir}"

if [[ ! -s "${state_gmt}" ]]; then
  log "Skipping ${TISSUE}: no curated state GMT at ${state_gmt}"
  printf 'tissue\treason\tstate_gmt\n%s\tmissing_state_gmt\t%s\n' "${TISSUE}" "${state_gmt}" > "${out_dir}/skipped.tsv"
  continue
fi

log "=== ${TISSUE} (${S3_DATASET_NAME}) ==="

if [[ -s "${scoring_done}" && ( "${RUN_EXPRESSION}" != "1" || -s "${expression_done}" ) ]]; then
  log "Skipping ${TISSUE}: requested outputs already complete"
  continue
fi

if [[ ! -s "${tissue_dir}/norm_counts.tsv.gz" || ! -s "${tissue_dir}/sample_metadata.tsv.gz" ]]; then
  log "Downloading ${TISSUE} from ${S3_PREFIX}/${S3_DATASET_NAME}/"
  "${AWS_CMD}" s3 cp "${S3_PREFIX}/${S3_DATASET_NAME}/" "${tissue_dir}/" --recursive \
    > "${logs_dir}/download.stdout.log" \
    2> "${logs_dir}/download.stderr.log"
else
  log "Download already present for ${TISSUE}; skipping S3 copy"
fi

if [[ ! -s "${tissue_dir}/norm_counts.tsv.gz" ]]; then
  echo "Missing downloaded norm_counts.tsv.gz for ${TISSUE}" >&2
  exit 1
fi
if [[ ! -s "${tissue_dir}/sample_metadata.tsv.gz" ]]; then
  echo "Missing downloaded sample_metadata.tsv.gz for ${TISSUE}" >&2
  exit 1
fi

metadata="${input_dir}/metadata.tsv.gz"
expression_query="${input_dir}/marker_query_expression.tsv.gz"
marker_genes="${input_dir}/state_and_qc_marker_genes.txt"

if [[ ! -s "${metadata}" ]]; then
  log "Preparing normalized metadata for ${TISSUE}"
  TISSUE="${TISSUE}" DAT_TISSUE="${dat_tissue}" METADATA_IN="${tissue_dir}/sample_metadata.tsv.gz" METADATA_OUT="${metadata}" \
    "${PYTHON_CMD}" - <<'PY'
import os
import pandas as pd

tissue = os.environ["TISSUE"]
dat_tissue = os.environ["DAT_TISSUE"]
metadata_in = os.environ["METADATA_IN"]
metadata_out = os.environ["METADATA_OUT"]
meta = pd.read_csv(metadata_in, sep="\t", compression="infer", low_memory=False)
if "cell_id" not in meta.columns:
    if "ID" in meta.columns:
        meta["cell_id"] = meta["ID"].astype(str)
    elif "barcode" in meta.columns:
        meta["cell_id"] = meta["barcode"].astype(str)
    else:
        raise SystemExit("Metadata has no cell ID column")
if "map_id" not in meta.columns:
    meta["map_id"] = meta.get("DI:Dataset", f"fnih_{tissue}")
if "tissue" not in meta.columns:
    meta["tissue"] = dat_tissue
else:
    meta["tissue"] = dat_tissue
if "annotated_cell_type" not in meta.columns:
    if "cell_type__kp" in meta.columns:
        meta["annotated_cell_type"] = meta["cell_type__kp"].astype(str)
    elif "Cell_Type" in meta.columns:
        meta["annotated_cell_type"] = meta["Cell_Type"].astype(str)
    elif "Annotations:celltype" in meta.columns:
        meta["annotated_cell_type"] = meta["Annotations:celltype"].astype(str)
    else:
        raise SystemExit("Metadata has no usable cell-type column")
if "donor_id" not in meta.columns:
    if "SI:Subject_ID" in meta.columns:
        meta["donor_id"] = meta["SI:Subject_ID"].astype(str)
    elif "biosample_id" in meta.columns:
        meta["donor_id"] = meta["biosample_id"].astype(str)
    else:
        meta["donor_id"] = meta["cell_id"].astype(str)
if "sample_id" not in meta.columns:
    if "biosample_id" in meta.columns:
        meta["sample_id"] = meta["biosample_id"].astype(str)
    else:
        meta["sample_id"] = meta["donor_id"].astype(str)
meta.to_csv(metadata_out, sep="\t", index=False, compression="gzip")
print(f"Wrote {metadata_out} with {len(meta)} cells")
PY
fi

if [[ ! -s "${marker_genes}" ]]; then
  log "Collecting state/QC marker genes for query expression matrix"
  STATE_GMT="${state_gmt}" QC_GMT="${QC_GMT}" MARKER_GENES_OUT="${marker_genes}" "${PYTHON_CMD}" - <<'PY'
import os
from pathlib import Path

genes = []
for env_name in ["STATE_GMT", "QC_GMT"]:
    path = Path(os.environ[env_name])
    if not path.exists():
        continue
    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            genes.extend(parts[2:])
seen = []
for gene in genes:
    gene = gene.strip()
    if gene and gene not in seen:
        seen.append(gene)
Path(os.environ["MARKER_GENES_OUT"]).write_text("\n".join(seen) + "\n", encoding="utf-8")
print(f"Wrote {len(seen)} marker genes")
PY
fi

if [[ ! -s "${rank_dir}/matrix.mtx.gz" ]]; then
  log "Converting dense normalized matrix to sparse 10x-like input for ${TISSUE}"
  "${PYTHON_CMD}" ${SCORING_ROOT}/scripts/convert_expression_tsv_to_sparse_10x.py \
    --matrix-tsv "${tissue_dir}/norm_counts.tsv.gz" \
    --out-dir "${rank_dir}" \
    --value-type "${VALUE_TYPE}" \
    > "${logs_dir}/convert_to_sparse.stdout.log" \
    2> "${logs_dir}/convert_to_sparse.stderr.log"
else
  log "Sparse converted matrix already present for ${TISSUE}; skipping conversion"
fi

if [[ ! -s "${expression_query}" ]]; then
  log "Extracting marker-only expression query table for ${TISSUE}"
  NORM_COUNTS="${tissue_dir}/norm_counts.tsv.gz" MARKER_GENES="${marker_genes}" QUERY_OUT="${expression_query}" \
    "${PYTHON_CMD}" - <<'PY'
import gzip
import os
from pathlib import Path

norm_counts = Path(os.environ["NORM_COUNTS"])
markers = {line.strip() for line in Path(os.environ["MARKER_GENES"]).read_text().splitlines() if line.strip()}
out_path = Path(os.environ["QUERY_OUT"])
open_in = gzip.open if str(norm_counts).endswith(".gz") else open
with open_in(norm_counts, "rt", encoding="utf-8") as inp, gzip.open(out_path, "wt", encoding="utf-8") as out:
    header = inp.readline().rstrip("\n").split("\t")
    cells = header[1:]
    out.write("cell_id\tgene\texpression\n")
    n_rows = 0
    n_values = 0
    for line in inp:
        parts = line.rstrip("\n").split("\t")
        gene = parts[0]
        if gene not in markers:
            continue
        n_rows += 1
        for cell, value in zip(cells, parts[1:]):
            if value and value != "0" and value.lower() != "nan":
                out.write(f"{cell}\t{gene}\t{value}\n")
                n_values += 1
print(f"Wrote {n_rows} marker genes and {n_values} nonzero values to {out_path}")
PY
fi

if [[ ! -s "${scoring_done}" ]]; then
  log "Running state scoring for ${TISSUE}"
  "${PYTHON_CMD}" ${SCORING_ROOT}/scripts/run_cmdkp_state_scoring.py \
    --rank-10x-dir "${rank_dir}" \
    --rank-value-type "${RANK_VALUE_TYPE}" \
    --expression-kind log1p_normalized \
    --expression-matrix "${expression_query}" \
    --cell-metadata "${metadata}" \
    --states-gmt "${state_gmt}" \
    --state-manifest "${STATE_MANIFEST}" \
    --qc-states-gmt "${QC_GMT}" \
    --map-id-col map_id \
    --tissue-col tissue \
    --cell-type-col annotated_cell_type \
    --donor-col donor_id \
    --sample-col sample_id \
    --progress-every-cells "${PROGRESS_EVERY_CELLS}" \
    --legacy-selected-gene-summaries skip \
    --allow-acceptance-failures \
    ${maybe_allow_small_rank_arg:+"${maybe_allow_small_rank_arg[@]}"} \
    --out-dir "${scoring_dir}" \
    > "${logs_dir}/state_scoring.stdout.log" \
    2> "${logs_dir}/state_scoring.stderr.log"
  date '+%Y-%m-%d %H:%M:%S' > "${scoring_done}"
else
  log "State scoring already complete for ${TISSUE}"
fi

if [[ "${RUN_EXPRESSION}" == "1" && ! -s "${expression_done}" ]]; then
  log "Running normalized expression summaries for ${TISSUE}"
  "${PYTHON_CMD}" ${SCORING_ROOT}/scripts/summarize_state_expression.py \
    --raw-10x-dir "${rank_dir}" \
    --expression-value-type "${EXPRESSION_VALUE_TYPE}" \
    --metadata "${metadata}" \
    --cell-state-activity "${scoring_dir}/cell_state_activity.tsv.gz" \
    --states-gmt "${state_gmt}" \
    --parent-group-cols tissue,annotated_cell_type \
    --cell-type-col annotated_cell_type \
    --donor-col donor_id \
    --donor-expression-genes none \
    --no-write-donor-state-expression \
    --out-dir "${expression_dir}" \
    > "${logs_dir}/state_expression.stdout.log" \
    2> "${logs_dir}/state_expression.stderr.log"
  date '+%Y-%m-%d %H:%M:%S' > "${expression_done}"
elif [[ "${RUN_EXPRESSION}" == "1" ]]; then
  log "Expression summaries already complete for ${TISSUE}"
fi

if [[ "${CLEANUP_LARGE_FILES}" == "1" && -s "${scoring_done}" && ( "${RUN_EXPRESSION}" != "1" || -s "${expression_done}" ) ]]; then
  log "Cleaning large staged matrix files for ${TISSUE}"
  rm -f "${tissue_dir}/norm_counts.tsv.gz"
  rm -f "${rank_dir}/matrix.mtx.gz"
fi

log "Finished ${TISSUE}"
