#!/usr/bin/env bash
set -euo pipefail

# API-minimal cell-state data pipeline for one tissue.
#
# Required configuration, passed as environment variables:
#   TISSUE_ROOT              Directory for this tissue run. Final outputs go to ${TISSUE_ROOT}/out; intermediates go to ${TISSUE_ROOT}/tmp.
#   TISSUE_ID                Stable tissue id, e.g. pancreas.
#   EXPRESSION_TSV           Wide expression TSV/TSV.GZ. First column is gene, remaining columns are cells.
#   EXPRESSION_VALUE_TYPE    Matrix value type: raw_counts, linear_cp10k, log1p_cp10k, linear_normalized, log1p_normalized, scaled, or auto.
#   METADATA                 Cell metadata TSV/TSV.GZ.
#   STATES_GMT               Curated biological state GMT for this tissue.
#   STATE_MANIFEST           Manifest TSV with at least state_name, tissue, cell_type, and state_class.
#
# Optional configuration:
#   CELL_ID_COL              Metadata cell ID column. Default: cell_id.
#   TISSUE_COL               Metadata tissue column. Default: tissue.
#   CELL_TYPE_COL            Metadata cell type column. Default: cell_type.
#   DONOR_COL                Metadata donor column. Default: donor_id.
#   SAMPLE_COL               Metadata sample/library column. Default: sample_id.
#   MAP_ID_COL               Metadata map/dataset column. Default: map_id.
#   DATASET                 Portal dataset label for program endpoints. Default: scRNA.
#   MODEL                   Portal model label for program endpoints. Default: cell_state_de.
#   PROGRAM_ROOT             Directory containing cell-type program directories with gene_loadings.tsv and optional cell_scores.tsv. Nested tissue/cell-type layouts are supported.
#   PROGRAM_CELL_TYPE_MAP    Optional TSV with columns program_dir and cell_type. If absent, program directory names are snake-cased.
#   TOP_PROGRAM_GENES        Number of top genes per factor used as program signatures. Default: 100.
#   QC_GMT                  Optional QC GMT. If absent, no QC states are scored.
#   CELL_STATE_PIGEAN        Optional precomputed curated-state PIGEAN result TSV/TSV.GZ. If unset, PIGEAN is run per cell type.
#   PROGRAM_PIGEAN           Optional precomputed program PIGEAN result TSV/TSV.GZ. If unset, PIGEAN is run per cell type.
#   PIGEAN_ENABLE            Run PIGEAN when precomputed outputs are not supplied. Default: 1.
#   PIGEAN_PYTHONPATH        PIGEAN source path. Default: /Users/flannick/codex-workspace/analysis/pigean_optimize/pigean/src.
#   PIGEAN_MULTI_Y_IN        Multi-trait gene stats file. Default: ../resources/pigean/data/large/all.gene_stats.large.gt1.out.gz.
#   PIGEAN_TRAIT_BLACKLIST   Trait blacklist path, or auto to generate HP_/exomes_ blacklist from PIGEAN_MULTI_Y_IN. Default: auto.
#   PIGEAN_GENE_UNIVERSE     Gene universe file. Default: ../resources/pigean/data/reference/NCBI37.3.plink.gene.loc.
#   PYTHON_CMD               Python executable. Default: ../.venv/bin/python.
#   CELL_STATE_DE_DIR        Repo directory. Default: directory above this script.
#   CELL_SAMPLE_FRACTION      Random cell fraction to keep while converting EXPRESSION_TSV. Default: 1.0.
#   CELL_SAMPLE_SEED          Random seed for CELL_SAMPLE_FRACTION. Default: 1.
#   CELL_INCLUDE              Optional one-column cell ID list to keep during conversion.
#   ALLOW_SMALL_RANK_UNIVERSE Set to 1 to pass --allow-small-rank-universe for smoke tests. Default: 0.
#
# Final files written under ${TISSUE_ROOT}/out:
#   cell_state_expression.tsv.gz
#   cell_type_expression.tsv.gz
#   program_expression.tsv.gz
#   program_factor_metadata.tsv.gz
#   program_gene_loadings.tsv.gz
#   program_pigean_trait_results.tsv.gz
#   cell_state_pigean_trait_results.tsv.gz
#   program_state_heatmap.tsv.gz
#
# Intermediate files written under ${TISSUE_ROOT}/tmp only:
#   rank_10x/                                Sparse matrix converted from EXPRESSION_TSV.
#   metadata.tsv.gz                          Normalized metadata with pipeline columns.
#   minimal_expression.tsv.gz                One dummy gene expression table required by the scorer while sparse ranking provides real scores.
#   combined_signatures.gmt                  Curated states plus program signatures.
#   combined_signature_manifest.tsv          Scope/class manifest for combined signatures.
#   signature_kind.tsv                       Signature kind lookup.
#   program_source_manifest.tsv              Per-cell-type program input manifest.
#   program_inputs/                          Program loadings/activity staged for matching.
#   split_state_gmts/                        Curated state GMTs split by cell type.
#   split_program_gmts/                      Program signature GMTs split by cell type.
#   pigean/                                  Per-cell-type PIGEAN runs and combined tables.
#   scoring/                                 API-minimal activity output and run summary.
#   expression/                              API-minimal state/cell-type expression output.
#   program_state_matches/                   API-minimal program-state heatmap and labels.
#   logs/                                    Command stdout/stderr logs.

PYTHON_CMD="${PYTHON_CMD:-../.venv/bin/python}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CELL_STATE_DE_DIR="${CELL_STATE_DE_DIR:-$(cd "${SCRIPT_DIR}/.." && pwd)}"
TISSUE_ROOT="${TISSUE_ROOT:?set TISSUE_ROOT}"
TISSUE_ID="${TISSUE_ID:?set TISSUE_ID}"
EXPRESSION_TSV="${EXPRESSION_TSV:?set EXPRESSION_TSV}"
EXPRESSION_VALUE_TYPE="${EXPRESSION_VALUE_TYPE:?set EXPRESSION_VALUE_TYPE}"
METADATA="${METADATA:?set METADATA}"
STATES_GMT="${STATES_GMT:?set STATES_GMT}"
STATE_MANIFEST="${STATE_MANIFEST:?set STATE_MANIFEST}"
CELL_ID_COL="${CELL_ID_COL:-cell_id}"
TISSUE_COL="${TISSUE_COL:-tissue}"
CELL_TYPE_COL="${CELL_TYPE_COL:-cell_type}"
DONOR_COL="${DONOR_COL:-donor_id}"
SAMPLE_COL="${SAMPLE_COL:-sample_id}"
MAP_ID_COL="${MAP_ID_COL:-map_id}"
DATASET="${DATASET:-scRNA}"
MODEL="${MODEL:-cell_state_de}"
PROGRAM_ROOT="${PROGRAM_ROOT:-}"
PROGRAM_CELL_TYPE_MAP="${PROGRAM_CELL_TYPE_MAP:-}"
TOP_PROGRAM_GENES="${TOP_PROGRAM_GENES:-100}"
QC_GMT="${QC_GMT:-}"
CELL_STATE_PIGEAN="${CELL_STATE_PIGEAN:-}"
PROGRAM_PIGEAN="${PROGRAM_PIGEAN:-}"
PROGRESS_EVERY_CELLS="${PROGRESS_EVERY_CELLS:-10000}"
CELL_SAMPLE_FRACTION="${CELL_SAMPLE_FRACTION:-1.0}"
CELL_SAMPLE_SEED="${CELL_SAMPLE_SEED:-1}"
CELL_INCLUDE="${CELL_INCLUDE:-}"
ALLOW_SMALL_RANK_UNIVERSE="${ALLOW_SMALL_RANK_UNIVERSE:-0}"
PIGEAN_ENABLE="${PIGEAN_ENABLE:-1}"
PIGEAN_PYTHONPATH="${PIGEAN_PYTHONPATH:-/Users/flannick/codex-workspace/analysis/pigean_optimize/pigean/src}"
PIGEAN_MULTI_Y_IN="${PIGEAN_MULTI_Y_IN:-../resources/pigean/data/large/all.gene_stats.large.gt1.out.gz}"
PIGEAN_TRAIT_BLACKLIST="${PIGEAN_TRAIT_BLACKLIST:-auto}"
PIGEAN_GENE_UNIVERSE="${PIGEAN_GENE_UNIVERSE:-../resources/pigean/data/reference/NCBI37.3.plink.gene.loc}"
DEFAULT_STATE_WEIGHT_TYPE="${DEFAULT_STATE_WEIGHT_TYPE:-gradient_percentile_squared}"
FORCE="${FORCE:-0}"

OUT_DIR="${TISSUE_ROOT}/out"
TMP_DIR="${TISSUE_ROOT}/tmp"
LOG_DIR="${TMP_DIR}/logs"
mkdir -p "${OUT_DIR}" "${TMP_DIR}" "${LOG_DIR}" "${TMP_DIR}/program_inputs" "${TMP_DIR}/split_state_gmts" "${TMP_DIR}/split_program_gmts" "${TMP_DIR}/program_state_matches" "${TMP_DIR}/pigean"

cat > "${TMP_DIR}/run_parameters.txt" <<EOF
TISSUE_ROOT=${TISSUE_ROOT}
TISSUE_ID=${TISSUE_ID}
EXPRESSION_TSV=${EXPRESSION_TSV}
EXPRESSION_VALUE_TYPE=${EXPRESSION_VALUE_TYPE}
METADATA=${METADATA}
STATES_GMT=${STATES_GMT}
STATE_MANIFEST=${STATE_MANIFEST}
CELL_ID_COL=${CELL_ID_COL}
TISSUE_COL=${TISSUE_COL}
CELL_TYPE_COL=${CELL_TYPE_COL}
DONOR_COL=${DONOR_COL}
SAMPLE_COL=${SAMPLE_COL}
MAP_ID_COL=${MAP_ID_COL}
DATASET=${DATASET}
MODEL=${MODEL}
PROGRAM_ROOT=${PROGRAM_ROOT}
PROGRAM_CELL_TYPE_MAP=${PROGRAM_CELL_TYPE_MAP}
TOP_PROGRAM_GENES=${TOP_PROGRAM_GENES}
QC_GMT=${QC_GMT}
CELL_STATE_PIGEAN=${CELL_STATE_PIGEAN}
PROGRAM_PIGEAN=${PROGRAM_PIGEAN}
CELL_SAMPLE_FRACTION=${CELL_SAMPLE_FRACTION}
CELL_SAMPLE_SEED=${CELL_SAMPLE_SEED}
CELL_INCLUDE=${CELL_INCLUDE}
ALLOW_SMALL_RANK_UNIVERSE=${ALLOW_SMALL_RANK_UNIVERSE}
PIGEAN_ENABLE=${PIGEAN_ENABLE}
PIGEAN_MULTI_Y_IN=${PIGEAN_MULTI_Y_IN}
PIGEAN_TRAIT_BLACKLIST=${PIGEAN_TRAIT_BLACKLIST}
PIGEAN_GENE_UNIVERSE=${PIGEAN_GENE_UNIVERSE}
DEFAULT_STATE_WEIGHT_TYPE=${DEFAULT_STATE_WEIGHT_TYPE}
FORCE=${FORCE}
EOF

if [[ ! -s "${TMP_DIR}/rank_10x/matrix.mtx.gz" ]]; then
  echo "[$(date)] Converting expression TSV to sparse 10x under ${TMP_DIR}/rank_10x" | tee "${LOG_DIR}/pipeline.progress.log"
  INCLUDE_ARGS=()
  if [[ -n "${CELL_INCLUDE}" ]]; then
    INCLUDE_ARGS=(--cell-include "${CELL_INCLUDE}")
  fi
  "${PYTHON_CMD}" "${CELL_STATE_DE_DIR}/scripts/convert_expression_tsv_to_sparse_10x.py" \
    --matrix-tsv "${EXPRESSION_TSV}" \
    --out-dir "${TMP_DIR}/rank_10x" \
    --orientation gene_by_cell \
    --value-type "${EXPRESSION_VALUE_TYPE}" \
    --cell-sample-fraction "${CELL_SAMPLE_FRACTION}" \
    --cell-sample-seed "${CELL_SAMPLE_SEED}" \
    ${INCLUDE_ARGS[@]+"${INCLUDE_ARGS[@]}"} \
    > >(tee "${LOG_DIR}/convert_matrix.stdout.log") \
    2> >(tee "${LOG_DIR}/convert_matrix.stderr.log" >&2)
else
  echo "[$(date)] Reusing existing sparse 10x matrix at ${TMP_DIR}/rank_10x" | tee "${LOG_DIR}/pipeline.progress.log"
fi

echo "[$(date)] Preparing normalized metadata" | tee -a "${LOG_DIR}/pipeline.progress.log"
TISSUE_ID="${TISSUE_ID}" METADATA="${METADATA}" TMP_DIR="${TMP_DIR}" CELL_ID_COL="${CELL_ID_COL}" TISSUE_COL="${TISSUE_COL}" CELL_TYPE_COL="${CELL_TYPE_COL}" DONOR_COL="${DONOR_COL}" SAMPLE_COL="${SAMPLE_COL}" MAP_ID_COL="${MAP_ID_COL}" \
"${PYTHON_CMD}" - <<'PY'
import os
from pathlib import Path
import pandas as pd

def pick(frame, col, default):
    if col in frame.columns:
        return frame[col]
    return default

metadata = pd.read_csv(os.environ["METADATA"], sep="\t", compression="infer", low_memory=False)
barcodes_path = Path(os.environ["TMP_DIR"]) / "rank_10x" / "barcodes.tsv.gz"
if barcodes_path.exists():
    keep_cells = pd.read_csv(barcodes_path, sep="\t", header=None, compression="infer").iloc[:, 0].astype(str)
else:
    keep_cells = None
out = pd.DataFrame()
out["cell_id"] = metadata[os.environ["CELL_ID_COL"]].astype(str)
if keep_cells is not None:
    keep_set = set(keep_cells)
    metadata = metadata.loc[out["cell_id"].isin(keep_set)].copy()
    out = pd.DataFrame()
    out["cell_id"] = metadata[os.environ["CELL_ID_COL"]].astype(str)
out["map_id"] = pick(metadata, os.environ["MAP_ID_COL"], os.environ["TISSUE_ID"]).astype(str) if os.environ["MAP_ID_COL"] in metadata.columns else os.environ["TISSUE_ID"]
out["tissue"] = pick(metadata, os.environ["TISSUE_COL"], os.environ["TISSUE_ID"]).astype(str) if os.environ["TISSUE_COL"] in metadata.columns else os.environ["TISSUE_ID"]
out["cell_type"] = pick(metadata, os.environ["CELL_TYPE_COL"], "unknown").astype(str) if os.environ["CELL_TYPE_COL"] in metadata.columns else "unknown"
out["annotated_cell_type"] = out["cell_type"]
out["donor_id"] = pick(metadata, os.environ["DONOR_COL"], "unknown").astype(str) if os.environ["DONOR_COL"] in metadata.columns else "unknown"
out["sample_id"] = pick(metadata, os.environ["SAMPLE_COL"], out["donor_id"]).astype(str) if os.environ["SAMPLE_COL"] in metadata.columns else out["donor_id"]
Path(os.environ["TMP_DIR"]).mkdir(parents=True, exist_ok=True)
out.to_csv(Path(os.environ["TMP_DIR"]) / "metadata.tsv.gz", sep="\t", index=False, compression="gzip")
pd.DataFrame({"cell_id": out["cell_id"], "gene": "DUMMY_GENE", "expression": 0.0}).to_csv(Path(os.environ["TMP_DIR"]) / "minimal_expression.tsv.gz", sep="\t", index=False, compression="gzip")
PY

echo "[$(date)] Building combined curated-state and program GMT" | tee -a "${LOG_DIR}/pipeline.progress.log"
TISSUE_ID="${TISSUE_ID}" TMP_DIR="${TMP_DIR}" STATES_GMT="${STATES_GMT}" STATE_MANIFEST="${STATE_MANIFEST}" PROGRAM_ROOT="${PROGRAM_ROOT}" PROGRAM_CELL_TYPE_MAP="${PROGRAM_CELL_TYPE_MAP}" TOP_PROGRAM_GENES="${TOP_PROGRAM_GENES}" \
"${PYTHON_CMD}" - <<'PY'
import os, re
from pathlib import Path
import pandas as pd

def norm(value):
    value = str(value).strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return re.sub(r"_+", "_", value).strip("_")

def read_gmt(path):
    rows=[]
    with open(path) as handle:
        for line in handle:
            parts=line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                rows.append((parts[0], parts[1], [g for g in parts[2:] if g]))
    return rows

def write_gmt(rows, path):
    with open(path, "w") as handle:
        for name, desc, genes in rows:
            handle.write("\t".join([name, desc] + list(dict.fromkeys(map(str, genes)))) + "\n")

def normalize_cell_id(cell_id):
    text = str(cell_id)
    if "_" in text:
        tail = text.split("_", 1)[1]
        if "-" in tail:
            return tail
    return text

tmp = Path(os.environ["TMP_DIR"])
tissue = os.environ["TISSUE_ID"]
state_manifest = pd.read_csv(os.environ["STATE_MANIFEST"], sep="\t", low_memory=False)
if "state_id" in state_manifest.columns and "state_name" not in state_manifest.columns:
    state_manifest = state_manifest.rename(columns={"state_id": "state_name"})
if "cell_type_id" in state_manifest.columns and "cell_type" not in state_manifest.columns:
    state_manifest = state_manifest.rename(columns={"cell_type_id": "cell_type"})
if "tissue_id" in state_manifest.columns and "tissue" not in state_manifest.columns:
    state_manifest = state_manifest.rename(columns={"tissue_id": "tissue"})
state_manifest = state_manifest[state_manifest["tissue"].astype(str).eq(tissue)].copy()
state_ids = set(state_manifest["state_name"].astype(str))
curated_rows = [row for row in read_gmt(Path(os.environ["STATES_GMT"])) if row[0] in state_ids]
manifest_rows = []
for row in state_manifest.itertuples(index=False):
    manifest_rows.append({
        "state_name": str(getattr(row, "state_name")),
        "tissue": str(getattr(row, "tissue")),
        "cell_type": str(getattr(row, "cell_type")),
        "state_class": str(getattr(row, "state_class", "unknown")),
        "is_composite_required": str(getattr(row, "is_composite_required", "false")).lower(),
        "signature_kind": "curated_state",
    })
program_rows=[]; program_manifest=[]; program_source=[]
cell_type_map={}
map_path=os.environ.get("PROGRAM_CELL_TYPE_MAP", "")
if map_path and Path(map_path).exists():
    m=pd.read_csv(map_path, sep="\t")
    for r in m.itertuples(index=False):
        cell_type_map[str(getattr(r,"program_dir"))]=str(getattr(r,"cell_type"))
program_root=os.environ.get("PROGRAM_ROOT", "")
program_dirs=[]
if program_root:
    root_path=Path(program_root)
    if not root_path.exists():
        raise SystemExit(f"PROGRAM_ROOT was set but does not exist: {program_root}")
    program_dirs=sorted({p.parent for p in root_path.rglob("gene_loadings.tsv")})
    if not program_dirs:
        raise SystemExit(f"PROGRAM_ROOT was set but no gene_loadings.tsv files were found under: {program_root}")
    print(f"Discovered {len(program_dirs)} program directories under {program_root}", flush=True)
    for pdir in program_dirs:
        loadings_path=pdir/"gene_loadings.tsv"
        cell_type=cell_type_map.get(pdir.name, norm(pdir.name))
        loadings=pd.read_csv(loadings_path, sep="\t", index_col=0)
        renamed={}
        for factor in loadings.columns:
            factor_id=norm(str(factor).replace("Factor_", "factor_"))
            state_name=f"{tissue}_{cell_type}_program_{factor_id}"
            renamed[factor]=state_name
            top=loadings[factor].sort_values(ascending=False).head(int(os.environ["TOP_PROGRAM_GENES"]))
            genes=[str(g) for g,v in top.items() if pd.notna(v) and float(v) > 0]
            if genes:
                program_rows.append((state_name, f"type=program;cell_type={cell_type};source={pdir.name}", genes))
                program_manifest.append({"state_name": state_name, "tissue": tissue, "cell_type": cell_type, "state_class": "broad_function_gradient", "is_composite_required": "false", "signature_kind": "program"})
        loadings_out=tmp/"program_inputs"/f"{cell_type}.program_loadings.tsv.gz"
        loadings.rename(columns=renamed).reset_index(names="gene").to_csv(loadings_out, sep="\t", index=False, compression="gzip")
        scores_path=pdir/"cell_scores.tsv"
        activity_out=""
        if scores_path.exists():
            scores=pd.read_csv(scores_path, sep="\t", index_col=0).rename(columns=renamed)
            scores.index=[normalize_cell_id(x) for x in scores.index]
            activity_out=str(tmp/"program_inputs"/f"{cell_type}.program_cell_activity.tsv.gz")
            scores.reset_index(names="cell_id").to_csv(activity_out, sep="\t", index=False, compression="gzip")
        program_source.append({"cell_type": cell_type, "program_dir": pdir.name, "program_loadings": str(loadings_out), "program_cell_activity": activity_out})
combined_rows=curated_rows+program_rows
write_gmt(combined_rows, tmp/"combined_signatures.gmt")
pd.DataFrame(manifest_rows+program_manifest).to_csv(tmp/"combined_signature_manifest.tsv", sep="\t", index=False)
pd.DataFrame([{"state_name": x["state_name"], "cell_type": x["cell_type"], "signature_kind": x["signature_kind"]} for x in manifest_rows+program_manifest]).to_csv(tmp/"signature_kind.tsv", sep="\t", index=False)
pd.DataFrame(program_source, columns=["cell_type", "program_dir", "program_loadings", "program_cell_activity"]).to_csv(tmp/"program_source_manifest.tsv", sep="\t", index=False)
for cell_type, rows in state_manifest.groupby("cell_type"):
    ids=set(rows["state_name"].astype(str))
    write_gmt([row for row in curated_rows if row[0] in ids], tmp/"split_state_gmts"/f"{cell_type}.gmt")
program_by_cell_type = {}
for row in program_rows:
    match = re.match(rf"^{re.escape(tissue)}_(.+)_program_factor_", row[0])
    if match:
        program_by_cell_type.setdefault(match.group(1), []).append(row)
for cell_type, rows in program_by_cell_type.items():
    write_gmt(rows, tmp/"split_program_gmts"/f"{cell_type}.gmt")
PY

NO_QC_GMT="${TMP_DIR}/no_qc_states.gmt"
printf "__no_qc__\tplaceholder\t__NO_QC_GENE__\n" > "${NO_QC_GMT}"
QC_ARGS=(--qc-gmt "${NO_QC_GMT}")
if [[ -n "${QC_GMT}" && -s "${QC_GMT}" ]]; then
  QC_ARGS=(--qc-gmt "${QC_GMT}")
fi
SMALL_RANK_ARGS=()
if [[ "${ALLOW_SMALL_RANK_UNIVERSE}" == "1" || "${ALLOW_SMALL_RANK_UNIVERSE}" == "true" ]]; then
  SMALL_RANK_ARGS=(--allow-small-rank-universe)
fi
if [[ "${FORCE}" == "1" || ! -s "${TMP_DIR}/scoring/cell_state_activity.tsv.gz" ]]; then
  echo "[$(date)] Scoring cell states/programs" | tee -a "${LOG_DIR}/pipeline.progress.log"
  "${PYTHON_CMD}" "${CELL_STATE_DE_DIR}/scripts/run_cmdkp_state_scoring.py" \
    --rank-10x-dir "${TMP_DIR}/rank_10x" \
    --rank-value-type "${EXPRESSION_VALUE_TYPE}" \
    --expression-matrix "${TMP_DIR}/minimal_expression.tsv.gz" \
    --expression-kind linear_normalized \
    --cell-metadata "${TMP_DIR}/metadata.tsv.gz" \
    --states-gmt "${TMP_DIR}/combined_signatures.gmt" \
    --state-manifest "${TMP_DIR}/combined_signature_manifest.tsv" \
    --require-state-manifest \
    "${QC_ARGS[@]}" \
    --map-id-col map_id --tissue-col tissue --cell-type-col annotated_cell_type --donor-col donor_id --sample-col sample_id \
    --progress-every-cells "${PROGRESS_EVERY_CELLS}" \
    --legacy-selected-gene-summaries skip \
    --api-minimal-output \
    ${SMALL_RANK_ARGS[@]+"${SMALL_RANK_ARGS[@]}"} \
    --allow-acceptance-failures \
    --out-dir "${TMP_DIR}/scoring" \
    > >(tee "${LOG_DIR}/scoring.stdout.log") \
    2> >(tee "${LOG_DIR}/scoring.stderr.log" >&2)
else
  echo "[$(date)] Reusing existing scoring output at ${TMP_DIR}/scoring/cell_state_activity.tsv.gz" | tee -a "${LOG_DIR}/pipeline.progress.log"
fi

if [[ "${FORCE}" == "1" || ! -s "${TMP_DIR}/expression/all_gene_state_expression_specificity_cp10k.tsv.gz" || ! -s "${TMP_DIR}/expression/all_gene_cell_type_expression_cp10k.tsv.gz" ]]; then
  echo "[$(date)] Summarizing expression" | tee -a "${LOG_DIR}/pipeline.progress.log"
  "${PYTHON_CMD}" "${CELL_STATE_DE_DIR}/scripts/summarize_state_expression.py" \
    --raw-10x-dir "${TMP_DIR}/rank_10x" \
    --expression-value-type "${EXPRESSION_VALUE_TYPE}" \
    --metadata "${TMP_DIR}/metadata.tsv.gz" \
    --cell-state-activity "${TMP_DIR}/scoring/cell_state_activity.tsv.gz" \
    --states-gmt "${TMP_DIR}/combined_signatures.gmt" \
    --parent-group-cols tissue,annotated_cell_type \
    --cell-type-col annotated_cell_type \
    --donor-col donor_id \
    --donor-expression-genes none \
    --no-write-donor-state-expression \
    --api-minimal-output \
    --out-dir "${TMP_DIR}/expression" \
    > >(tee "${LOG_DIR}/expression.stdout.log") \
    2> >(tee "${LOG_DIR}/expression.stderr.log" >&2)
else
  echo "[$(date)] Reusing existing expression summaries under ${TMP_DIR}/expression" | tee -a "${LOG_DIR}/pipeline.progress.log"
fi

echo "[$(date)] Splitting curated-state and program expression/activity outputs" | tee -a "${LOG_DIR}/pipeline.progress.log"
TMP_DIR="${TMP_DIR}" "${PYTHON_CMD}" - <<'PY'
import os
from pathlib import Path
import pandas as pd
tmp=Path(os.environ["TMP_DIR"])
kind=pd.read_csv(tmp/"signature_kind.tsv", sep="\t")
kind_lookup=kind[["state_name", "signature_kind"]].drop_duplicates()
expr=pd.read_csv(tmp/"expression"/"all_gene_state_expression_specificity_cp10k.tsv.gz", sep="\t").merge(kind_lookup, on="state_name", how="left")
curated=expr[expr["signature_kind"].eq("curated_state")].copy()
program=expr[expr["signature_kind"].eq("program")].copy()
curated.to_csv(tmp/"expression"/"curated_state_expression.tsv.gz", sep="\t", index=False, compression="gzip")
program.to_csv(tmp/"expression"/"program_expression.tsv.gz", sep="\t", index=False, compression="gzip")
activity=pd.read_csv(tmp/"scoring"/"cell_state_activity.tsv.gz", sep="\t").merge(kind_lookup, on="state_name", how="left")
if "cell_type" not in activity.columns and "annotated_cell_type" in activity.columns:
    activity["cell_type"] = activity["annotated_cell_type"]
for cell_type, group in activity[activity["signature_kind"].eq("curated_state")].groupby("cell_type", dropna=True):
    d=tmp/"scoring"/"by_cell_type"/str(cell_type); d.mkdir(parents=True, exist_ok=True)
    group.to_csv(d/"curated_state_activity.tsv.gz", sep="\t", index=False, compression="gzip")
for cell_type, group in curated.groupby("cell_type", dropna=True):
    d=tmp/"expression"/"by_cell_type"/str(cell_type); d.mkdir(parents=True, exist_ok=True)
    group.to_csv(d/"curated_state_expression.tsv.gz", sep="\t", index=False, compression="gzip")
PY

if [[ -s "${TMP_DIR}/program_source_manifest.tsv" ]]; then
  tail -n +2 "${TMP_DIR}/program_source_manifest.tsv" | while IFS=$'\t' read -r CELL_TYPE_VALUE PROGRAM_DIR_VALUE PROGRAM_LOADINGS PROGRAM_CELL_ACTIVITY; do
    [[ -z "${CELL_TYPE_VALUE}" ]] && continue
    STATE_GMT_BY_CT="${TMP_DIR}/split_state_gmts/${CELL_TYPE_VALUE}.gmt"
    [[ -s "${STATE_GMT_BY_CT}" ]] || continue
    ACTIVITY_FILE="${TMP_DIR}/scoring/by_cell_type/${CELL_TYPE_VALUE}/curated_state_activity.tsv.gz"
    STATE_EXPR_FILE="${TMP_DIR}/expression/by_cell_type/${CELL_TYPE_VALUE}/curated_state_expression.tsv.gz"
    if [[ ! -s "${ACTIVITY_FILE}" || ! -s "${STATE_EXPR_FILE}" ]]; then
      echo "[$(date)] Skipping program-state matching for ${CELL_TYPE_VALUE}: no sampled cells or expression rows for this cell type" | tee "${LOG_DIR}/match_${CELL_TYPE_VALUE}.skipped.log"
      continue
    fi
    MATCH_OUT="${TMP_DIR}/program_state_matches/${CELL_TYPE_VALUE}"
    mkdir -p "${MATCH_OUT}"
    ACTIVITY_ARG=()
    if [[ -n "${PROGRAM_CELL_ACTIVITY}" && -s "${PROGRAM_CELL_ACTIVITY}" ]]; then
      ACTIVITY_ARG=(--program-cell-activity "${PROGRAM_CELL_ACTIVITY}")
    fi
    MATCH_NEEDED=1
    if [[ "${FORCE}" != "1" && -s "${MATCH_OUT}/program_state_heatmap_long.tsv.gz" ]]; then
      MATCH_NEEDED=$(MATCH_OUT="${MATCH_OUT}" PROGRAM_CELL_ACTIVITY="${PROGRAM_CELL_ACTIVITY}" "${PYTHON_CMD}" - <<'PYCHECK'
import os
from pathlib import Path
import pandas as pd
p = Path(os.environ["MATCH_OUT"]) / "program_state_heatmap_long.tsv.gz"
activity = os.environ.get("PROGRAM_CELL_ACTIVITY", "")
try:
    df = pd.read_csv(p, sep="\t", compression="infer", usecols=["correlation"])
    has_corr = df["correlation"].notna().any()
except Exception:
    has_corr = False
print("0" if (not activity or has_corr) else "1")
PYCHECK
)
    fi
    if [[ "${MATCH_NEEDED}" == "1" ]]; then
      echo "[$(date)] Matching programs to curated states for ${CELL_TYPE_VALUE}" | tee -a "${LOG_DIR}/pipeline.progress.log"
      "${PYTHON_CMD}" "${CELL_STATE_DE_DIR}/scripts/match_programs_to_cell_states.py" \
        --program-loadings "${PROGRAM_LOADINGS}" \
        --state-gmt "${STATE_GMT_BY_CT}" \
        "${ACTIVITY_ARG[@]}" \
        --cell-state-activity "${ACTIVITY_FILE}" \
        --state-expression "${STATE_EXPR_FILE}" \
        --metadata "${TMP_DIR}/metadata.tsv.gz" \
        --tissue "${TISSUE_ID}" \
        --cell-type "${CELL_TYPE_VALUE}" \
        --gsea-permutations 1000 \
        --api-minimal-output \
        --out-dir "${MATCH_OUT}" \
        > >(tee "${LOG_DIR}/match_${CELL_TYPE_VALUE}.stdout.log") \
        2> >(tee "${LOG_DIR}/match_${CELL_TYPE_VALUE}.stderr.log" >&2)
    else
      echo "[$(date)] Reusing existing program-state match for ${CELL_TYPE_VALUE}" | tee -a "${LOG_DIR}/pipeline.progress.log"
    fi
  done
fi


CELL_STATE_PIGEAN_EFFECTIVE="${CELL_STATE_PIGEAN}"
PROGRAM_PIGEAN_EFFECTIVE="${PROGRAM_PIGEAN}"
if [[ "${PIGEAN_ENABLE}" == "1" || "${PIGEAN_ENABLE}" == "true" ]]; then
  if [[ -z "${CELL_STATE_PIGEAN_EFFECTIVE}" ]]; then
    CELL_STATE_PIGEAN_EFFECTIVE="${TMP_DIR}/pigean/curated/combined_cell_state_pigean.tsv.gz"
    echo "[$(date)] Running/reusing curated-state PIGEAN by cell type" | tee -a "${LOG_DIR}/pipeline.progress.log"
    "${PYTHON_CMD}" "${CELL_STATE_DE_DIR}/scripts/run_api_pigean.py" \
      --gmt-dir "${TMP_DIR}/split_state_gmts" \
      --out-dir "${TMP_DIR}/pigean/curated/by_cell_type" \
      --combined-out "${CELL_STATE_PIGEAN_EFFECTIVE}" \
      --kind curated \
      --tissue "${TISSUE_ID}" \
      --dataset "${DATASET}" \
      --model "${MODEL}" \
      --python "${PYTHON_CMD}" \
      --pythonpath "${PIGEAN_PYTHONPATH}" \
      --multi-y-in "${PIGEAN_MULTI_Y_IN}" \
      --trait-blacklist-in "${PIGEAN_TRAIT_BLACKLIST}" \
      --gene-universe-in "${PIGEAN_GENE_UNIVERSE}" \
      > >(tee "${LOG_DIR}/pigean_curated.stdout.log") \
      2> >(tee "${LOG_DIR}/pigean_curated.stderr.log" >&2)
  fi
  if [[ -z "${PROGRAM_PIGEAN_EFFECTIVE}" && -d "${TMP_DIR}/split_program_gmts" ]]; then
    PROGRAM_PIGEAN_EFFECTIVE="${TMP_DIR}/pigean/program/combined_program_pigean.tsv.gz"
    echo "[$(date)] Running/reusing program PIGEAN by cell type" | tee -a "${LOG_DIR}/pipeline.progress.log"
    "${PYTHON_CMD}" "${CELL_STATE_DE_DIR}/scripts/run_api_pigean.py" \
      --gmt-dir "${TMP_DIR}/split_program_gmts" \
      --out-dir "${TMP_DIR}/pigean/program/by_cell_type" \
      --combined-out "${PROGRAM_PIGEAN_EFFECTIVE}" \
      --kind program \
      --tissue "${TISSUE_ID}" \
      --dataset "${DATASET}" \
      --model "${MODEL}" \
      --python "${PYTHON_CMD}" \
      --pythonpath "${PIGEAN_PYTHONPATH}" \
      --multi-y-in "${PIGEAN_MULTI_Y_IN}" \
      --trait-blacklist-in "${PIGEAN_TRAIT_BLACKLIST}" \
      --gene-universe-in "${PIGEAN_GENE_UNIVERSE}" \
      > >(tee "${LOG_DIR}/pigean_program.stdout.log") \
      2> >(tee "${LOG_DIR}/pigean_program.stderr.log" >&2)
  fi
fi

PORTAL_DONE="${OUT_DIR}/.api_build_complete"
if [[ "${FORCE}" == "1" || ! -s "${PORTAL_DONE}" || ! -s "${OUT_DIR}/cell_state_expression.tsv.gz" || ! -s "${OUT_DIR}/program_expression.tsv.gz" || ! -s "${OUT_DIR}/program_state_heatmap.tsv.gz" ]]; then
  echo "[$(date)] Building final API-ready data tables in ${OUT_DIR}" | tee -a "${LOG_DIR}/pipeline.progress.log"
  "${PYTHON_CMD}" "${CELL_STATE_DE_DIR}/scripts/build_portal_api_data_tables.py" \
    --out-dir "${OUT_DIR}" \
    --tissue "${TISSUE_ID}" \
    --dataset "${DATASET}" \
    --model "${MODEL}" \
    --cell-state-expression "${TMP_DIR}/expression/curated_state_expression.tsv.gz" \
    --program-expression "${TMP_DIR}/expression/program_expression.tsv.gz" \
    --cell-type-expression "${TMP_DIR}/expression/all_gene_cell_type_expression_cp10k.tsv.gz" \
    --program-loadings-manifest "${TMP_DIR}/program_source_manifest.tsv" \
    --program-match-dir "${TMP_DIR}/program_state_matches" \
    --cell-state-pigean "${CELL_STATE_PIGEAN_EFFECTIVE}" \
    --program-pigean "${PROGRAM_PIGEAN_EFFECTIVE}" \
    --default-state-weight-type "${DEFAULT_STATE_WEIGHT_TYPE}" \
    > >(tee "${LOG_DIR}/build_portal.stdout.log") \
    2> >(tee "${LOG_DIR}/build_portal.stderr.log" >&2)
  date > "${PORTAL_DONE}"
else
  echo "[$(date)] Reusing final API-ready data tables in ${OUT_DIR}" | tee -a "${LOG_DIR}/pipeline.progress.log"
fi

echo "[$(date)] Done. Final API files are in ${OUT_DIR}" | tee -a "${LOG_DIR}/pipeline.progress.log"
