#!/usr/bin/env python3
"""Summarize all-gene CP10K expression across continuous cell-state weights."""

from __future__ import annotations

import warnings
import argparse
import gzip
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread

warnings.filterwarnings("ignore", message="Mean of empty slice", category=RuntimeWarning)

from matrix_value_types import (
    VALUE_TYPES,
    expression_unit_for_value_type,
    infer_matrix_value_type,
    linearize_expression_matrix,
    load_value_type_from_metadata,
    resolve_value_type,
    specificity_label_for_value_type,
)


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def write_table(frame: pd.DataFrame, path: Path, output_format: str = "tsv") -> None:
    frame.to_csv(path, sep="\t", index=False, compression="infer")
    if output_format in {"parquet", "both"}:
        parquet_path = path.with_suffix("").with_suffix(".parquet") if path.name.endswith(".tsv.gz") else path.with_suffix(".parquet")
        try:
            frame.to_parquet(parquet_path, index=False)
        except ImportError as exc:
            raise SystemExit("Parquet output requires pyarrow or fastparquet") from exc


def read_one_column(path: Path) -> list[str]:
    with open_text(path) as handle:
        return [line.strip().split("\t")[0] for line in handle if line.strip()]


def read_features(path: Path) -> list[str]:
    genes = []
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            genes.append(parts[1] if len(parts) > 1 and parts[1] else parts[0])
    return genes


def resolve(directory: Path, names: list[str]) -> Path:
    for name in names:
        path = directory / name
        if path.exists():
            return path
    raise SystemExit(f"Could not find any of {', '.join(names)} in {directory}")


def load_10x(directory: Path) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    genes = read_features(resolve(directory, ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"]))
    cells = read_one_column(resolve(directory, ["barcodes.tsv.gz", "barcodes.tsv"]))
    mat = mmread(resolve(directory, ["matrix.mtx.gz", "matrix.mtx"])).tocsr()
    if mat.shape == (len(genes), len(cells)):
        mat = mat.T.tocsr()
    elif mat.shape != (len(cells), len(genes)):
        raise SystemExit(f"10x matrix shape {mat.shape} does not match feature/cell files")
    return mat.astype(float).tocsr(), cells, genes


def matrix_path_10x(directory: Path) -> Path:
    return resolve(directory, ["matrix.mtx.gz", "matrix.mtx"])


def iter_matrix_market_entries(path: Path):
    with open_text(path) as handle:
        dims = None
        for line in handle:
            if line.startswith("%"):
                continue
            parts = line.strip().split()
            if not parts:
                continue
            if dims is None:
                dims = tuple(int(x) for x in parts[:3])
                yield ("dims", dims)
            else:
                yield ("entry", (int(parts[0]) - 1, int(parts[1]) - 1, float(parts[2])))



def add_generic_expression_columns(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    aliases = {
        "cell_type_mean_cp10k": "cell_type_mean_expression",
        "log1p_cell_type_mean_cp10k": "log1p_cell_type_mean_expression",
        "mean_cp10k_all_context": "mean_expression_all_context",
        "log1p_mean_cp10k_all_context": "log1p_mean_expression_all_context",
        "mean_cp10k_other_context": "mean_expression_other_context",
        "log1p_mean_cp10k_other_context": "log1p_mean_expression_other_context",
        "mean_cp10k_all_parent": "mean_expression_all_parent",
        "log1p_mean_cp10k_all_parent": "log1p_mean_expression_all_parent",
        "weighted_mean_cp10k": "weighted_mean_expression",
        "log1p_weighted_mean_cp10k": "log1p_weighted_mean_expression",
        "donor_balanced_weighted_mean_cp10k": "donor_balanced_weighted_mean_expression",
        "log1p_donor_balanced_weighted_mean_cp10k": "log1p_donor_balanced_weighted_mean_expression",
    }
    for old, new in aliases.items():
        if old in frame.columns and new not in frame.columns:
            frame[new] = frame[old]
    return frame

def build_cell_type_expression_streaming_10x(
    raw_10x_dir: Path,
    cell_totals_path: str,
    metadata_path: Path,
    query_genes_path: str,
    parent_cols: list[str],
    cell_type_col: str,
    pseudocount: float,
    min_mean_for_log2fc: float,
    cell_type_filter: set[str] | None,
    expression_value_type: str,
    expression_unit: str,
) -> tuple[pd.DataFrame, int, int]:
    genes = read_features(resolve(raw_10x_dir, ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"]))
    cells = read_one_column(resolve(raw_10x_dir, ["barcodes.tsv.gz", "barcodes.tsv"]))
    metadata = read_table(metadata_path).drop_duplicates("cell_id").set_index("cell_id").reindex(cells)
    if "map_id" not in metadata.columns:
        metadata["map_id"] = "all"
    if cell_type_col not in metadata.columns and "annotated_cell_type" in metadata.columns:
        cell_type_col = "annotated_cell_type"
    if cell_type_col not in metadata.columns:
        metadata[cell_type_col] = "all"
    if "tissue" not in metadata.columns:
        metadata["tissue"] = "all"
    parent_cols = [c if c in metadata.columns else ("annotated_cell_type" if c == "cell_type" and "annotated_cell_type" in metadata.columns else c) for c in parent_cols]
    for col in parent_cols:
        if col not in metadata.columns:
            metadata[col] = "all"
    context_cols = [col for col in parent_cols if col != cell_type_col] or ["map_id", "tissue"]
    for col in context_cols:
        if col not in metadata.columns:
            metadata[col] = "all"
    query_genes = set(read_one_column(Path(query_genes_path))) if query_genes_path else set(genes)
    totals = None
    if expression_value_type == "raw_counts":
        totals_path = Path(cell_totals_path) if cell_totals_path else raw_10x_dir / "cell_total_counts.tsv.gz"
        if totals_path.exists():
            totals_frame = read_table(totals_path).set_index("cell_id")
            total_col = "total_counts" if "total_counts" in totals_frame.columns else totals_frame.columns[0]
            totals = totals_frame[total_col].reindex(cells).to_numpy(float)
        else:
            totals = np.zeros(len(cells), dtype=float)
            feature_by_cell = None
            for kind, payload in iter_matrix_market_entries(matrix_path_10x(raw_10x_dir)):
                if kind == "dims":
                    n_rows, n_cols, _ = payload
                    if n_rows == len(genes) and n_cols == len(cells):
                        feature_by_cell = True
                    elif n_rows == len(cells) and n_cols == len(genes):
                        feature_by_cell = False
                    else:
                        raise SystemExit(f"Matrix shape {(n_rows, n_cols)} does not match features/cells")
                    continue
                row, col, value = payload
                cell_idx = col if feature_by_cell else row
                totals[cell_idx] += value
        if np.any(~np.isfinite(totals)) or np.any(totals <= 0):
            raise SystemExit("Raw count totals are missing or nonpositive")
    context_key = metadata[context_cols].astype(str).agg("|".join, axis=1).to_numpy()
    cell_type_key = metadata[parent_cols].astype(str).agg("|".join, axis=1).to_numpy()
    cell_type_labels = metadata[cell_type_col].astype(str).to_numpy()
    selected_group_values = []
    for group_value in pd.unique(cell_type_key):
        first_idx = int(np.flatnonzero(cell_type_key == group_value)[0])
        if cell_type_filter and cell_type_labels[first_idx] not in cell_type_filter:
            continue
        selected_group_values.append(group_value)
    group_index = {value: i for i, value in enumerate(selected_group_values)}
    selected_contexts = {context_key[int(np.flatnonzero(cell_type_key == value)[0])] for value in selected_group_values}
    context_index = {value: i for i, value in enumerate(sorted(selected_contexts))}
    n_groups = len(group_index)
    n_contexts = len(context_index)
    n_genes = len(genes)
    group_sum = np.zeros((n_groups, n_genes), dtype=float)
    group_det = np.zeros((n_groups, n_genes), dtype=float)
    context_sum = np.zeros((n_contexts, n_genes), dtype=float)
    context_det = np.zeros((n_contexts, n_genes), dtype=float)
    group_cell_n = np.zeros(n_groups, dtype=int)
    context_cell_n = np.zeros(n_contexts, dtype=int)
    cell_to_group = np.full(len(cells), -1, dtype=int)
    cell_to_context = np.full(len(cells), -1, dtype=int)
    for i, value in enumerate(cell_type_key):
        if value in group_index:
            gi = group_index[value]
            cell_to_group[i] = gi
            group_cell_n[gi] += 1
        cv = context_key[i]
        if cv in context_index:
            ci = context_index[cv]
            cell_to_context[i] = ci
            context_cell_n[ci] += 1
    feature_by_cell = None
    for kind, payload in iter_matrix_market_entries(matrix_path_10x(raw_10x_dir)):
        if kind == "dims":
            n_rows, n_cols, _ = payload
            if n_rows == len(genes) and n_cols == len(cells):
                feature_by_cell = True
            elif n_rows == len(cells) and n_cols == len(genes):
                feature_by_cell = False
            else:
                raise SystemExit(f"Matrix shape {(n_rows, n_cols)} does not match features/cells")
            continue
        row, col, value = payload
        gene_idx, cell_idx = (row, col) if feature_by_cell else (col, row)
        if expression_value_type == "raw_counts":
            scale_value = value * 10000.0 / totals[cell_idx]
        elif expression_value_type in {"log1p_cp10k", "log1p_normalized"}:
            scale_value = float(np.expm1(value))
        else:
            scale_value = value
        gi = cell_to_group[cell_idx]
        if gi >= 0:
            group_sum[gi, gene_idx] += scale_value
            group_det[gi, gene_idx] += 1.0
        ci = cell_to_context[cell_idx]
        if ci >= 0:
            context_sum[ci, gene_idx] += scale_value
            context_det[ci, gene_idx] += 1.0
    rows = []
    value_to_context = {value: context_key[int(np.flatnonzero(cell_type_key == value)[0])] for value in selected_group_values}
    for group_value, gi in group_index.items():
        idx = np.flatnonzero(cell_type_key == group_value)
        group_meta = metadata.iloc[idx][["map_id", "tissue", cell_type_col]].astype(str).mode(dropna=False).iloc[0]
        context_value = value_to_context[group_value]
        ci = context_index[context_value]
        cell_type_mean = group_sum[gi] / max(group_cell_n[gi], 1)
        cell_type_pct = group_det[gi] / max(group_cell_n[gi], 1)
        context_mean = context_sum[ci] / max(context_cell_n[ci], 1)
        context_pct = context_det[ci] / max(context_cell_n[ci], 1)
        n_other = int(context_cell_n[ci] - group_cell_n[gi])
        if n_other > 0:
            other_mean = (context_sum[ci] - group_sum[gi]) / n_other
            other_pct = (context_det[ci] - group_det[gi]) / n_other
        else:
            other_mean = np.full_like(cell_type_mean, np.nan)
            other_pct = np.full_like(cell_type_pct, np.nan)
        with np.errstate(invalid="ignore", divide="ignore"):
            log2fc = np.log2((cell_type_mean + pseudocount) / (context_mean + pseudocount))
            log2fc_other = np.log2((cell_type_mean + pseudocount) / (other_mean + pseudocount))
        low_mean = (cell_type_mean < min_mean_for_log2fc) & (context_mean < min_mean_for_log2fc)
        log2fc[low_mean] = np.nan
        low_other = (cell_type_mean < min_mean_for_log2fc) & (other_mean < min_mean_for_log2fc)
        log2fc_other[low_other] = np.nan
        for gene, mean_value, pct_value, bg_mean, bg_pct, other_bg_mean, other_bg_pct, fc, fc_other in zip(genes, cell_type_mean, cell_type_pct, context_mean, context_pct, other_mean, other_pct, log2fc, log2fc_other):
            if gene not in query_genes:
                continue
            rows.append({
                "map_id": group_meta["map_id"],
                "tissue": group_meta["tissue"],
                "annotated_cell_type": group_meta[cell_type_col],
                "gene": gene,
                "cell_type_mean_cp10k": float(mean_value),
                "log1p_cell_type_mean_cp10k": float(np.log1p(mean_value)),
                "cell_type_pct_detected": float(pct_value),
                "mean_cp10k_all_context": float(bg_mean),
                "log1p_mean_cp10k_all_context": float(np.log1p(bg_mean)),
                "pct_detected_all_context": float(bg_pct),
                "log2fc_cell_type_vs_all_context": float(fc) if np.isfinite(fc) else np.nan,
                "mean_cp10k_other_context": float(other_bg_mean) if np.isfinite(other_bg_mean) else np.nan,
                "log1p_mean_cp10k_other_context": float(np.log1p(other_bg_mean)) if np.isfinite(other_bg_mean) else np.nan,
                "pct_detected_other_context": float(other_bg_pct) if np.isfinite(other_bg_pct) else np.nan,
                "log2fc_cell_type_vs_other_context": float(fc_other) if np.isfinite(fc_other) else np.nan,
                "n_cells_cell_type": int(group_cell_n[gi]),
                "n_cells_context": int(context_cell_n[ci]),
                "n_cells_other_context": n_other,
                "cell_type_group_key": group_value,
                "context_group_key": context_value,
                "expression_unit": expression_unit,
                "expression_result_scope": "cell_type_summary",
            })
    return pd.DataFrame(rows), len(cells), len(genes)


def bh_fdr(values: pd.Series) -> pd.Series:
    p = pd.to_numeric(values, errors="coerce")
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.notna()
    if valid.sum() == 0:
        return out
    ranked = p.loc[valid].sort_values()
    n = len(ranked)
    q = ranked.to_numpy(float) * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out.loc[ranked.index] = np.clip(q, 0, 1)
    return out


def read_gmt(path: Path) -> dict[str, set[str]]:
    if not path:
        return {}
    out: dict[str, set[str]] = {}
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                out[parts[0]] = set(parts[2:])
    return out


def weighted_mean(mat: sparse.csr_matrix, weights: np.ndarray) -> np.ndarray:
    denom = float(weights.sum())
    if denom <= 0:
        return np.full(mat.shape[1], np.nan)
    return np.asarray(weights @ mat / denom).ravel()


def weighted_pct(detected: sparse.csr_matrix, weights: np.ndarray) -> np.ndarray:
    return weighted_mean(detected, weights)


def effective_n(weights: np.ndarray) -> float:
    denom = float(np.square(weights).sum())
    return float(weights.sum() ** 2 / denom) if denom > 0 else np.nan


def sparse_square(mat: sparse.csr_matrix) -> sparse.csr_matrix:
    out = mat.copy()
    out.data = np.square(out.data)
    return out


def weighted_vs_parent_p_values(
    state_mean: np.ndarray,
    state_second: np.ndarray,
    parent_mean: np.ndarray,
    parent_second: np.ndarray,
    n_parent: int,
    eff_n: float,
) -> tuple[np.ndarray, np.ndarray]:
    state_var = np.maximum(state_second - np.square(state_mean), 0)
    parent_var = np.maximum(parent_second - np.square(parent_mean), 0)
    if not np.isfinite(eff_n) or eff_n <= 1 or n_parent <= 1:
        return np.full(len(state_mean), np.nan), np.full(len(state_mean), np.nan)
    se = np.sqrt((state_var / eff_n) + (parent_var / n_parent))
    diff = state_mean - parent_mean
    z = np.full(len(state_mean), np.nan)
    valid = np.isfinite(se) & (se > 0)
    z[valid] = diff[valid] / se[valid]
    p = np.full(len(state_mean), np.nan)
    p[valid] = 2 * stats.norm.sf(np.abs(z[valid]))
    return z, p


def expression_scope(row: pd.Series) -> str:
    if row.get("state_class") == "composite_required":
        return "component_score_summary_not_true_state"
    if row.get("threshold_status") in {"no_signal", "aucell_iqr_below_minimum"}:
        return "low_signal_flagged"
    return "state_weighted_summary"


def build_cell_type_expression(
    cp10k: sparse.csr_matrix,
    detected: sparse.csr_matrix,
    metadata: pd.DataFrame,
    cells: list[str],
    genes: list[str],
    query_genes: set[str],
    parent_cols: list[str],
    cell_type_col: str,
    pseudocount: float,
    min_mean_for_log2fc: float,
    cell_type_filter: set[str] | None = None,
    expression_unit: str = "CP10K",
) -> pd.DataFrame:
    rows = []
    if not parent_cols:
        parent_cols = ["map_id", "tissue", cell_type_col]
    context_cols = [col for col in parent_cols if col != cell_type_col]
    if not context_cols:
        context_cols = ["map_id", "tissue"]
    for col in context_cols:
        if col not in metadata.columns:
            metadata[col] = "all"
    context_key = metadata[context_cols].astype(str).agg("|".join, axis=1)
    cell_type_key = metadata[parent_cols].astype(str).agg("|".join, axis=1)
    context_values = context_key.to_numpy()
    cell_type_values = cell_type_key.to_numpy()
    for group_value in pd.unique(cell_type_values):
        idx = np.flatnonzero(cell_type_values == group_value)
        if len(idx) == 0:
            continue
        group_meta = metadata.iloc[idx][["map_id", "tissue", cell_type_col]].astype(str).mode(dropna=False).iloc[0]
        if cell_type_filter and group_meta[cell_type_col] not in cell_type_filter:
            continue
        context_value = context_values[idx[0]]
        context_idx = np.flatnonzero(context_values == context_value)
        cell_type_mean = np.asarray(cp10k[idx, :].mean(axis=0)).ravel()
        cell_type_pct = np.asarray(detected[idx, :].mean(axis=0)).ravel()
        context_mean = np.asarray(cp10k[context_idx, :].mean(axis=0)).ravel()
        context_pct = np.asarray(detected[context_idx, :].mean(axis=0)).ravel()
        other_idx = np.setdiff1d(context_idx, idx, assume_unique=False)
        if len(other_idx) > 0:
            other_mean = np.asarray(cp10k[other_idx, :].mean(axis=0)).ravel()
            other_pct = np.asarray(detected[other_idx, :].mean(axis=0)).ravel()
        else:
            other_mean = np.full_like(cell_type_mean, np.nan)
            other_pct = np.full_like(cell_type_pct, np.nan)
        with np.errstate(invalid="ignore", divide="ignore"):
            log2fc = np.log2((cell_type_mean + pseudocount) / (context_mean + pseudocount))
            log2fc_other = np.log2((cell_type_mean + pseudocount) / (other_mean + pseudocount))
        low_mean = (cell_type_mean < min_mean_for_log2fc) & (context_mean < min_mean_for_log2fc)
        log2fc[low_mean] = np.nan
        low_other = (cell_type_mean < min_mean_for_log2fc) & (other_mean < min_mean_for_log2fc)
        log2fc_other[low_other] = np.nan
        for gene, mean_value, pct_value, bg_mean, bg_pct, other_bg_mean, other_bg_pct, fc, fc_other in zip(genes, cell_type_mean, cell_type_pct, context_mean, context_pct, other_mean, other_pct, log2fc, log2fc_other):
            if gene not in query_genes:
                continue
            rows.append(
                {
                    "map_id": group_meta["map_id"],
                    "tissue": group_meta["tissue"],
                    "annotated_cell_type": group_meta[cell_type_col],
                    "gene": gene,
                    "cell_type_mean_cp10k": float(mean_value),
                    "log1p_cell_type_mean_cp10k": float(np.log1p(mean_value)),
                    "cell_type_pct_detected": float(pct_value),
                    "mean_cp10k_all_context": float(bg_mean),
                    "log1p_mean_cp10k_all_context": float(np.log1p(bg_mean)),
                    "pct_detected_all_context": float(bg_pct),
                    "log2fc_cell_type_vs_all_context": float(fc) if np.isfinite(fc) else np.nan,
                    "mean_cp10k_other_context": float(other_bg_mean) if np.isfinite(other_bg_mean) else np.nan,
                    "log1p_mean_cp10k_other_context": float(np.log1p(other_bg_mean)) if np.isfinite(other_bg_mean) else np.nan,
                    "pct_detected_other_context": float(other_bg_pct) if np.isfinite(other_bg_pct) else np.nan,
                    "log2fc_cell_type_vs_other_context": float(fc_other) if np.isfinite(fc_other) else np.nan,
                    "n_cells_cell_type": int(len(idx)),
                    "n_cells_context": int(len(context_idx)),
                    "n_cells_other_context": int(len(other_idx)),
                    "cell_type_group_key": group_value,
                    "context_group_key": context_value,
                    "expression_unit": expression_unit,
                    "expression_result_scope": "cell_type_summary",
                }
            )
    return pd.DataFrame(rows)


def api_minimal_expression(frame: pd.DataFrame, kind: str) -> pd.DataFrame:
    if frame.empty:
        base = ["gene", "tissue", "cell_type"]
        if kind == "state":
            base += ["state_name"]
        if kind == "state":
            base += ["state_weight_type"]
        base += ["weighted_mean_expression", "log10_cpk", "log2fc_weighted_vs_all_parent", "p_value"]
        return pd.DataFrame(columns=base)
    out = frame.copy()
    if "annotated_cell_type" in out.columns and "cell_type" not in out.columns:
        out["cell_type"] = out["annotated_cell_type"]
    if "weighted_mean_cp10k" in out.columns:
        cpk = pd.to_numeric(out["weighted_mean_cp10k"], errors="coerce")
    elif "weighted_mean_expression" in out.columns:
        cpk = pd.to_numeric(out["weighted_mean_expression"], errors="coerce")
    elif "cell_type_mean_cp10k" in out.columns:
        cpk = pd.to_numeric(out["cell_type_mean_cp10k"], errors="coerce")
    elif "cell_type_mean_expression" in out.columns:
        cpk = pd.to_numeric(out["cell_type_mean_expression"], errors="coerce")
    else:
        cpk = pd.Series(np.nan, index=out.index)
    out["log10_cpk"] = np.log10(cpk.fillna(0).clip(lower=0) + 1.0)
    if "weighted_mean_expression" not in out.columns:
        out["weighted_mean_expression"] = cpk
    if "log2fc_weighted_vs_all_parent" not in out.columns:
        out["log2fc_weighted_vs_all_parent"] = out.get("log2fc_cell_type_vs_other_context", np.nan)
    if "p_value" not in out.columns:
        out["p_value"] = np.nan
    cols = ["gene", "tissue", "cell_type"]
    if kind == "state":
        cols += ["state_name"]
    if kind == "state":
        cols += ["state_weight_type"]
    cols += ["weighted_mean_expression", "log10_cpk", "log2fc_weighted_vs_all_parent", "p_value"]
    for col in cols:
        if col not in out.columns:
            out[col] = "" if col in {"gene", "tissue", "cell_type", "state_name", "state_weight_type"} else np.nan
    return out[cols]

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-10x-dir", type=Path, required=True)
    ap.add_argument("--expression-value-type", choices=sorted(VALUE_TYPES), default="raw_counts", help="Value type for --raw-10x-dir; use auto to infer from matrix metadata/distribution")
    ap.add_argument("--cell-totals", default="")
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--cell-state-activity", type=Path)
    ap.add_argument("--states-gmt", default="")
    ap.add_argument("--parent-group-cols", default="map_id,tissue,cell_type")
    ap.add_argument("--donor-col", default="donor_id")
    ap.add_argument("--cell-type-col", default="cell_type")
    ap.add_argument("--pseudocount", type=float, default=0.05)
    ap.add_argument("--min-mean-for-log2fc", type=float, default=0.01)
    ap.add_argument("--query-genes", default="")
    ap.add_argument("--output-format", choices=["tsv", "parquet", "both"], default="tsv")
    ap.add_argument("--write-donor-state-expression", action=argparse.BooleanOptionalAction, default=True)
    ap.add_argument("--donor-expression-genes", choices=["query", "all", "none"], default="all")
    ap.add_argument("--cell-type-expression-only", action="store_true", help="Only write all_gene_cell_type_expression_cp10k and summary outputs.")
    ap.add_argument("--cell-type-expression-cell-types", default="", help="Optional comma-separated cell types to keep in the cell-type expression table.")
    ap.add_argument("--api-minimal-output", action="store_true", help="Write compact expression outputs with only columns needed to build portal data APIs.")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    parent_cols_arg = [c.strip() for c in args.parent_group_cols.split(",") if c.strip()]
    cell_type_filter_arg = {x.strip() for x in args.cell_type_expression_cell_types.split(",") if x.strip()} or None
    raw_probe = None
    cells = None
    genes = None
    if args.cell_type_expression_only:
        if args.expression_value_type == "auto":
            metadata_value = load_value_type_from_metadata(args.raw_10x_dir)
            if not metadata_value:
                raise SystemExit("--expression-value-type auto in --cell-type-expression-only mode requires matrix_value_type.json metadata; pass an explicit value type to avoid loading the full matrix")
            expression_value_type = metadata_value
            value_type_report = {"requested_value_type": "auto", "resolved_value_type": expression_value_type, "resolution_source": "matrix_value_type_metadata"}
        else:
            expression_value_type = args.expression_value_type
            value_type_report = {"requested_value_type": args.expression_value_type, "resolved_value_type": expression_value_type, "resolution_source": "explicit"}
    else:
        raw_probe, cells, genes = load_10x(args.raw_10x_dir)
        expression_value_type, value_type_report = resolve_value_type(args.expression_value_type, raw_probe, metadata_dir=args.raw_10x_dir, context="expression")
    expression_unit = expression_unit_for_value_type(expression_value_type)
    if args.cell_type_expression_only:
        cell_type_expression, n_cells, n_genes = build_cell_type_expression_streaming_10x(
            args.raw_10x_dir,
            args.cell_totals,
            args.metadata,
            args.query_genes,
            parent_cols_arg,
            args.cell_type_col,
            args.pseudocount,
            args.min_mean_for_log2fc,
            cell_type_filter_arg,
            expression_value_type,
            expression_unit,
        )
        cell_type_expression = add_generic_expression_columns(cell_type_expression)
        if not cell_type_expression.empty:
            cell_type_expression["expression_unit"] = expression_unit
        cell_type_out = api_minimal_expression(cell_type_expression, "cell_type") if args.api_minimal_output else cell_type_expression
        write_table(cell_type_out, args.out_dir / "all_gene_cell_type_expression_cp10k.tsv.gz", args.output_format)
        summary = {
            "raw_10x_dir": str(args.raw_10x_dir),
            "n_cells": n_cells,
            "n_genes": n_genes,
            "n_cell_type_expression_rows": int(len(cell_type_expression)),
            "cell_type_expression": str(args.out_dir / "all_gene_cell_type_expression_cp10k.tsv.gz"),
            "cell_type_filter": sorted(cell_type_filter_arg) if cell_type_filter_arg else None,
            "expression_unit": expression_unit,
            "expression_value_type": expression_value_type,
            "value_type_inference": value_type_report,
            "output_format": args.output_format,
            "api_minimal_output": args.api_minimal_output,
            "timestamp": datetime.now().isoformat(timespec="seconds"),
            "mode": "cell_type_expression_only_streaming",
        }
        (args.out_dir / "state_expression_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
        return

    counts = raw_probe
    assert cells is not None and genes is not None
    totals = None
    if expression_value_type == "raw_counts":
        totals = np.asarray(counts.sum(axis=1)).ravel()
        if args.cell_totals:
            total_frame = read_table(Path(args.cell_totals)).set_index("cell_id")
            if "total_counts" in total_frame.columns:
                totals = total_frame["total_counts"].reindex(cells).to_numpy(float)
    cp10k = linearize_expression_matrix(counts, expression_value_type, totals=totals)
    cp10k_sq = sparse_square(cp10k)
    detected = counts.copy()
    detected.data = np.ones_like(detected.data, dtype=float)
    detected = detected.astype(float).tocsr()

    metadata = read_table(args.metadata).drop_duplicates("cell_id").set_index("cell_id").reindex(cells)
    if "map_id" not in metadata.columns:
        metadata["map_id"] = "all"
    if args.cell_type_col not in metadata.columns and "annotated_cell_type" in metadata.columns:
        args.cell_type_col = "annotated_cell_type"
    if args.cell_type_col not in metadata.columns:
        metadata[args.cell_type_col] = "all"
    if "tissue" not in metadata.columns:
        metadata["tissue"] = "all"
    if args.donor_col not in metadata.columns:
        metadata[args.donor_col] = metadata.index
    parent_cols = [c.strip() for c in args.parent_group_cols.split(",") if c.strip()]
    parent_cols = [c if c in metadata.columns else ("annotated_cell_type" if c == "cell_type" and "annotated_cell_type" in metadata.columns else c) for c in parent_cols]
    for col in parent_cols:
        if col not in metadata.columns:
            metadata[col] = "all"
    parent_key = metadata[parent_cols].astype(str).agg("|".join, axis=1)

    query_genes = set(read_one_column(Path(args.query_genes))) if args.query_genes else set(genes)
    donor_gene_set = query_genes if args.donor_expression_genes == "query" else (set(genes) if args.donor_expression_genes == "all" else set())
    cell_type_filter = {x.strip() for x in args.cell_type_expression_cell_types.split(",") if x.strip()} or None

    cell_type_expression = build_cell_type_expression(
        cp10k,
        detected,
        metadata,
        cells,
        genes,
        query_genes,
        parent_cols,
        args.cell_type_col,
        args.pseudocount,
        args.min_mean_for_log2fc,
        cell_type_filter,
        expression_unit,
    )

    if args.cell_state_activity is None:
        raise SystemExit("--cell-state-activity is required unless --cell-type-expression-only is set")

    activity = read_table(args.cell_state_activity)
    activity = activity.loc[activity["state_type"].eq("biological") & activity["cell_id"].isin(cells)].copy()
    meta_for_activity = metadata[["map_id", "tissue", args.cell_type_col]].rename(columns={args.cell_type_col: "annotated_cell_type"}).reset_index()
    meta_for_activity = meta_for_activity.rename(columns={"index": "cell_id"})
    for col in ["map_id", "tissue", "annotated_cell_type"]:
        if col not in activity.columns:
            activity[col] = np.nan
    activity = activity.merge(meta_for_activity, on="cell_id", how="left", suffixes=("", "_metadata"))
    for col in ["map_id", "tissue", "annotated_cell_type"]:
        metadata_col = f"{col}_metadata"
        activity[col] = activity[col].replace("", np.nan).fillna(activity[metadata_col]).fillna("all")
        activity = activity.drop(columns=[metadata_col])

    cell_pos = pd.Series(np.arange(len(cells)), index=cells)
    markers = read_gmt(Path(args.states_gmt)) if args.states_gmt else {}

    all_parent_rows = []
    for parent in pd.unique(parent_key):
        idx = np.flatnonzero(parent_key.to_numpy() == parent)
        parent_meta = metadata.iloc[idx][["map_id", "tissue", args.cell_type_col]].astype(str).mode(dropna=False).iloc[0]
        parent_mean = np.asarray(cp10k[idx, :].mean(axis=0)).ravel()
        parent_pct = np.asarray(detected[idx, :].mean(axis=0)).ravel()
        for gene, mean_value, pct_value in zip(genes, parent_mean, parent_pct):
            if gene in query_genes:
                all_parent_rows.append(
                    {
                        "map_id": parent_meta["map_id"],
                        "tissue": parent_meta["tissue"],
                        "annotated_cell_type": parent_meta[args.cell_type_col],
                        "expression_result_scope": parent,
                        "gene": gene,
                        "mean_cp10k_all_parent": float(mean_value),
                        "log1p_mean_cp10k_all_parent": float(np.log1p(mean_value)),
                        "pct_detected_all_parent": float(pct_value),
                        "n_cells_parent": int(len(idx)),
                    }
                )
    all_parent = pd.DataFrame(all_parent_rows)

    expr_rows = []
    spec_rows = []
    donor_level_rows = []
    weight_specs = [
        ("gradient_percentile_squared", "state_activity_weight_gradient"),
        ("high_tail_percentile_90_100", "state_activity_weight_hightail"),
    ]
    state_groups: dict[tuple[str, str, str, str], pd.DataFrame] = {}
    state_meta: dict[tuple[str, str, str, str], pd.Series] = {}
    state_parent: dict[tuple[str, str, str, str], str] = {}
    group_cols = ["map_id", "tissue", "annotated_cell_type", "state_name"]
    for state_key, group in activity.groupby(group_cols, sort=False):
        map_id, tissue, annotated_cell_type, state = state_key
        positions = cell_pos.reindex(group["cell_id"]).dropna().astype(int).to_numpy()
        if len(positions):
            candidate_parents = pd.unique(parent_key.iloc[positions])
            if len(candidate_parents) != 1:
                raise SystemExit(
                    "State activity group spans multiple parent groups; rerun with parent-key columns "
                    f"or split input. Offending key: {state_key}"
                )
            parent = candidate_parents[0]
        else:
            parent = "|".join([str(tissue), str(annotated_cell_type)])
        state_groups[state_key] = group
        state_meta[state_key] = group.iloc[0]
        state_parent[state_key] = parent

    for parent in pd.unique(pd.Series(state_parent.values())):
        states_for_parent = [state_key for state_key, value in state_parent.items() if value == parent]
        parent_idx = np.flatnonzero(parent_key.to_numpy() == parent)
        if len(parent_idx) == 0:
            continue
        parent_meta = metadata.iloc[parent_idx][["map_id", "tissue", args.cell_type_col]].astype(str).mode(dropna=False).iloc[0]
        parent_mean = np.asarray(cp10k[parent_idx, :].mean(axis=0)).ravel()
        parent_second = np.asarray(cp10k_sq[parent_idx, :].mean(axis=0)).ravel()
        parent_pct = np.asarray(detected[parent_idx, :].mean(axis=0)).ravel()
        donor_ids = metadata.iloc[parent_idx][args.donor_col].astype(str).to_numpy()
        donor_levels = pd.Index(pd.unique(donor_ids))
        columns = []
        col_meta = []
        for state_key in states_for_parent:
            group = state_groups[state_key]
            state = state_key[3]
            positions = cell_pos.reindex(group["cell_id"]).dropna().astype(int).to_numpy()
            for weight_type, weight_col in weight_specs:
                weights = pd.to_numeric(group[weight_col], errors="coerce").fillna(0.0).to_numpy(float)
                keep = weights > 0
                if keep.any():
                    columns.append(sparse.csr_matrix((weights[keep], (positions[keep], np.zeros(int(keep.sum())))), shape=(len(cells), 1)))
                else:
                    columns.append(sparse.csr_matrix((len(cells), 1), dtype=float))
                col_meta.append((state_key, state, weight_type))
        if not columns:
            continue
        weight_matrix = sparse.hstack(columns, format="csr")
        denom = np.asarray(weight_matrix.sum(axis=0)).ravel()
        denom_safe = np.where(denom > 0, denom, np.nan)
        state_sums = weight_matrix.T @ cp10k
        state_detected = weight_matrix.T @ detected
        state_second_sums = weight_matrix.T @ cp10k_sq
        state_means = np.asarray(state_sums.multiply(1 / denom_safe[:, None]).todense())
        state_pcts = np.asarray(state_detected.multiply(1 / denom_safe[:, None]).todense())
        state_seconds = np.asarray(state_second_sums.multiply(1 / denom_safe[:, None]).todense())
        weight_square_sums = np.asarray(weight_matrix.power(2).sum(axis=0)).ravel()
        eff_n_by_col = np.full(len(denom), np.nan)
        valid_eff = weight_square_sums > 0
        eff_n_by_col[valid_eff] = denom[valid_eff] ** 2 / weight_square_sums[valid_eff]
        median_cells = float(pd.Series(donor_ids).value_counts().median()) if len(donor_ids) else np.nan

        donor_balanced_mean_sum = np.zeros_like(state_means)
        donor_balanced_pct_sum = np.zeros_like(state_pcts)
        donor_balanced_count = np.zeros((len(col_meta), 1), dtype=float)
        donor_counts_with_weight = np.zeros(len(col_meta), dtype=int)
        donor_weight_sums_by_col = [[] for _ in col_meta]
        for donor in donor_levels:
            didx = parent_idx[donor_ids == donor]
            if len(didx) == 0:
                continue
            donor_weights = weight_matrix[didx, :]
            donor_denoms = np.asarray(donor_weights.sum(axis=0)).ravel()
            donor_valid = donor_denoms > 0
            donor_counts_with_weight += donor_valid.astype(int)
            for col_idx, value in enumerate(donor_denoms):
                donor_weight_sums_by_col[col_idx].append(float(value))
            donor_safe = np.where(donor_valid, donor_denoms, np.nan)
            donor_mean = np.asarray((donor_weights.T @ cp10k[didx, :]).multiply(1 / donor_safe[:, None]).todense())
            donor_pct = np.asarray((donor_weights.T @ detected[didx, :]).multiply(1 / donor_safe[:, None]).todense())
            donor_mean_filled = np.where(np.isfinite(donor_mean), donor_mean, 0.0)
            donor_pct_filled = np.where(np.isfinite(donor_pct), donor_pct, 0.0)
            donor_balanced_mean_sum += donor_mean_filled
            donor_balanced_pct_sum += donor_pct_filled
            donor_balanced_count += donor_valid.astype(float)[:, None]
            for col_idx, (state_key, state, weight_type) in enumerate(col_meta):
                if not donor_valid[col_idx]:
                    continue
                for gidx, gene in enumerate(genes):
                    if gene in donor_gene_set:
                        donor_level_rows.append(
                            {
                                "donor_id": donor,
                                "map_id": state_key[0],
                                "tissue": state_key[1],
                                "annotated_cell_type": state_key[2],
                                "expression_result_scope": parent,
                                "state_name": state,
                                "state_weight_type": weight_type,
                                "gene": gene,
                                "weighted_mean_cp10k": float(donor_mean[col_idx, gidx]),
                                "log1p_weighted_mean_cp10k": float(np.log1p(donor_mean[col_idx, gidx])) if np.isfinite(donor_mean[col_idx, gidx]) else np.nan,
                                "weighted_pct_detected": float(donor_pct[col_idx, gidx]),
                                "sum_state_weight": float(donor_denoms[col_idx]),
                            }
                        )

        for col_idx, (state_key, state, weight_type) in enumerate(col_meta):
            meta = state_meta[state_key]
            state_mean = state_means[col_idx, :]
            state_second = state_seconds[col_idx, :]
            state_pct = state_pcts[col_idx, :]
            with np.errstate(invalid="ignore", divide="ignore"):
                donor_balanced = donor_balanced_mean_sum[col_idx, :] / donor_balanced_count[col_idx, 0]
                donor_balanced_pct = donor_balanced_pct_sum[col_idx, :] / donor_balanced_count[col_idx, 0]
            log2fc = np.log2((state_mean + args.pseudocount) / (parent_mean + args.pseudocount))
            low_mean = (state_mean < args.min_mean_for_log2fc) & (parent_mean < args.min_mean_for_log2fc)
            log2fc[low_mean] = np.nan
            specificity_z, specificity_p = weighted_vs_parent_p_values(
                state_mean,
                state_second,
                parent_mean,
                parent_second,
                len(parent_idx),
                float(eff_n_by_col[col_idx]),
            )
            for gidx, gene in enumerate(genes):
                if gene not in query_genes:
                    continue
                marker_set = markers.get(state, set())
                is_marker = gene in marker_set
                loo_reason = "all_gene_mode_not_recomputed" if is_marker else "not_marker"
                row = {
                    "map_id": state_key[0],
                    "tissue": state_key[1],
                    "annotated_cell_type": state_key[2],
                    "gene": gene,
                    "state_name": state,
                    "state_weight_type": weight_type,
                    "weighted_mean_cp10k": float(state_mean[gidx]),
                    "log1p_weighted_mean_cp10k": float(np.log1p(state_mean[gidx])) if np.isfinite(state_mean[gidx]) else np.nan,
                    "weighted_pct_detected": float(state_pct[gidx]),
                    "donor_balanced_weighted_mean_cp10k": float(donor_balanced[gidx]),
                    "log1p_donor_balanced_weighted_mean_cp10k": float(np.log1p(donor_balanced[gidx])) if np.isfinite(donor_balanced[gidx]) else np.nan,
                    "donor_balanced_weighted_pct_detected": float(donor_balanced_pct[gidx]),
                    "mean_cp10k_all_parent": float(parent_mean[gidx]),
                    "log1p_mean_cp10k_all_parent": float(np.log1p(parent_mean[gidx])),
                    "pct_detected_all_parent": float(parent_pct[gidx]),
                    "log2fc_weighted_vs_all_parent": float(log2fc[gidx]) if np.isfinite(log2fc[gidx]) else np.nan,
                    "p_value": float(specificity_p[gidx]) if np.isfinite(specificity_p[gidx]) else np.nan,
                    "specificity_test": specificity_label_for_value_type(expression_value_type),
                    "spearman_rho": np.nan,
                    "specificity_z": float(specificity_z[gidx]) if np.isfinite(specificity_z[gidx]) else np.nan,
                    "n_cells": int(len(parent_idx)),
                    "sum_state_weight": float(denom[col_idx]),
                    "effective_n_cells": float(eff_n_by_col[col_idx]) if np.isfinite(eff_n_by_col[col_idx]) else np.nan,
                    "n_donors_with_weight": int(donor_counts_with_weight[col_idx]),
                    "median_cells_per_donor": median_cells,
                    "state_class": meta.get("state_class", "unknown"),
                    "threshold_status": meta.get("threshold_status", "unknown"),
                    "expression_result_scope": expression_scope(meta),
                    "parent_group_key": parent,
                    "leave_one_gene_out_used": False,
                    "leave_one_gene_out_reason": loo_reason,
                    "n_markers_after_leave_one_out": max(len(marker_set) - 1, 0) if is_marker else len(marker_set),
                }
                expr_rows.append(row)
                spec_rows.append(row)

    expression = pd.DataFrame(expr_rows)
    specificity = pd.DataFrame(spec_rows)
    if not specificity.empty:
        specificity["q_value"] = bh_fdr(specificity["p_value"])
        specificity["q_by_state"] = specificity.groupby(["map_id", "tissue", "annotated_cell_type", "state_name"], group_keys=False)["p_value"].apply(bh_fdr)
        specificity["q_by_gene"] = specificity.groupby("gene", group_keys=False)["p_value"].apply(bh_fdr)
        expression = expression.merge(
            specificity[["map_id", "tissue", "annotated_cell_type", "gene", "state_name", "state_weight_type", "q_value", "q_by_state", "q_by_gene"]],
            on=["map_id", "tissue", "annotated_cell_type", "gene", "state_name", "state_weight_type"],
            how="left",
        )

    all_parent = add_generic_expression_columns(all_parent)
    cell_type_expression = add_generic_expression_columns(cell_type_expression)
    expression = add_generic_expression_columns(expression)
    specificity = add_generic_expression_columns(specificity)
    if not all_parent.empty:
        all_parent["expression_unit"] = expression_unit
    if not cell_type_expression.empty:
        cell_type_expression["expression_unit"] = expression_unit
    if not expression.empty:
        expression["expression_unit"] = expression_unit
    if not specificity.empty:
        specificity["expression_unit"] = expression_unit
    if args.api_minimal_output:
        write_table(api_minimal_expression(cell_type_expression, "cell_type"), args.out_dir / "all_gene_cell_type_expression_cp10k.tsv.gz", args.output_format)
        write_table(api_minimal_expression(expression, "state"), args.out_dir / "all_gene_state_expression_specificity_cp10k.tsv.gz", args.output_format)
    else:
        write_table(all_parent, args.out_dir / "all_gene_all_parent_cp10k.tsv.gz", args.output_format)
        write_table(cell_type_expression, args.out_dir / "all_gene_cell_type_expression_cp10k.tsv.gz", args.output_format)
        write_table(expression, args.out_dir / "all_gene_state_expression_cp10k.tsv.gz", args.output_format)
        write_table(specificity, args.out_dir / "all_gene_state_specificity_cp10k.tsv.gz", args.output_format)
        write_table(expression, args.out_dir / "all_gene_state_expression_specificity_cp10k.tsv.gz", args.output_format)
    if args.write_donor_state_expression and not args.api_minimal_output:
        donor_level = add_generic_expression_columns(pd.DataFrame(donor_level_rows))
        if not donor_level.empty:
            donor_level["expression_unit"] = expression_unit
        write_table(donor_level, args.out_dir / "donor_state_weighted_expression_cp10k.tsv.gz", args.output_format)
    summary = {
        "raw_10x_dir": str(args.raw_10x_dir),
        "cell_state_activity": str(args.cell_state_activity),
        "n_cells": len(cells),
        "n_genes": len(genes),
        "n_states": int(activity["state_name"].nunique()),
        "n_cell_type_expression_rows": int(len(cell_type_expression)),
        "cell_type_expression": str(args.out_dir / "all_gene_cell_type_expression_cp10k.tsv.gz"),
        "weight_types": [x[0] for x in weight_specs],
        "expression_unit": expression_unit,
        "expression_value_type": expression_value_type,
        "value_type_inference": value_type_report,
        "specificity_test": specificity_label_for_value_type(expression_value_type),
        "donor_state_weighted_expression": str(args.out_dir / "donor_state_weighted_expression_cp10k.tsv.gz") if args.write_donor_state_expression else None,
        "donor_expression_genes": args.donor_expression_genes,
        "output_format": args.output_format,
        "api_minimal_output": args.api_minimal_output,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    (args.out_dir / "state_expression_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
