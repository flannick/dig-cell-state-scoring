#!/usr/bin/env python3
"""Summarize all-gene CP10K expression across continuous cell-state weights."""

from __future__ import annotations

import argparse
import gzip
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def write_table(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, compression="infer")


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


def expression_scope(row: pd.Series) -> str:
    if row.get("state_class") == "composite_required":
        return "component_score_summary_not_true_state"
    if row.get("threshold_status") in {"no_signal", "aucell_iqr_below_minimum"}:
        return "low_signal_flagged"
    return "state_weighted_summary"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-10x-dir", type=Path, required=True)
    ap.add_argument("--cell-totals", default="")
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--cell-state-activity", type=Path, required=True)
    ap.add_argument("--states-gmt", default="")
    ap.add_argument("--parent-group-cols", default="tissue,cell_type")
    ap.add_argument("--donor-col", default="donor_id")
    ap.add_argument("--cell-type-col", default="cell_type")
    ap.add_argument("--pseudocount", type=float, default=0.05)
    ap.add_argument("--min-mean-for-log2fc", type=float, default=0.01)
    ap.add_argument("--query-genes", default="")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    counts, cells, genes = load_10x(args.raw_10x_dir)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    if args.cell_totals:
        total_frame = read_table(Path(args.cell_totals)).set_index("cell_id")
        if "total_counts" in total_frame.columns:
            totals = total_frame["total_counts"].reindex(cells).to_numpy(float)
    if np.any(~np.isfinite(totals)) or np.any(totals <= 0):
        raise SystemExit("Raw count totals are missing or nonpositive")
    cp10k = counts.multiply(10000.0 / totals[:, None]).tocsr()
    detected = counts.copy()
    detected.data = np.ones_like(detected.data, dtype=float)
    detected = detected.astype(float).tocsr()

    metadata = read_table(args.metadata).drop_duplicates("cell_id").set_index("cell_id").reindex(cells)
    if args.cell_type_col not in metadata.columns and "annotated_cell_type" in metadata.columns:
        args.cell_type_col = "annotated_cell_type"
    if args.donor_col not in metadata.columns:
        metadata[args.donor_col] = metadata.index
    parent_cols = [c.strip() for c in args.parent_group_cols.split(",") if c.strip()]
    parent_cols = [c if c in metadata.columns else ("annotated_cell_type" if c == "cell_type" and "annotated_cell_type" in metadata.columns else c) for c in parent_cols]
    for col in parent_cols:
        if col not in metadata.columns:
            metadata[col] = "all"
    parent_key = metadata[parent_cols].astype(str).agg("|".join, axis=1)

    activity = read_table(args.cell_state_activity)
    activity = activity.loc[activity["state_type"].eq("biological") & activity["cell_id"].isin(cells)].copy()
    cell_pos = pd.Series(np.arange(len(cells)), index=cells)
    markers = read_gmt(Path(args.states_gmt)) if args.states_gmt else {}
    query_genes = set(read_one_column(Path(args.query_genes))) if args.query_genes else set(genes)

    all_parent_rows = []
    for parent in pd.unique(parent_key):
        idx = np.flatnonzero(parent_key.to_numpy() == parent)
        parent_mean = np.asarray(cp10k[idx, :].mean(axis=0)).ravel()
        parent_pct = np.asarray(detected[idx, :].mean(axis=0)).ravel()
        for gene, mean_value, pct_value in zip(genes, parent_mean, parent_pct):
            if gene in query_genes:
                all_parent_rows.append(
                    {
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
    weight_specs = [
        ("gradient_percentile_squared", "state_activity_weight_gradient"),
        ("high_tail_percentile_90_100", "state_activity_weight_hightail"),
    ]
    for state, group in activity.groupby("state_name", sort=False):
        meta = group.iloc[0]
        positions = cell_pos.reindex(group["cell_id"]).dropna().astype(int).to_numpy()
        parent = parent_key.iloc[positions].mode().iloc[0] if len(positions) else "all"
        parent_idx = np.flatnonzero(parent_key.to_numpy() == parent)
        parent_mean = np.asarray(cp10k[parent_idx, :].mean(axis=0)).ravel()
        parent_pct = np.asarray(detected[parent_idx, :].mean(axis=0)).ravel()
        donor_ids = metadata.iloc[parent_idx][args.donor_col].astype(str).to_numpy()
        donor_levels = pd.Index(pd.unique(donor_ids))
        parent_donor_means = []
        for donor in donor_levels:
            didx = parent_idx[donor_ids == donor]
            parent_donor_means.append(np.asarray(cp10k[didx, :].mean(axis=0)).ravel())
        parent_donor_means = np.vstack(parent_donor_means) if len(parent_donor_means) else np.empty((0, len(genes)))

        for weight_type, weight_col in weight_specs:
            weights = np.zeros(len(cells), dtype=float)
            weights[positions] = pd.to_numeric(group[weight_col], errors="coerce").fillna(0.0).to_numpy(float)
            state_mean = weighted_mean(cp10k, weights)
            state_pct = weighted_pct(detected, weights)
            donor_weighted = []
            donor_weighted_pct = []
            donor_weight_sums = []
            for donor in donor_levels:
                didx = parent_idx[donor_ids == donor]
                dw = weights[didx]
                donor_weight_sums.append(float(dw.sum()))
                donor_weighted.append(weighted_mean(cp10k[didx, :], dw))
                donor_weighted_pct.append(weighted_pct(detected[didx, :], dw))
            donor_weighted = np.vstack(donor_weighted) if len(donor_weighted) else np.empty((0, len(genes)))
            donor_weighted_pct = np.vstack(donor_weighted_pct) if len(donor_weighted_pct) else np.empty((0, len(genes)))
            donor_balanced = np.nanmean(donor_weighted, axis=0) if donor_weighted.size else np.full(len(genes), np.nan)
            donor_balanced_pct = np.nanmean(donor_weighted_pct, axis=0) if donor_weighted_pct.size else np.full(len(genes), np.nan)
            log2fc = np.log2((state_mean + args.pseudocount) / (parent_mean + args.pseudocount))
            low_mean = (state_mean < args.min_mean_for_log2fc) & (parent_mean < args.min_mean_for_log2fc)
            log2fc[low_mean] = np.nan
            n_donors_with_weight = int(np.sum(np.asarray(donor_weight_sums) > 0))
            median_cells = float(pd.Series(donor_ids).value_counts().median()) if len(donor_ids) else np.nan

            for gidx, gene in enumerate(genes):
                if gene not in query_genes:
                    continue
                rho, pval = (np.nan, np.nan)
                if np.isfinite(weights).sum() > 2 and np.unique(weights).size > 1:
                    rho, pval = stats.spearmanr(np.log1p(cp10k[:, gidx].toarray().ravel()), weights, nan_policy="omit")
                marker_set = markers.get(state, set())
                is_marker = gene in marker_set
                loo_reason = "not_marker"
                if is_marker:
                    loo_reason = "all_gene_mode_not_recomputed"
                row = {
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
                    "p_value": float(pval) if np.isfinite(pval) else np.nan,
                    "specificity_test": "spearman_log1p_cp10k_vs_state_weight",
                    "spearman_rho": float(rho) if np.isfinite(rho) else np.nan,
                    "n_cells": int(len(parent_idx)),
                    "sum_state_weight": float(weights.sum()),
                    "effective_n_cells": effective_n(weights),
                    "n_donors_with_weight": n_donors_with_weight,
                    "median_cells_per_donor": median_cells,
                    "state_class": meta.get("state_class", "unknown"),
                    "threshold_status": meta.get("threshold_status", "unknown"),
                    "expression_result_scope": expression_scope(meta),
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
        specificity["q_by_state"] = specificity.groupby("state_name", group_keys=False)["p_value"].apply(bh_fdr)
        specificity["q_by_gene"] = specificity.groupby("gene", group_keys=False)["p_value"].apply(bh_fdr)
        expression = expression.merge(
            specificity[["gene", "state_name", "state_weight_type", "q_value", "q_by_state", "q_by_gene"]],
            on=["gene", "state_name", "state_weight_type"],
            how="left",
        )

    write_table(all_parent, args.out_dir / "all_gene_all_parent_cp10k.tsv.gz")
    write_table(expression, args.out_dir / "all_gene_state_expression_cp10k.tsv.gz")
    write_table(specificity, args.out_dir / "all_gene_state_specificity_cp10k.tsv.gz")
    write_table(expression, args.out_dir / "all_gene_state_expression_specificity_cp10k.tsv.gz")
    summary = {
        "raw_10x_dir": str(args.raw_10x_dir),
        "cell_state_activity": str(args.cell_state_activity),
        "n_cells": len(cells),
        "n_genes": len(genes),
        "n_states": int(activity["state_name"].nunique()),
        "weight_types": [x[0] for x in weight_specs],
        "expression_unit": "CP10K",
        "specificity_test": "spearman_log1p_cp10k_vs_state_weight",
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    (args.out_dir / "state_expression_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
