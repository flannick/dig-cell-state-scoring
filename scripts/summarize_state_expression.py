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


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-10x-dir", type=Path, required=True)
    ap.add_argument("--cell-totals", default="")
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--cell-state-activity", type=Path, required=True)
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
    query_genes = set(read_one_column(Path(args.query_genes))) if args.query_genes else set(genes)
    donor_gene_set = query_genes if args.donor_expression_genes == "query" else (set(genes) if args.donor_expression_genes == "all" else set())

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

    for parent in pd.unique(list(state_parent.values())):
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
                    "specificity_test": "weighted_vs_parent_mean_cp10k_normal_approximation_screening",
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

    write_table(all_parent, args.out_dir / "all_gene_all_parent_cp10k.tsv.gz", args.output_format)
    write_table(expression, args.out_dir / "all_gene_state_expression_cp10k.tsv.gz", args.output_format)
    write_table(specificity, args.out_dir / "all_gene_state_specificity_cp10k.tsv.gz", args.output_format)
    write_table(expression, args.out_dir / "all_gene_state_expression_specificity_cp10k.tsv.gz", args.output_format)
    if args.write_donor_state_expression:
        write_table(pd.DataFrame(donor_level_rows), args.out_dir / "donor_state_weighted_expression_cp10k.tsv.gz", args.output_format)
    summary = {
        "raw_10x_dir": str(args.raw_10x_dir),
        "cell_state_activity": str(args.cell_state_activity),
        "n_cells": len(cells),
        "n_genes": len(genes),
        "n_states": int(activity["state_name"].nunique()),
        "weight_types": [x[0] for x in weight_specs],
        "expression_unit": "CP10K",
        "specificity_test": "weighted_vs_parent_mean_cp10k_normal_approximation_screening",
        "donor_state_weighted_expression": str(args.out_dir / "donor_state_weighted_expression_cp10k.tsv.gz") if args.write_donor_state_expression else None,
        "donor_expression_genes": args.donor_expression_genes,
        "output_format": args.output_format,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    (args.out_dir / "state_expression_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
