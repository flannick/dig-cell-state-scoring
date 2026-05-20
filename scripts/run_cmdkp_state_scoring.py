#!/usr/bin/env python3
"""Run the general CMDKP cell-state scoring and calling workflow.

The runner consumes an exported cell x gene expression table plus cell metadata.
It scores biological state GMTs and the auxiliary bad-cell QC GMT with a local
UCell-style rank statistic, calibrates state thresholds, and writes the standard
CMDKP output bundle.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from datetime import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import yaml
from scipy import stats
from scipy.ndimage import gaussian_filter1d


PROCESS_TOKENS = {
    "upr",
    "er_stress",
    "interferon",
    "ifn",
    "mhc",
    "oxidative",
    "hypoxia",
    "proliferation",
    "proliferative",
    "cell_cycle",
    "cycling",
    "fibrosis",
    "remodeling",
}
COMPOSITE_TOKENS = {
    "dedifferentiation",
    "dedifferentiated",
    "low_identity",
    "disallowed",
    "senescence",
    "senescence_like",
    "doublet",
    "contamination",
}
HEMOGLOBIN_PREFIXES = ("HBA", "HBB", "HBD", "HBG", "HBM", "HBQ")


@dataclass
class GeneSet:
    name: str
    description: str
    genes: list[str]
    meta: dict[str, str]


def norm_name(value: str) -> str:
    value = str(value).strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return re.sub(r"_+", "_", value).strip("_")


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def write_table(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, compression="infer")


def empty_table(columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame({column: pd.Series(dtype="object") for column in columns})


def output_path(args: argparse.Namespace, attr: str, default_name: str) -> Path:
    explicit = getattr(args, attr)
    if explicit:
        return Path(explicit)
    return Path(args.out_dir) / default_name


def read_expression_input(path: str) -> pd.DataFrame:
    frame = read_table(path)
    if {"cell_id", "gene", "expression"}.issubset(frame.columns):
        return frame[["cell_id", "gene", "expression"]].copy()
    if "cell_id" not in frame.columns:
        first = frame.columns[0]
        frame = frame.rename(columns={first: "cell_id"})
    long = frame.melt(id_vars=["cell_id"], var_name="gene", value_name="expression")
    long = long.loc[pd.to_numeric(long["expression"], errors="coerce").fillna(0) != 0].copy()
    return long


def load_threshold_yaml(path: str) -> dict[str, float]:
    if not path:
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        obj = yaml.safe_load(handle) or {}
    if isinstance(obj, dict) and "states" in obj:
        obj = obj["states"]
    out: dict[str, float] = {}
    if isinstance(obj, dict):
        for key, value in obj.items():
            if isinstance(value, dict):
                value = value.get("threshold", value.get("state_probability", value.get("probability")))
            if value is not None:
                out[str(key)] = float(value)
    return out


def read_gene_list(path: str) -> list[str]:
    if not path:
        return []
    genes: list[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            genes.append(value.split()[0])
    return list(dict.fromkeys(genes))


def resolve_query_genes(args: argparse.Namespace, matrix: pd.DataFrame) -> list[str]:
    requested = read_gene_list(args.query_genes) if args.query_genes else []
    if args.query_gene:
        requested.extend(g.strip() for g in args.query_gene.split(",") if g.strip())
    if not requested:
        return list(matrix.columns)
    requested = list(dict.fromkeys(requested))
    present = [gene for gene in requested if gene in matrix.columns]
    missing = [gene for gene in requested if gene not in matrix.columns]
    if missing:
        print(f"Warning: {len(missing)} query genes were not found in the expression matrix", file=sys.stderr)
    if not present:
        raise SystemExit("No query genes were found in the expression matrix")
    return present


def bh_fdr(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce")
    out = pd.Series(np.nan, index=p.index)
    valid = p.notna()
    if valid.sum() == 0:
        return out
    order = p[valid].sort_values().index
    ranked = p.loc[order].to_numpy()
    n = len(ranked)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out.loc[order] = np.clip(q, 0, 1)
    return out


def parse_description(desc: str) -> dict[str, str]:
    meta: dict[str, str] = {}
    for part in str(desc).split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            meta[key.strip()] = value.strip()
    return meta


def read_gmt(path: str, regex: str = "") -> list[GeneSet]:
    pattern = re.compile(regex) if regex else None
    sets: list[GeneSet] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            if pattern and not pattern.search(name):
                continue
            genes = list(dict.fromkeys(g for g in parts[2:] if g))
            sets.append(GeneSet(name=name, description=parts[1], genes=genes, meta=parse_description(parts[1])))
    if not sets:
        raise SystemExit(f"No GMT rows found in {path}")
    return sets


def harmonize_expression(
    expression: pd.DataFrame,
    expression_kind: str,
    gene_map_path: str,
    duplicate_method: str,
) -> tuple[pd.DataFrame, dict[str, str]]:
    required = {"cell_id", "gene", "expression"}
    missing = sorted(required - set(expression.columns))
    if missing:
        raise SystemExit(f"Expression table is missing required column(s): {', '.join(missing)}")
    expr = expression[["cell_id", "gene", "expression"]].copy()
    expr["expression"] = pd.to_numeric(expr["expression"], errors="coerce").fillna(0.0)
    mapping_info = {"gene_id_type": "assumed_hgnc_symbol", "duplicate_collapse_method": "none"}

    if gene_map_path:
        gene_map = read_table(gene_map_path)
        for col in ["gene_id", "gene_symbol"]:
            if col not in gene_map.columns:
                raise SystemExit(f"Gene map is missing required column: {col}")
        gene_map = gene_map[["gene_id", "gene_symbol"]].dropna().drop_duplicates("gene_id")
        expr = expr.merge(gene_map, left_on="gene", right_on="gene_id", how="left")
        expr["gene"] = expr["gene_symbol"].fillna(expr["gene"])
        expr = expr.drop(columns=["gene_id", "gene_symbol"])
        mapping_info["gene_id_type"] = "mapped_to_hgnc_symbol"

    expr["gene"] = expr["gene"].astype(str)
    if duplicate_method == "auto":
        duplicate_method = "sum" if expression_kind == "raw_counts" else "max"
    mapping_info["duplicate_collapse_method"] = duplicate_method
    if duplicate_method == "sum":
        expr = expr.groupby(["cell_id", "gene"], as_index=False)["expression"].sum()
    elif duplicate_method == "mean":
        expr = expr.groupby(["cell_id", "gene"], as_index=False)["expression"].mean()
    elif duplicate_method == "max":
        expr = expr.groupby(["cell_id", "gene"], as_index=False)["expression"].max()
    else:
        raise SystemExit("--duplicate-collapse must be auto, sum, mean, or max")
    return expr, mapping_info


def expression_matrix(expression: pd.DataFrame, cells: pd.Series) -> pd.DataFrame:
    matrix = expression.pivot(index="cell_id", columns="gene", values="expression")
    matrix = matrix.reindex(cells).fillna(0.0)
    return matrix.astype(float)


def rank_matrix_for_ucell(matrix: pd.DataFrame, max_rank: int) -> pd.DataFrame:
    ranks = matrix.rank(axis=1, ascending=False, method="average")
    ranks = ranks.where(ranks <= max_rank, max_rank + 1)
    return ranks


def ucell_score_from_ranks(ranks: pd.DataFrame, genes: list[str], max_rank: int) -> pd.Series:
    present = [g for g in genes if g in ranks.columns]
    n = len(present)
    if n == 0:
        return pd.Series(np.nan, index=ranks.index)
    max_u = n * max_rank - (n * (n + 1)) / 2
    if max_u <= 0:
        return pd.Series(np.nan, index=ranks.index)
    rank_sum = ranks[present].sum(axis=1)
    u_stat = rank_sum - (n * (n + 1)) / 2
    return (1 - (u_stat / max_u)).clip(0, 1)


def marker_info(gene_set: GeneSet, matrix_genes: Iterable[str]) -> dict[str, object]:
    genes = set(matrix_genes)
    present = [g for g in gene_set.genes if g in genes]
    missing = [g for g in gene_set.genes if g not in genes]
    return {
        "markers_requested": ";".join(gene_set.genes),
        "markers_present": ";".join(present),
        "markers_missing": ";".join(missing),
        "n_markers_requested": len(gene_set.genes),
        "n_markers_present": len(present),
        "marker_coverage_fraction": len(present) / len(gene_set.genes) if gene_set.genes else np.nan,
        "present_list": present,
    }


def infer_state_scope(state_name: str, tissues: Iterable[str], cell_types: Iterable[str]) -> tuple[str, str, str]:
    norm_state = norm_name(state_name)
    tissue_match = ""
    for tissue in sorted({norm_name(t) for t in tissues}, key=len, reverse=True):
        if norm_state.startswith(tissue + "_"):
            tissue_match = tissue
            remainder = norm_state[len(tissue) + 1 :]
            break
    else:
        parts = norm_state.split("_", 1)
        tissue_match = parts[0]
        remainder = parts[1] if len(parts) > 1 else ""

    cell_match = ""
    state_suffix = remainder
    for cell_type in sorted({norm_name(c) for c in cell_types}, key=len, reverse=True):
        if remainder.startswith(cell_type + "_"):
            cell_match = cell_type
            state_suffix = remainder[len(cell_type) + 1 :]
            break
    return tissue_match, cell_match, state_suffix


def state_kind(state_name: str) -> str:
    name = norm_name(state_name)
    if any(token in name for token in PROCESS_TOKENS):
        return "process"
    return "biological"


def is_composite_state(state_name: str) -> bool:
    name = norm_name(state_name)
    return any(token in name for token in COMPOSITE_TOKENS)


def score_biological_states(
    matrix: pd.DataFrame,
    ranks: pd.DataFrame,
    metadata: pd.DataFrame,
    gene_sets: list[GeneSet],
    args: argparse.Namespace,
) -> pd.DataFrame:
    rows = []
    tissues = metadata[args.tissue_col].unique()
    cell_types = metadata[args.cell_type_col].unique()
    group_keys = [args.map_id_col, args.tissue_col, args.cell_type_col]
    for gene_set in gene_sets:
        scope_tissue, scope_cell_type, _ = infer_state_scope(gene_set.name, tissues, cell_types)
        info = marker_info(gene_set, matrix.columns)
        score = ucell_score_from_ranks(ranks, info["present_list"], args.max_rank)
        score_frame = metadata[["cell_id", args.map_id_col, args.tissue_col, args.cell_type_col]].copy()
        score_frame["state_name"] = gene_set.name
        score_frame["ucell_score"] = score.reindex(score_frame["cell_id"]).to_numpy()
        score_frame["scope_tissue"] = scope_tissue
        score_frame["scope_cell_type"] = scope_cell_type
        score_frame["state_kind"] = state_kind(gene_set.name)
        score_frame["is_composite_state"] = is_composite_state(gene_set.name)
        for key, value in info.items():
            if key != "present_list":
                score_frame[key] = value

        relevant = (metadata[args.tissue_col].map(norm_name) == scope_tissue)
        if scope_cell_type:
            relevant &= metadata[args.cell_type_col].map(norm_name) == scope_cell_type
        score_frame = score_frame.loc[relevant.to_numpy()].copy()
        if score_frame.empty:
            continue
        score_frame["score_percentile_within_calibration_group"] = (
            score_frame.groupby(group_keys + ["state_name"])["ucell_score"].rank(pct=True, method="average")
        )
        rows.append(score_frame)
    if not rows:
        raise SystemExit("No biological state scores were relevant to the provided metadata")
    out = pd.concat(rows, ignore_index=True)
    return out.rename(
        columns={
            args.map_id_col: "map_id",
            args.tissue_col: "tissue",
            args.cell_type_col: "annotated_cell_type",
        }
    )


def score_qc_signatures(
    matrix: pd.DataFrame,
    ranks: pd.DataFrame,
    metadata: pd.DataFrame,
    qc_sets: list[GeneSet],
    args: argparse.Namespace,
) -> pd.DataFrame:
    rows = []
    for gene_set in qc_sets:
        info = marker_info(gene_set, matrix.columns)
        score = ucell_score_from_ranks(ranks, info["present_list"], args.max_rank)
        frame = metadata[["cell_id", args.map_id_col, args.tissue_col, args.cell_type_col, args.sample_col]].copy()
        frame["qc_signature_name"] = gene_set.name
        frame["ucell_score"] = score.reindex(frame["cell_id"]).to_numpy()
        frame["qc_tier"] = gene_set.meta.get("tier", "")
        frame["qc_category"] = gene_set.meta.get("category", "")
        frame["score_percentile_within_sample"] = (
            frame.groupby([args.map_id_col, args.sample_col, "qc_signature_name"])["ucell_score"].rank(pct=True, method="average")
        )
        for key, value in info.items():
            if key != "present_list":
                frame[key] = value
        rows.append(frame)
    out = pd.concat(rows, ignore_index=True)
    return out.rename(
        columns={
            args.map_id_col: "map_id",
            args.tissue_col: "tissue",
            args.cell_type_col: "annotated_cell_type",
        }
    )


def mad(series: pd.Series) -> float:
    med = series.median()
    return float((series - med).abs().median())


def adaptive_high_flag(values: pd.Series) -> pd.Series:
    med = values.median()
    spread = mad(values)
    return values > med + 3 * spread if spread > 0 else pd.Series(False, index=values.index)


def adaptive_low_flag(values: pd.Series) -> pd.Series:
    med = values.median()
    spread = mad(values)
    return values < med - 3 * spread if spread > 0 else pd.Series(False, index=values.index)


def qc_metrics(matrix: pd.DataFrame, metadata: pd.DataFrame, qc_scores: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    genes = pd.Series(matrix.columns, index=matrix.columns)
    total_counts = matrix.sum(axis=1)
    n_genes = (matrix > 0).sum(axis=1)
    mito = matrix.loc[:, genes.str.startswith("MT-").to_numpy()].sum(axis=1) if genes.str.startswith("MT-").any() else 0
    ribo = matrix.loc[:, genes.str.match(r"RP[LS]").to_numpy()].sum(axis=1) if genes.str.match(r"RP[LS]").any() else 0
    hemo_mask = genes.str.startswith(HEMOGLOBIN_PREFIXES)
    hemo = matrix.loc[:, hemo_mask.to_numpy()].sum(axis=1) if hemo_mask.any() else 0
    malat = matrix["MALAT1"] if "MALAT1" in matrix.columns else 0
    out = metadata[["cell_id", args.map_id_col, args.tissue_col, args.cell_type_col, args.sample_col]].copy()
    out["total_counts"] = total_counts.reindex(out["cell_id"]).to_numpy()
    out["n_genes_detected"] = n_genes.reindex(out["cell_id"]).to_numpy()
    denom = out["total_counts"].replace(0, np.nan)
    out["percent_mitochondrial"] = 100 * np.asarray(mito.reindex(out["cell_id"]) if hasattr(mito, "reindex") else 0) / denom
    out["percent_ribosomal"] = 100 * np.asarray(ribo.reindex(out["cell_id"]) if hasattr(ribo, "reindex") else 0) / denom
    out["percent_hemoglobin"] = 100 * np.asarray(hemo.reindex(out["cell_id"]) if hasattr(hemo, "reindex") else 0) / denom
    out["percent_malat1"] = 100 * np.asarray(malat.reindex(out["cell_id"]) if hasattr(malat, "reindex") else 0) / denom
    if args.doublet_col and args.doublet_col in metadata.columns:
        out["doublet_score"] = pd.to_numeric(metadata[args.doublet_col], errors="coerce")
    else:
        out["doublet_score"] = np.nan

    out["hard_exclusion_flag"] = False
    out["review_flag"] = False
    reasons = {cell: [] for cell in out["cell_id"]}
    for _, group in out.groupby([args.map_id_col, args.sample_col], sort=False):
        idx = group.index
        high_total = adaptive_high_flag(group["total_counts"])
        low_total = adaptive_low_flag(group["total_counts"])
        high_genes = adaptive_high_flag(group["n_genes_detected"])
        low_genes = adaptive_low_flag(group["n_genes_detected"])
        high_mito = adaptive_high_flag(group["percent_mitochondrial"].fillna(0))
        high_ribo = adaptive_high_flag(group["percent_ribosomal"].fillna(0))
        high_hemo = adaptive_high_flag(group["percent_hemoglobin"].fillna(0))
        out.loc[idx[high_mito | high_hemo], "hard_exclusion_flag"] = True
        out.loc[idx[high_total | low_total | high_genes | low_genes | high_ribo], "review_flag"] = True
        for flag, label in [
            (high_mito, "high_percent_mitochondrial"),
            (high_hemo, "high_percent_hemoglobin"),
            (high_ribo, "high_percent_ribosomal"),
            (high_total, "high_total_counts"),
            (low_total, "low_total_counts"),
            (high_genes, "high_n_genes_detected"),
            (low_genes, "low_n_genes_detected"),
        ]:
            for cell in group.loc[flag, "cell_id"]:
                reasons[cell].append(label)

    high_qc = qc_scores.loc[qc_scores["score_percentile_within_sample"] >= args.qc_extreme_percentile].copy()
    hard_qc = high_qc["qc_tier"].str.contains("hard_exclude", na=False)
    review_qc = high_qc["qc_tier"].str.contains("review", na=False)
    hard_cells = set(high_qc.loc[hard_qc, "cell_id"])
    review_cells = set(high_qc.loc[review_qc, "cell_id"])
    out.loc[out["cell_id"].isin(hard_cells), "hard_exclusion_flag"] = True
    out.loc[out["cell_id"].isin(review_cells), "review_flag"] = True
    for _, row in high_qc.loc[hard_qc | review_qc].iterrows():
        reasons[row["cell_id"]].append(row["qc_signature_name"])

    identity = high_qc.loc[high_qc["qc_signature_name"].str.contains("identity", case=False, na=False)]
    ambient = high_qc.loc[high_qc["qc_signature_name"].str.contains("ambient", case=False, na=False)]
    top_identity = identity.sort_values("ucell_score").groupby("cell_id").tail(1).set_index("cell_id")["qc_signature_name"]
    top_ambient = ambient.sort_values("ucell_score").groupby("cell_id").tail(1).set_index("cell_id")["qc_signature_name"]
    out["bad_cell_reason"] = out["cell_id"].map(lambda c: ";".join(sorted(set(reasons[c]))) if reasons[c] else "none")
    out["top_offtarget_identity_signature"] = out["cell_id"].map(top_identity).fillna("none")
    out["top_ambient_signature"] = out["cell_id"].map(top_ambient).fillna("none")
    out["parent_identity_score"] = np.nan
    return out.rename(
        columns={
            args.map_id_col: "map_id",
            args.tissue_col: "tissue",
            args.cell_type_col: "annotated_cell_type",
        }
    )


def random_gene_sets(present: list[str], genes: list[str], bin_key: pd.Series, n_sets: int) -> list[list[str]]:
    out = []
    avoid = set(present)
    for _ in range(n_sets):
        sampled = []
        for gene in present:
            same_bin = [g for g in genes if bin_key.get(g) == bin_key.get(gene) and g not in avoid]
            pool = same_bin or [g for g in genes if g not in avoid] or genes
            sampled.append(str(np.random.choice(pool)))
        out.append(sampled)
    return out


def stable_seed(*parts: object) -> int:
    key = "|".join(str(part) for part in parts)
    digest = hashlib.md5(key.encode("utf-8")).hexdigest()
    return int(digest[:8], 16)


def select_null_calibration_cells(cells: list[str], state_name: str, args: argparse.Namespace) -> list[str]:
    unique_cells = sorted(dict.fromkeys(cells))
    max_cells = int(args.null_max_cells)
    if max_cells <= 0 or len(unique_cells) <= max_cells:
        return unique_cells
    rng = np.random.default_rng(args.random_seed + stable_seed("null_calibration", state_name))
    selected = rng.choice(unique_cells, size=max_cells, replace=False)
    return sorted(str(cell) for cell in selected)


def state_seed(args: argparse.Namespace, label: str, state_name: str) -> int:
    return int((args.random_seed + stable_seed(label, state_name)) % (2**32 - 1))


def null_scores_for_state(
    matrix: pd.DataFrame,
    ranks: pd.DataFrame,
    cells: list[str],
    present: list[str],
    args: argparse.Namespace,
    state_name: str,
) -> tuple[list[float], int]:
    if not present:
        return [], 0
    calibration_cells = select_null_calibration_cells(cells, state_name, args)
    sub_matrix = matrix.loc[calibration_cells]
    sub_ranks = ranks.loc[calibration_cells]
    means = sub_matrix.mean(axis=0)
    detected = (sub_matrix > 0).mean(axis=0)
    bin_key = make_bins(means, args.expression_bins).astype(str) + ":" + make_bins(detected, args.detection_bins).astype(str)
    np.random.seed(state_seed(args, "null_sets", state_name))
    random_sets = random_gene_sets(present, list(matrix.columns), bin_key, args.null_n)
    values: list[float] = []
    for genes in random_sets:
        values.extend(ucell_score_from_ranks(sub_ranks, genes, args.max_rank).dropna().tolist())
    return values, len(calibration_cells)


def probability_from_null(real_scores: pd.Series, null_scores: list[float], bins: int = 101) -> tuple[pd.Series, str]:
    real = pd.to_numeric(real_scores, errors="coerce")
    null = pd.Series(null_scores, dtype=float).dropna()
    if real.notna().sum() < 3 or len(null) < 10 or real.nunique(dropna=True) < 2:
        if len(null) == 0:
            return pd.Series(0.0, index=real_scores.index), "empirical_null_cdf_fallback"
        probs = real.apply(lambda x: np.nan if pd.isna(x) else float((null < x).mean()))
        probs = probs.fillna(0).clip(0, 1)
        return monotone_probabilities(real, probs), "empirical_null_cdf_fallback"

    try:
        edges = np.linspace(0, 1, bins)
        centers = (edges[:-1] + edges[1:]) / 2
        null_hist, _ = np.histogram(null.clip(0, 1), bins=edges, density=True)
        obs_hist, _ = np.histogram(real.dropna().clip(0, 1), bins=edges, density=True)
        null_density = gaussian_filter1d(null_hist.astype(float), sigma=1.0) + 1e-8
        obs_density = gaussian_filter1d(obs_hist.astype(float), sigma=1.0) + 1e-8
        f0 = np.interp(real.fillna(0).clip(0, 1), centers, null_density)
        fobs = np.interp(real.fillna(0).clip(0, 1), centers, obs_density)
        lfdr = np.minimum(1.0, f0 / fobs)
        probs = pd.Series(1.0 - lfdr, index=real.index).where(real.notna(), np.nan).fillna(0).clip(0, 1)
        return monotone_probabilities(real, probs), "smoothed_histogram_lfdr"
    except Exception:
        probs = real.apply(lambda x: np.nan if pd.isna(x) else float((null < x).mean()))
        probs = probs.fillna(0).clip(0, 1)
        return monotone_probabilities(real, probs), "empirical_null_cdf_fallback"


def monotone_probabilities(scores: pd.Series, probs: pd.Series) -> pd.Series:
    frame = pd.DataFrame({"score": scores, "prob": probs}).sort_values("score", kind="mergesort")
    frame["prob"] = frame["prob"].cummax()
    return frame.sort_index()["prob"].clip(0, 1)


def make_bins(values: pd.Series, n_bins: int) -> pd.Series:
    if values.nunique(dropna=True) <= 1:
        return pd.Series(1, index=values.index)
    try:
        return pd.qcut(values.rank(method="first"), q=min(n_bins, values.nunique()), labels=False, duplicates="drop").astype(int)
    except ValueError:
        return pd.Series(1, index=values.index)


def calibrate_thresholds(
    scores: pd.DataFrame,
    matrix: pd.DataFrame,
    ranks: pd.DataFrame,
    biological_sets: dict[str, GeneSet],
    args: argparse.Namespace,
) -> pd.DataFrame:
    np.random.seed(args.random_seed)
    rows = []
    for keys, group in scores.groupby(["map_id", "tissue", "annotated_cell_type", "state_name"], sort=False):
        map_id, tissue, cell_type, state_name = keys
        gene_set = biological_sets[state_name]
        info = marker_info(gene_set, matrix.columns)
        score_iqr = group["ucell_score"].quantile(0.75) - group["ucell_score"].quantile(0.25)
        base = {
            "map_id": map_id,
            "tissue": tissue,
            "annotated_cell_type": cell_type,
            "state_name": state_name,
            "mixture_threshold": np.nan,
            "n_cells_in_calibration_group": len(group),
            "n_null_calibration_cells": np.nan,
            "null_calibration_max_cells": args.null_max_cells,
            "score_iqr": score_iqr,
            "n_markers_requested": info["n_markers_requested"],
            "n_markers_present": info["n_markers_present"],
            "marker_coverage_fraction": info["marker_coverage_fraction"],
        }
        if len(group) < args.min_calibration_cells:
            rows.append({**base, "threshold_method": "insufficient_cells", "threshold_value": np.nan, "null_95_threshold": np.nan, "null_99_threshold": np.nan, "threshold_status": "insufficient_cells"})
            continue
        if score_iqr < args.min_score_iqr:
            rows.append({**base, "threshold_method": "low_dynamic_range", "threshold_value": np.nan, "null_95_threshold": np.nan, "null_99_threshold": np.nan, "threshold_status": "low_dynamic_range"})
            continue
        present = info["present_list"]
        null_scores, n_null_calibration_cells = null_scores_for_state(
            matrix, ranks, group["cell_id"].tolist(), present, args, state_name
        )
        null_95 = float(np.quantile(null_scores, 0.95)) if null_scores else np.nan
        null_99 = float(np.quantile(null_scores, 0.99)) if null_scores else np.nan
        rows.append(
            {
                **base,
                "threshold_method": "matched_random_gene_set_null99",
                "threshold_value": null_99,
                "null_95_threshold": null_95,
                "null_99_threshold": null_99,
                "n_null_calibration_cells": n_null_calibration_cells,
                "threshold_status": "ok",
            }
        )
    return pd.DataFrame(rows)


def compute_probabilities(
    scores: pd.DataFrame,
    matrix: pd.DataFrame,
    ranks: pd.DataFrame,
    gene_sets: dict[str, GeneSet],
    args: argparse.Namespace,
    state_col: str,
    state_type: str,
) -> pd.DataFrame:
    rows = []
    score_state_col = "state_name" if state_col == "state_name" else "qc_signature_name"
    for state_name, group in scores.groupby(score_state_col, sort=False):
        gene_set = gene_sets[state_name]
        info = marker_info(gene_set, matrix.columns)
        cells = group["cell_id"].tolist()
        null_values, n_null_calibration_cells = null_scores_for_state(
            matrix, ranks, cells, info["present_list"], args, state_name
        )
        probs, method = probability_from_null(group["ucell_score"], null_values)
        out = group[["cell_id"]].copy()
        if {"map_id", "tissue", "annotated_cell_type"}.issubset(group.columns):
            out[["map_id", "tissue", "annotated_cell_type"]] = group[["map_id", "tissue", "annotated_cell_type"]]
        out["state_type"] = state_type
        out["state_name"] = state_name
        out["ucell_score"] = group["ucell_score"].to_numpy()
        out["state_probability"] = probs.to_numpy()
        out["probability_method"] = method
        out["n_null_sets"] = args.null_n
        out["n_null_calibration_cells"] = n_null_calibration_cells
        out["null_calibration_max_cells"] = args.null_max_cells
        out["marker_coverage_fraction"] = info["marker_coverage_fraction"]
        out["n_markers_present"] = info["n_markers_present"]
        out["n_markers_total"] = info["n_markers_requested"]
        rows.append(out)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def build_ucell_scores_output(bio_scores: pd.DataFrame, qc_scores: pd.DataFrame) -> pd.DataFrame:
    bio = bio_scores.rename(columns={"state_name": "state_name", "n_markers_requested": "n_markers_total"}).copy()
    bio["state_type"] = "biological"
    bio_out = bio[
        [
            "cell_id",
            "state_type",
            "state_name",
            "ucell_score",
            "n_markers_present",
            "n_markers_total",
            "marker_coverage_fraction",
            "markers_present",
            "markers_missing",
        ]
    ]
    qc = qc_scores.rename(columns={"qc_signature_name": "state_name", "n_markers_requested": "n_markers_total"}).copy()
    qc["state_type"] = "qc"
    if "markers_present" not in qc.columns:
        qc["markers_present"] = ""
    if "markers_missing" not in qc.columns:
        qc["markers_missing"] = ""
    qc_out = qc[
        [
            "cell_id",
            "state_type",
            "state_name",
            "ucell_score",
            "n_markers_present",
            "n_markers_total",
            "marker_coverage_fraction",
            "markers_present",
            "markers_missing",
        ]
    ]
    return pd.concat([bio_out, qc_out], ignore_index=True)


def confidence_from_margin(margin: float, exploratory: bool = False) -> str:
    if exploratory:
        return "exploratory"
    if pd.isna(margin):
        return "low"
    if margin >= 0.05:
        return "high"
    if margin >= 0:
        return "medium"
    return "low"


def hard_assignments_from_probabilities(
    probabilities: pd.DataFrame,
    state_thresholds: dict[str, float],
    qc_thresholds: dict[str, float],
    args: argparse.Namespace,
) -> pd.DataFrame:
    frame = probabilities.copy()
    state_threshold = frame["state_name"].map(state_thresholds)
    qc_threshold = frame["state_name"].map(qc_thresholds)
    frame["threshold"] = np.where(
        frame["state_type"].eq("qc"),
        qc_threshold.fillna(args.default_qc_threshold),
        state_threshold.fillna(args.default_state_threshold),
    )
    frame["threshold_source"] = np.where(
        frame["state_type"].eq("qc"),
        np.where(qc_threshold.notna(), "yaml", "default_qc"),
        np.where(state_threshold.notna(), "yaml", "default_biological"),
    )
    frame["marker_coverage_pass"] = (frame["n_markers_present"] >= args.min_markers_present) & (
        frame["marker_coverage_fraction"] >= args.min_marker_coverage
    )
    frame["hard_call"] = frame["marker_coverage_pass"] & (frame["state_probability"] >= frame["threshold"])
    frame["reason"] = np.select(
        [
            ~frame["marker_coverage_pass"],
            frame["hard_call"],
        ],
        [
            "insufficient_marker_coverage",
            "probability_above_threshold",
        ],
        default="probability_below_threshold",
    )
    return frame[
        [
            "cell_id",
            "state_type",
            "state_name",
            "state_probability",
            "threshold",
            "hard_call",
            "threshold_source",
            "marker_coverage_pass",
            "reason",
        ]
    ].copy()


def qc_exclusions(probabilities: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    qc = probabilities.loc[probabilities["state_type"] == "qc"].copy()
    cells = pd.DataFrame({"cell_id": sorted(probabilities["cell_id"].unique())})
    if args.exclude_qc_above is None:
        cells["excluded"] = False
        cells["exclusion_reason"] = "qc_exclusion_not_requested"
        cells["triggering_qc_states"] = "none"
        cells["max_qc_probability"] = cells["cell_id"].map(qc.groupby("cell_id")["state_probability"].max()).fillna(0.0)
        return cells
    selected = set(args.exclude_qc_states.split(",")) if args.exclude_qc_states else set(qc["state_name"].unique())
    qc = qc.loc[qc["state_name"].isin(selected)].copy()
    triggers = qc.loc[qc["state_probability"] >= args.exclude_qc_above]
    grouped = triggers.groupby("cell_id")["state_name"].agg(lambda x: ";".join(sorted(set(x))))
    cells["triggering_qc_states"] = cells["cell_id"].map(grouped).fillna("none")
    cells["excluded"] = cells["triggering_qc_states"].ne("none")
    cells["exclusion_reason"] = np.where(cells["excluded"], "qc_probability_above_threshold", "not_excluded")
    cells["max_qc_probability"] = cells["cell_id"].map(qc.groupby("cell_id")["state_probability"].max()).fillna(0.0)
    return cells


def state_probability_lookup(probabilities: pd.DataFrame) -> dict[tuple[str, str], float]:
    bio = probabilities.loc[probabilities["state_type"] == "biological"]
    return {(row.cell_id, row.state_name): float(row.state_probability) for row in bio.itertuples(index=False)}


def loo_probability_for_gene(
    gene: str,
    state_name: str,
    gene_set: GeneSet,
    matrix: pd.DataFrame,
    ranks: pd.DataFrame,
    args: argparse.Namespace,
) -> pd.Series | None:
    if gene not in gene_set.genes:
        return None
    loo_genes = [g for g in gene_set.genes if g != gene]
    present = [g for g in loo_genes if g in matrix.columns]
    if not present:
        return pd.Series(0.0, index=matrix.index)
    score = ucell_score_from_ranks(ranks, present, args.max_rank)
    null_values, _ = null_scores_for_state(matrix, ranks, list(matrix.index), present, args, f"{state_name}__loo__{gene}")
    probs, _ = probability_from_null(score, null_values)
    return probs


EXPECTED_EXPRESSION_COLUMNS = [
    "gene",
    "state_name",
    "state_type",
    "weighted_mean_expression",
    "weighted_detection_fraction",
    "whole_cell_type_mean_expression",
    "log2_weighted_vs_all",
    "n_cells",
    "sum_state_weight",
    "leave_one_gene_out_used",
]


def expected_expression_summary(
    matrix: pd.DataFrame,
    probabilities: pd.DataFrame,
    biological_sets: dict[str, GeneSet],
    ranks: pd.DataFrame,
    args: argparse.Namespace,
    query_genes: list[str],
) -> pd.DataFrame:
    rows = []
    if args.mode == "hard":
        return empty_table(EXPECTED_EXPRESSION_COLUMNS)
    all_mean = matrix[query_genes].mean(axis=0)
    bio_probs = probabilities.loc[probabilities["state_type"] == "biological"].copy()
    for state_name, group in bio_probs.groupby("state_name", sort=False):
        weights = group.set_index("cell_id")["state_probability"].reindex(matrix.index).fillna(0.0)
        for gene in query_genes:
            used_loo = False
            gene_weights = weights
            if args.leave_one_gene_out and state_name in biological_sets:
                loo = loo_probability_for_gene(gene, state_name, biological_sets[state_name], matrix, ranks, args)
                if loo is not None:
                    gene_weights = loo.reindex(matrix.index).fillna(0.0)
                    used_loo = True
            denom = gene_weights.sum()
            if denom > 0:
                weighted_mean = float((gene_weights * matrix[gene]).sum() / denom)
                weighted_detection = float((gene_weights * (matrix[gene] > 0)).sum() / denom)
            else:
                weighted_mean = np.nan
                weighted_detection = np.nan
            rows.append(
                {
                    "gene": gene,
                    "state_name": state_name,
                    "state_type": "biological",
                    "weighted_mean_expression": weighted_mean,
                    "weighted_detection_fraction": weighted_detection,
                    "whole_cell_type_mean_expression": float(all_mean[gene]),
                    "log2_weighted_vs_all": float(np.log2((weighted_mean + 1e-6) / (all_mean[gene] + 1e-6))) if pd.notna(weighted_mean) else np.nan,
                    "n_cells": int(len(matrix)),
                    "sum_state_weight": float(denom),
                    "leave_one_gene_out_used": used_loo,
                }
            )
    return pd.DataFrame(rows, columns=EXPECTED_EXPRESSION_COLUMNS)


HARD_EXPRESSION_COLUMNS = [
    "gene",
    "state_name",
    "state_type",
    "mean_expression_state_positive",
    "mean_expression_state_negative",
    "detection_fraction_positive",
    "detection_fraction_negative",
    "log2_fc_positive_vs_negative",
    "p_value",
    "q_value",
    "n_positive_cells",
    "n_negative_cells",
]


def hard_expression_summary(matrix: pd.DataFrame, hard: pd.DataFrame, args: argparse.Namespace, query_genes: list[str]) -> pd.DataFrame:
    rows = []
    if args.mode == "expected":
        return empty_table(HARD_EXPRESSION_COLUMNS)
    bio_hard = hard.loc[hard["state_type"] == "biological"].copy()
    for state_name, group in bio_hard.groupby("state_name", sort=False):
        calls = group.set_index("cell_id")["hard_call"].reindex(matrix.index).fillna(False).astype(bool)
        pos = matrix.loc[calls]
        neg = matrix.loc[~calls]
        for gene in query_genes:
            if len(pos) > 0 and len(neg) > 0:
                stat = stats.mannwhitneyu(pos[gene], neg[gene], alternative="two-sided").pvalue
            else:
                stat = np.nan
            rows.append(
                {
                    "gene": gene,
                    "state_name": state_name,
                    "state_type": "biological",
                    "mean_expression_state_positive": float(pos[gene].mean()) if len(pos) else np.nan,
                    "mean_expression_state_negative": float(neg[gene].mean()) if len(neg) else np.nan,
                    "detection_fraction_positive": float((pos[gene] > 0).mean()) if len(pos) else np.nan,
                    "detection_fraction_negative": float((neg[gene] > 0).mean()) if len(neg) else np.nan,
                    "log2_fc_positive_vs_negative": float(np.log2((pos[gene].mean() + 1e-6) / (neg[gene].mean() + 1e-6))) if len(pos) and len(neg) else np.nan,
                    "p_value": stat,
                    "n_positive_cells": int(len(pos)),
                    "n_negative_cells": int(len(neg)),
                }
            )
    out = pd.DataFrame(rows, columns=[c for c in HARD_EXPRESSION_COLUMNS if c != "q_value"])
    out["q_value"] = bh_fdr(out["p_value"]) if not out.empty else []
    return out[HARD_EXPRESSION_COLUMNS]


def donor_weighted_expression(matrix: pd.DataFrame, metadata: pd.DataFrame, weights: pd.Series, donor_col: str, query_genes: list[str]) -> pd.DataFrame:
    rows = []
    meta = metadata.set_index("cell_id").reindex(matrix.index)
    for donor, idx in meta.groupby(donor_col).groups.items():
        donor_cells = list(idx)
        w = weights.reindex(donor_cells).fillna(0.0)
        denom = w.sum()
        values = {"donor_id": donor, "sum_state_weight": float(denom)}
        if denom > 0:
            weighted = matrix.loc[donor_cells, query_genes].multiply(w, axis=0).sum(axis=0) / denom
        else:
            weighted = pd.Series(np.nan, index=query_genes)
        values.update(weighted.to_dict())
        rows.append(values)
    return pd.DataFrame(rows)


def fit_simple_lm(y: pd.Series, design: pd.DataFrame, phenotype: str) -> tuple[float, float, int]:
    frame = pd.concat([y.rename("y"), design], axis=1).dropna()
    if len(frame) < 3 or phenotype not in frame.columns or frame[phenotype].nunique() < 2:
        return np.nan, np.nan, len(frame)
    x = frame.drop(columns=["y"])
    x = pd.get_dummies(x, drop_first=True, dtype=float)
    x.insert(0, "intercept", 1.0)
    if phenotype not in x.columns:
        phen_cols = [c for c in x.columns if c.startswith(phenotype + "_")]
        if not phen_cols:
            return np.nan, np.nan, len(frame)
        coef_col = phen_cols[0]
    else:
        coef_col = phenotype
    X = x.to_numpy(float)
    Y = frame["y"].to_numpy(float)
    beta, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ beta
    dof = len(Y) - X.shape[1]
    if dof <= 0:
        return beta[list(x.columns).index(coef_col)], np.nan, len(frame)
    sigma2 = float((resid @ resid) / dof)
    cov = sigma2 * np.linalg.pinv(X.T @ X)
    idx = list(x.columns).index(coef_col)
    se = float(np.sqrt(cov[idx, idx])) if cov[idx, idx] >= 0 else np.nan
    p = 2 * stats.t.sf(abs(beta[idx] / se), dof) if se and np.isfinite(se) and se > 0 else np.nan
    return float(beta[idx]), float(p), len(frame)


def de_summaries(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    probabilities: pd.DataFrame,
    hard: pd.DataFrame,
    phenotypes: pd.DataFrame | None,
    args: argparse.Namespace,
    query_genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expected_cols = ["gene", "state_name", "phenotype", "coefficient", "coefficient_units", "p_value", "q_global", "q_by_trait", "q_by_gene", "n_donors", "sum_state_weight", "model_formula"]
    hard_cols = ["gene", "state_name", "phenotype", "coefficient", "coefficient_units", "p_value", "q_global", "q_by_trait", "q_by_gene", "n_donors", "n_state_positive_cells", "model_formula"]
    if args.mode == "expected":
        hard_requested = False
        expected_requested = True
    elif args.mode == "hard":
        hard_requested = True
        expected_requested = False
    else:
        hard_requested = True
        expected_requested = True
    if phenotypes is None or phenotypes.empty:
        return empty_table(expected_cols), empty_table(hard_cols)
    if args.donor_col not in phenotypes.columns:
        raise SystemExit(f"Phenotype table is missing donor column: {args.donor_col}")
    pheno_cols = [c for c in phenotypes.columns if c != args.donor_col]
    if not pheno_cols:
        return empty_table(expected_cols), empty_table(hard_cols)
    phenotype = pheno_cols[0]
    donor_design = phenotypes.set_index(args.donor_col)[pheno_cols]
    exp_rows = []
    hard_rows = []
    bio_probs = probabilities.loc[probabilities["state_type"] == "biological"]
    for state_name, group in bio_probs.groupby("state_name", sort=False):
        weights = group.set_index("cell_id")["state_probability"].reindex(matrix.index).fillna(0.0)
        if expected_requested:
            donor_expr = donor_weighted_expression(matrix, metadata, weights, args.donor_col, query_genes).set_index("donor_id")
            design = donor_design.reindex(donor_expr.index)
            for gene in query_genes:
                coef, p, n = fit_simple_lm(donor_expr[gene], design, phenotype)
                exp_rows.append(
                    {
                        "gene": gene,
                        "state_name": state_name,
                        "phenotype": phenotype,
                        "coefficient": coef,
                        "coefficient_units": "expression_per_phenotype_unit",
                        "p_value": p,
                        "n_donors": n,
                        "sum_state_weight": float(weights.sum()),
                        "model_formula": f"E_d[{gene}|{state_name}] ~ {phenotype}",
                    }
                )
        if hard_requested:
            calls = hard.loc[(hard["state_type"] == "biological") & (hard["state_name"] == state_name)].set_index("cell_id")["hard_call"].reindex(matrix.index).fillna(False).astype(bool)
            meta = metadata.set_index("cell_id").reindex(matrix.index)
            for gene in query_genes:
                donor_values = []
                for donor, idx in meta.groupby(args.donor_col).groups.items():
                    donor_cells = list(idx)
                    active_cells = [c for c in donor_cells if calls.loc[c]]
                    value = float(matrix.loc[active_cells, gene].mean()) if active_cells else np.nan
                    donor_values.append({"donor_id": donor, "value": value})
                donor_frame = pd.DataFrame(donor_values).set_index("donor_id")
                design = donor_design.reindex(donor_frame.index)
                coef, p, n = fit_simple_lm(donor_frame["value"], design, phenotype)
                hard_rows.append(
                    {
                        "gene": gene,
                        "state_name": state_name,
                        "phenotype": phenotype,
                        "coefficient": coef,
                        "coefficient_units": "expression_per_phenotype_unit",
                        "p_value": p,
                        "n_donors": n,
                        "n_state_positive_cells": int(calls.sum()),
                        "model_formula": f"mean_hard_positive[{gene}|{state_name}] ~ {phenotype}",
                    }
                )
    exp = pd.DataFrame(exp_rows, columns=[c for c in expected_cols if not c.startswith("q_")])
    hard_de = pd.DataFrame(hard_rows, columns=[c for c in hard_cols if not c.startswith("q_")])
    for frame in [exp, hard_de]:
        if frame.empty:
            continue
        frame["q_global"] = bh_fdr(frame["p_value"])
        frame["q_by_trait"] = frame.groupby("phenotype", group_keys=False)["p_value"].apply(bh_fdr)
        frame["q_by_gene"] = frame.groupby("gene", group_keys=False)["p_value"].apply(bh_fdr)
    return exp.reindex(columns=expected_cols), hard_de.reindex(columns=hard_cols)


def call_states(scores: pd.DataFrame, thresholds: pd.DataFrame, bad_flags: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    threshold_cols = [
        "map_id",
        "tissue",
        "annotated_cell_type",
        "state_name",
        "threshold_method",
        "threshold_value",
        "threshold_status",
        "score_iqr",
        "n_cells_in_calibration_group",
    ]
    frame = scores.merge(thresholds[threshold_cols], on=["map_id", "tissue", "annotated_cell_type", "state_name"], how="left")
    frame = frame.merge(bad_flags[["cell_id", "hard_exclusion_flag", "review_flag"]], on="cell_id", how="left")
    calls = []
    for _, row in frame.iterrows():
        requires_composite = bool(row["is_composite_state"])
        composite_rule = "none"
        if row["n_markers_present"] < args.min_markers_present or row["marker_coverage_fraction"] < args.min_marker_coverage:
            call = "insufficient_marker_coverage"
            reason = "insufficient_marker_coverage"
        elif row["threshold_status"] == "insufficient_cells":
            call = "not_called_insufficient_cells"
            reason = "calibration_group_below_min_cells"
        elif row["threshold_status"] == "low_dynamic_range":
            call = "not_called_low_dynamic_range"
            reason = "score_iqr_below_minimum"
        elif pd.isna(row["threshold_value"]):
            call = "not_called_missing_threshold"
            reason = "missing_threshold"
        elif row["ucell_score"] >= row["threshold_value"]:
            if bool(row["hard_exclusion_flag"]):
                call = "ambiguous_qc_flagged"
                reason = "score_high_but_hard_exclusion_qc_flag"
            elif requires_composite:
                call = "exploratory_marker_high"
                reason = "composite_state_requires_validation"
            else:
                call = "active"
                reason = "score_above_null99"
        else:
            call = "inactive"
            reason = "score_below_threshold"
        margin = row["ucell_score"] - row["threshold_value"] if not pd.isna(row["threshold_value"]) else np.nan
        calls.append(
            {
                "cell_id": row["cell_id"],
                "map_id": row["map_id"],
                "tissue": row["tissue"],
                "annotated_cell_type": row["annotated_cell_type"],
                "state_name": row["state_name"],
                "ucell_score": row["ucell_score"],
                "threshold_value": row["threshold_value"],
                "call": call,
                "confidence": confidence_from_margin(margin, call == "exploratory_marker_high"),
                "reason": reason,
                "hard_exclusion_flag": bool(row["hard_exclusion_flag"]),
                "review_flag": bool(row["review_flag"]),
                "composite_rule_used": composite_rule,
                "requires_composite_validation": requires_composite and composite_rule == "none",
                "state_kind": row["state_kind"],
            }
        )
    return pd.DataFrame(calls)


def multilabel_summary(calls: pd.DataFrame, bad_flags: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cell_id, group in calls.groupby("cell_id", sort=False):
        active_bio = sorted(group.loc[(group["call"] == "active") & (group["state_kind"] == "biological"), "state_name"].unique())
        process = sorted(group.loc[(group["call"] == "active") & (group["state_kind"] == "process"), "state_name"].unique())
        exploratory = sorted(group.loc[group["call"] == "exploratory_marker_high", "state_name"].unique())
        flags = bad_flags.loc[bad_flags["cell_id"] == cell_id].iloc[0]
        use_program = "exclude_hard_flagged" if flags["hard_exclusion_flag"] else "include"
        use_de = "exclude_hard_flagged" if flags["hard_exclusion_flag"] else "include"
        rows.append(
            {
                "cell_id": cell_id,
                "map_id": flags["map_id"],
                "tissue": flags["tissue"],
                "annotated_cell_type": flags["annotated_cell_type"],
                "active_biological_states": ";".join(active_bio) if active_bio else "none",
                "active_process_flags": ";".join(process) if process else "none",
                "exploratory_states": ";".join(exploratory) if exploratory else "none",
                "qc_flags": flags["bad_cell_reason"],
                "recommended_use_for_program_inference": use_program,
                "recommended_use_for_de": use_de,
                "summary_label": "; ".join(active_bio + process + exploratory) if (active_bio or process or exploratory) else "no_active_state_calls",
            }
        )
    return pd.DataFrame(rows)


def state_call_summary(calls: pd.DataFrame, metadata: pd.DataFrame, thresholds: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    active = calls["call"].isin(["active", "exploratory_marker_high"])
    frame = calls.copy()
    frame["is_active"] = active
    frame = frame.merge(metadata[["cell_id", args.donor_col]], on="cell_id", how="left")
    rows = []
    for keys, group in frame.groupby(["map_id", "tissue", "annotated_cell_type", "state_name"], sort=False):
        threshold = thresholds.set_index(["map_id", "tissue", "annotated_cell_type", "state_name"]).loc[keys]
        donor_active = group.groupby(args.donor_col)["is_active"].sum()
        rows.append(
            {
                "map_id": keys[0],
                "tissue": keys[1],
                "annotated_cell_type": keys[2],
                "state_name": keys[3],
                "n_cells_scored": len(group),
                "n_cells_active": int(group["is_active"].sum()),
                "active_fraction": float(group["is_active"].mean()) if len(group) else np.nan,
                "n_donors_with_active_cells": int((donor_active > 0).sum()),
                "median_active_cells_per_donor": float(donor_active.median()) if len(donor_active) else np.nan,
                "threshold_method": threshold["threshold_method"],
                "threshold_value": threshold["threshold_value"],
                "marker_coverage_fraction": threshold["marker_coverage_fraction"],
                "confidence_summary": ";".join(sorted(group["confidence"].astype(str).unique())),
            }
        )
    return pd.DataFrame(rows)


def acceptance_checks(summary: pd.DataFrame, scores: pd.DataFrame, qc_scores: pd.DataFrame, calls: pd.DataFrame) -> list[str]:
    failures = []
    for _, group in summary.groupby(["map_id", "tissue", "annotated_cell_type"], sort=False):
        if len(group) > 1 and group["active_fraction"].nunique(dropna=True) == 1:
            failures.append(
                "All biological states have the same active fraction for "
                f"{group.iloc[0]['map_id']}/{group.iloc[0]['tissue']}/{group.iloc[0]['annotated_cell_type']}"
            )
    if scores["marker_coverage_fraction"].isna().any():
        failures.append("Marker coverage is missing for at least one biological score row")
    if qc_scores.empty:
        failures.append("Auxiliary bad-cell QC signatures were not scored")
    if not {"hard_exclusion_flag", "review_flag"}.issubset(set(calls.columns)):
        failures.append("Hard-exclusion and review flags are not both present in calls")
    composite = calls.loc[calls["requires_composite_validation"]]
    if not composite.empty and not composite["call"].eq("exploratory_marker_high").all():
        failures.append("Composite states without rules were not labeled exploratory marker-high")
    return failures


def write_methods(out_dir: Path, args: argparse.Namespace, mapping_info: dict[str, str], failures: list[str]) -> None:
    lines = [
        "# CMDKP cell-state scoring method",
        "",
        f"- Input matrix type: `{args.expression_kind}`",
        "- UCell implementation: `local_ucell_style_rank_statistic`",
        "- UCell package version: `not_used`",
        f"- UCell parameters: `maxRank={args.max_rank}`, `ties.method=average`, `chunk.size=not_applicable_python_runner`, `missing_genes=skip`, `knn_smoothing=none`",
        f"- Biological GMT: `{args.biological_gmt}`",
        f"- Auxiliary bad-cell QC GMT: `{args.qc_gmt}`",
        f"- Gene ID handling: `{mapping_info['gene_id_type']}`",
        f"- Duplicate collapse method: `{mapping_info['duplicate_collapse_method']}`",
        f"- Marker coverage rules: score with at least 1 marker present; confident calls require `n_markers_present >= {args.min_markers_present}` and `marker_coverage_fraction >= {args.min_marker_coverage}`",
        f"- Thresholding/probability calibration method: matched random gene-set null, null95 and null99 reported, null99 used for legacy state calls, `null_n={args.null_n}`, `null_max_cells={args.null_max_cells}`",
        f"- Calibration group: `map_id + tissue + annotated_cell_type + state_name`; minimum cells `{args.min_calibration_cells}`; minimum score IQR `{args.min_score_iqr}`",
        f"- Summary query genes: `{args.query_genes or args.query_gene or 'all_expression_matrix_genes'}`",
        f"- Summary mode: `{args.mode}`",
        "- Bad-cell hard exclusions: technical metric outliers and QC signatures marked with hard-exclude tiers in the auxiliary GMT.",
        "- Bad-cell review flags: ribosomal/translation, heat-shock, immediate-early, and other review-tier QC signatures when extreme.",
        "- Scores are primary and calls are secondary. State calls are multi-label and not mutually exclusive.",
        "",
        "## Acceptance checks",
    ]
    if failures:
        lines.extend(f"- FAILED: {failure}" for failure in failures)
    else:
        lines.append("- All implemented acceptance checks passed.")
    (out_dir / "state_scoring_method.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", default="", help="Long or wide expression TSV/TSV.GZ")
    parser.add_argument("--expression-matrix", default="", help="Alias for --expression")
    parser.add_argument("--metadata", default="", help="Cell metadata TSV/TSV.GZ")
    parser.add_argument("--cell-metadata", default="", help="Alias for --metadata")
    parser.add_argument("--biological-gmt", default="")
    parser.add_argument("--states-gmt", default="", help="Alias for --biological-gmt")
    parser.add_argument("--qc-gmt", default="")
    parser.add_argument("--qc-states-gmt", default="", help="Alias for --qc-gmt")
    parser.add_argument("--out-dir", default="")
    parser.add_argument("--expression-kind", choices=["raw_counts", "log1p_normalized"], default="log1p_normalized")
    parser.add_argument("--phenotypes", default="")
    parser.add_argument("--state-thresholds-yaml", default="")
    parser.add_argument("--qc-thresholds-yaml", default="")
    parser.add_argument("--mode", choices=["expected", "hard", "both"], default="both")
    parser.add_argument("--exclude-qc-above", type=float, default=None)
    parser.add_argument("--exclude-qc-states", default="")
    parser.add_argument("--gene-map", default="")
    parser.add_argument("--duplicate-collapse", choices=["auto", "sum", "mean", "max"], default="auto")
    parser.add_argument("--map-id-col", default="map_id")
    parser.add_argument("--tissue-col", default="tissue")
    parser.add_argument("--cell-type-col", default="cell_type")
    parser.add_argument("--donor-col", default="donor_id")
    parser.add_argument("--sample-col", default="sample_id")
    parser.add_argument("--doublet-col", default="")
    parser.add_argument("--max-rank", type=int, default=1500)
    parser.add_argument("--null-n", type=int, default=1000)
    parser.add_argument("--null-max-cells", type=int, default=20000, help="Maximum cells used per state to estimate matched-null calibration backgrounds; <=0 uses all cells")
    parser.add_argument("--random-seed", type=int, default=1)
    parser.add_argument("--expression-bins", type=int, default=20)
    parser.add_argument("--detection-bins", type=int, default=5)
    parser.add_argument("--min-markers-present", type=int, default=5)
    parser.add_argument("--min-marker-coverage", type=float, default=0.5)
    parser.add_argument("--default-state-threshold", type=float, default=0.80)
    parser.add_argument("--default-qc-threshold", type=float, default=0.95)
    parser.add_argument("--leave-one-gene-out", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--query-genes", default="", help="Optional newline-delimited gene list restricting expression and DE summaries")
    parser.add_argument("--query-gene", default="", help="Optional comma-separated genes restricting expression and DE summaries")
    parser.add_argument("--min-calibration-cells", type=int, default=100)
    parser.add_argument("--min-score-iqr", type=float, default=0.01)
    parser.add_argument("--qc-extreme-percentile", type=float, default=0.99)
    parser.add_argument("--ucell-scores-out", default="")
    parser.add_argument("--cell-state-probabilities-out", default="")
    parser.add_argument("--cell-state-hard-assignments-out", default="")
    parser.add_argument("--qc-exclusions-out", default="")
    parser.add_argument("--expression-expected-assignments-out", default="")
    parser.add_argument("--expression-hard-assignments-out", default="")
    parser.add_argument("--de-expected-assignments-out", default="")
    parser.add_argument("--de-hard-assignments-out", default="")
    parser.add_argument("--state-summary-out", default="")
    parser.add_argument("--run-summary-out", default="")
    parser.add_argument("--allow-acceptance-failures", action="store_true")
    args = parser.parse_args()

    args.expression = args.expression or args.expression_matrix
    args.metadata = args.metadata or args.cell_metadata
    args.biological_gmt = args.biological_gmt or args.states_gmt
    args.qc_gmt = args.qc_gmt or args.qc_states_gmt or "out/qc/cmdkp_all_tissues_minimal_bad_cell_qc_signatures.gmt"
    if not args.expression or not args.metadata or not args.biological_gmt:
        raise SystemExit("--expression-matrix, --cell-metadata, and --states-gmt are required")
    if not args.out_dir:
        if not all(
            [
                args.ucell_scores_out,
                args.cell_state_probabilities_out,
                args.cell_state_hard_assignments_out,
                args.qc_exclusions_out,
                args.expression_expected_assignments_out,
                args.expression_hard_assignments_out,
                args.de_expected_assignments_out,
                args.de_hard_assignments_out,
                args.state_summary_out,
                args.run_summary_out,
            ]
        ):
            raise SystemExit("Provide --out-dir or all explicit output paths")
        args.out_dir = str(Path(args.ucell_scores_out).parent)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    metadata = read_table(args.metadata)
    if args.cell_type_col not in metadata.columns and "annotated_cell_type" in metadata.columns:
        args.cell_type_col = "annotated_cell_type"
    required_meta = ["cell_id", args.map_id_col, args.tissue_col, args.cell_type_col, args.donor_col]
    if args.sample_col not in metadata.columns:
        metadata[args.sample_col] = metadata[args.donor_col]
    missing = [c for c in required_meta if c not in metadata.columns]
    if missing:
        raise SystemExit(f"Metadata is missing required column(s): {', '.join(missing)}")
    metadata = metadata.drop_duplicates("cell_id").copy()
    expression, mapping_info = harmonize_expression(read_expression_input(args.expression), args.expression_kind, args.gene_map, args.duplicate_collapse)
    matrix = expression_matrix(expression, metadata["cell_id"])
    query_genes = resolve_query_genes(args, matrix)
    ranks = rank_matrix_for_ucell(matrix, args.max_rank)

    biological_sets = read_gmt(args.biological_gmt)
    qc_sets = read_gmt(args.qc_gmt)
    biological_scores = score_biological_states(matrix, ranks, metadata, biological_sets, args)
    qc_scores = score_qc_signatures(matrix, ranks, metadata, qc_sets, args)
    thresholds = calibrate_thresholds(biological_scores, matrix, ranks, {s.name: s for s in biological_sets}, args)
    bad_flags = qc_metrics(matrix, metadata, qc_scores, args)
    bio_prob = compute_probabilities(biological_scores, matrix, ranks, {s.name: s for s in biological_sets}, args, "state_name", "biological")
    qc_prob = compute_probabilities(qc_scores, matrix, ranks, {s.name: s for s in qc_sets}, args, "qc_signature_name", "qc")
    probabilities = pd.concat([bio_prob, qc_prob], ignore_index=True)
    state_thresholds = load_threshold_yaml(args.state_thresholds_yaml)
    qc_thresholds = load_threshold_yaml(args.qc_thresholds_yaml)
    hard = hard_assignments_from_probabilities(probabilities, state_thresholds, qc_thresholds, args)
    exclusions = qc_exclusions(probabilities, args)
    expression_expected = expected_expression_summary(matrix, probabilities, {s.name: s for s in biological_sets}, ranks, args, query_genes)
    expression_hard = hard_expression_summary(matrix, hard, args, query_genes)
    phenotypes = read_table(args.phenotypes) if args.phenotypes else None
    de_expected, de_hard = de_summaries(matrix, metadata, probabilities, hard, phenotypes, args, query_genes)

    calls = call_states(biological_scores, thresholds, bad_flags, args)
    multilabel = multilabel_summary(calls, bad_flags)
    summary = state_call_summary(calls, metadata, thresholds, args)
    state_summary = probabilities.groupby(["state_type", "state_name"], sort=False).agg(
        n_cells_scored=("cell_id", "nunique"),
        mean_ucell_score=("ucell_score", "mean"),
        mean_probability=("state_probability", "mean"),
        marker_coverage_fraction=("marker_coverage_fraction", "first"),
    ).reset_index()
    hard_counts = hard.groupby(["state_type", "state_name"], sort=False).agg(
        n_hard_assigned=("hard_call", "sum"),
        threshold=("threshold", "first"),
    ).reset_index()
    state_summary = state_summary.merge(hard_counts, on=["state_type", "state_name"], how="left")
    state_summary["hard_assigned_fraction"] = state_summary["n_hard_assigned"] / state_summary["n_cells_scored"]
    failures = acceptance_checks(summary, biological_scores, qc_scores, calls)

    ucell_scores = build_ucell_scores_output(biological_scores, qc_scores)
    write_table(ucell_scores, output_path(args, "ucell_scores_out", "ucell_scores.tsv.gz"))
    write_table(probabilities[["cell_id", "state_type", "state_name", "ucell_score", "state_probability", "probability_method", "n_null_sets", "n_null_calibration_cells", "null_calibration_max_cells", "marker_coverage_fraction"]], output_path(args, "cell_state_probabilities_out", "cell_state_probabilities.tsv.gz"))
    write_table(hard, output_path(args, "cell_state_hard_assignments_out", "cell_state_hard_assignments.tsv.gz"))
    write_table(exclusions, output_path(args, "qc_exclusions_out", "qc_exclusions.tsv.gz"))
    write_table(expression_expected, output_path(args, "expression_expected_assignments_out", "expression_expected_assignments.tsv.gz"))
    write_table(expression_hard, output_path(args, "expression_hard_assignments_out", "expression_hard_assignments.tsv.gz"))
    write_table(de_expected, output_path(args, "de_expected_assignments_out", "de_expected_assignments.tsv.gz"))
    write_table(de_hard, output_path(args, "de_hard_assignments_out", "de_hard_assignments.tsv.gz"))
    write_table(state_summary[["state_type", "state_name", "n_cells_scored", "mean_ucell_score", "mean_probability", "n_hard_assigned", "hard_assigned_fraction", "marker_coverage_fraction", "threshold"]], output_path(args, "state_summary_out", "state_summary.tsv.gz"))
    write_table(biological_scores.drop(columns=["scope_tissue", "scope_cell_type", "state_kind", "is_composite_state"]), out_dir / "cell_state_scores.tsv.gz")
    write_table(thresholds, out_dir / "cell_state_thresholds.tsv.gz")
    write_table(calls.drop(columns=["state_kind"]), out_dir / "cell_state_calls.tsv.gz")
    write_table(qc_scores.drop(columns=["qc_tier", "qc_category"]), out_dir / "qc_signature_scores.tsv.gz")
    write_table(bad_flags, out_dir / "bad_cell_qc_flags.tsv.gz")
    write_table(multilabel, out_dir / "cell_multilabel_state_summary.tsv.gz")
    write_table(summary, out_dir / "state_call_summary.tsv.gz")
    write_methods(out_dir, args, mapping_info, failures)

    run_summary = {
        "input_files": {
            "expression_matrix": args.expression,
            "cell_metadata": args.metadata,
            "states_gmt": args.biological_gmt,
            "qc_states_gmt": args.qc_gmt,
            "phenotypes": args.phenotypes or None,
        },
        "parameters": {
            "mode": args.mode,
            "max_rank": args.max_rank,
            "null_n": args.null_n,
            "null_max_cells": args.null_max_cells,
            "query_genes": args.query_genes or None,
            "query_gene": args.query_gene or None,
            "n_query_genes": len(query_genes),
            "default_state_threshold": args.default_state_threshold,
            "default_qc_threshold": args.default_qc_threshold,
            "leave_one_gene_out": args.leave_one_gene_out,
            "exclude_qc_above": args.exclude_qc_above,
        },
        "n_cells": int(matrix.shape[0]),
        "n_genes": int(matrix.shape[1]),
        "n_states": int(len(biological_sets)),
        "n_qc_states": int(len(qc_sets)),
        "n_excluded_cells": int(exclusions["excluded"].sum()),
        "software_versions": {"python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__},
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    output_path(args, "run_summary_out", "run_summary.json").write_text(json.dumps(run_summary, indent=2), encoding="utf-8")

    if failures and not args.allow_acceptance_failures:
        raise SystemExit("Acceptance checks failed:\n" + "\n".join(failures))


if __name__ == "__main__":
    main()
