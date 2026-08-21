#!/usr/bin/env python3
"""Match data-driven programs to curated cell states and QC signatures."""

from __future__ import annotations

import warnings
import argparse
import gzip
import hashlib
import json
import math
import platform
import re
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import scipy
from scipy import stats

warnings.filterwarnings("ignore", message="Mean of empty slice", category=RuntimeWarning)


PROGRAM_ALIASES = ("program_id", "program", "factor", "module", "component")
GENE_ALIASES = ("gene", "gene_symbol", "gene_name")
LOADING_ALIASES = ("loading", "weight", "value", "loading_score")
CELL_ALIASES = ("cell_id", "cell", "barcode")
ACTIVITY_ALIASES = ("program_activity", "usage", "score", "weight", "activity")
STATE_ACTIVITY_COLUMNS = ["state_activity_weight_gradient", "state_activity_weight_hightail", "aucell_score"]
MATCH_SUMMARY_COLUMNS = [
    "tissue",
    "cell_type",
    "program_id",
    "state_id",
    "state_label",
    "state_type",
    "n_program_genes",
    "n_state_markers",
    "n_state_markers_in_program_universe",
    "marker_coverage_fraction",
    "gsea_nes",
    "gsea_p",
    "gsea_q",
    "loading_auc",
    "loading_mwu_q",
    "leading_edge_genes",
    "top100_overlap_n",
    "top100_overlap_genes",
    "cell_spearman_r_gradient",
    "cell_spearman_q_gradient",
    "donor_spearman_r_gradient",
    "donor_spearman_q_gradient",
    "expression_score_spearman_r",
    "expression_score_spearman_q",
    "best_gene_level_score",
    "best_cell_level_score",
    "combined_match_score",
    "match_class",
    "interpretation",
    "qc_caveat",
]


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else path.open("r", encoding="utf-8")


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def write_table(frame: pd.DataFrame, path: Path, columns: list[str] | None = None) -> None:
    if columns is not None:
        for col in columns:
            if col not in frame.columns:
                frame[col] = np.nan
        frame = frame[columns]
    frame.to_csv(path, sep="\t", index=False, compression="infer")


def norm_col(value: str) -> str:
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", str(value).strip().lower())).strip("_")


def find_col(frame: pd.DataFrame, aliases: tuple[str, ...]) -> str | None:
    by_norm = {norm_col(col): col for col in frame.columns}
    for alias in aliases:
        col = by_norm.get(norm_col(alias))
        if col is not None:
            return col
    return None


def display_label(value: str) -> str:
    text = re.sub(r"\s+", " ", str(value).replace("_", " ")).strip()
    return text[:1].upper() + text[1:] if text else ""


def parse_bool(value: str | None, default: bool) -> bool:
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def parse_list(value: str, cast=str) -> list[Any]:
    return [cast(x.strip()) for x in str(value).split(",") if x.strip()]


def canonical_cell_id(value: Any) -> str:
    text = str(value).strip()
    if not text:
        return ""
    text = re.sub(r"_1$", "", text)
    if text.startswith("SRR"):
        text = text.replace("_", "-")
    return text


def read_gmt(path: Path, state_type: str) -> pd.DataFrame:
    rows = []
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            state_id = parts[0]
            markers = list(dict.fromkeys(g for g in parts[2:] if g))
            rows.append(
                {
                    "state_id": state_id,
                    "state_label": display_label(state_id),
                    "state_type": state_type,
                    "description": parts[1],
                    "markers": markers,
                    "n_state_markers": len(markers),
                }
            )
    return pd.DataFrame(rows)


def pair_seed(base_seed: int, program_id: str, state_id: str) -> int:
    digest = hashlib.sha256(f"{base_seed}|{program_id}|{state_id}".encode()).digest()
    return int.from_bytes(digest[:8], "little")


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


def add_q_values(frame: pd.DataFrame, p_col: str, q_col: str, group_cols: list[str] | None = None) -> pd.DataFrame:
    frame[q_col] = bh_fdr(frame[p_col])
    frame[f"{q_col}_global"] = frame[q_col]
    if group_cols:
        for _, idx in frame.groupby(group_cols, dropna=False).groups.items():
            frame.loc[idx, f"{q_col}_within_state_type"] = bh_fdr(frame.loc[idx, p_col])
    return frame


def load_program_loadings(path: Path) -> pd.DataFrame:
    frame = read_table(path)
    program_col = find_col(frame, PROGRAM_ALIASES)
    gene_col = find_col(frame, GENE_ALIASES)
    loading_col = find_col(frame, LOADING_ALIASES)
    if program_col and gene_col and loading_col:
        out = frame[[program_col, gene_col, loading_col]].rename(columns={program_col: "program_id", gene_col: "gene", loading_col: "loading"})
    else:
        gene_col = frame.columns[0]
        out = frame.melt(id_vars=[gene_col], var_name="program_id", value_name="loading").rename(columns={gene_col: "gene"})
    out["program_id"] = out["program_id"].astype(str).str.strip()
    out["gene"] = out["gene"].astype(str).str.strip()
    out["loading"] = pd.to_numeric(out["loading"], errors="coerce")
    out = out.loc[out["program_id"].ne("") & out["gene"].ne("") & out["loading"].notna()].copy()
    out = out.groupby(["program_id", "gene"], as_index=False)["loading"].max()
    return out


def load_program_cell_activity(path: Path) -> pd.DataFrame:
    frame = read_table(path)
    program_col = find_col(frame, PROGRAM_ALIASES)
    cell_col = find_col(frame, CELL_ALIASES)
    activity_col = find_col(frame, ACTIVITY_ALIASES)
    if program_col and cell_col and activity_col:
        out = frame[[cell_col, program_col, activity_col]].rename(columns={cell_col: "cell_id", program_col: "program_id", activity_col: "program_activity"})
    else:
        cell_col = frame.columns[0]
        out = frame.melt(id_vars=[cell_col], var_name="program_id", value_name="program_activity").rename(columns={cell_col: "cell_id"})
    out["cell_id"] = out["cell_id"].astype(str).str.strip()
    out["cell_id_match"] = out["cell_id"].map(canonical_cell_id)
    out["program_id"] = out["program_id"].astype(str).str.strip()
    out["program_activity"] = pd.to_numeric(out["program_activity"], errors="coerce")
    return out.loc[out["cell_id"].ne("") & out["cell_id_match"].ne("") & out["program_id"].ne("") & out["program_activity"].notna()].copy()


def gsea_es_for_hit_indices(loadings: np.ndarray, hit_indices: np.ndarray, weight: float) -> tuple[float, int]:
    n = len(loadings)
    nh = len(hit_indices)
    hit_mask = np.zeros(n, dtype=bool)
    hit_mask[hit_indices] = True
    hit_weights = np.abs(loadings[hit_indices]) ** weight
    hit_total = hit_weights.sum()
    if hit_total <= 0:
        hit_weights = np.ones(nh, dtype=float) / nh
    else:
        hit_weights = hit_weights / hit_total
    increments = np.full(n, -1.0 / (n - nh), dtype=float)
    increments[hit_indices] = hit_weights
    running = np.cumsum(increments)
    max_idx = int(np.argmax(running))
    return float(max(0.0, running[max_idx])), max_idx


def marker_enrichment_for_pair(
    program: str,
    loadings: pd.DataFrame,
    state: pd.Series,
    args: argparse.Namespace,
    base_seed: int,
    top_n_list: list[int],
    tissue: str,
    cell_type: str,
) -> dict[str, Any]:
    loadings = loadings.sort_values(["loading", "gene"], ascending=[False, True]).reset_index(drop=True)
    genes = loadings["gene"].tolist()
    loading_values = loadings["loading"].to_numpy(float)
    universe = set(genes)
    markers = list(state["markers"])
    marker_set = set(markers)
    hits = [g for g in markers if g in universe]
    missing = [g for g in markers if g not in universe]
    hit_positions = np.array([i for i, gene in enumerate(genes) if gene in set(hits)], dtype=int)
    n_program = len(genes)
    n_markers = len(markers)
    nh = len(hit_positions)
    coverage = nh / n_markers if n_markers else np.nan
    status = "ok"
    if nh < args.min_marker_overlap or (not np.isnan(coverage) and coverage < args.min_marker_coverage) or nh == n_program:
        status = "insufficient_marker_coverage"

    row: dict[str, Any] = {
        "tissue": tissue,
        "cell_type": cell_type,
        "program_id": program,
        "state_id": state["state_id"],
        "state_label": state["state_label"],
        "state_type": state["state_type"],
        "n_program_genes": n_program,
        "n_state_markers": n_markers,
        "n_state_markers_in_program_universe": nh,
        "marker_coverage_fraction": coverage,
        "missing_state_markers": ";".join(missing),
        "gsea_es": np.nan,
        "gsea_nes": np.nan,
        "gsea_p": np.nan,
        "loading_auc": np.nan,
        "loading_mwu_p": np.nan,
        "leading_edge_genes": "",
        "match_status": status,
    }

    for top_n in top_n_list:
        effective_n = min(top_n, n_program)
        top_genes = genes[:effective_n]
        overlap = [g for g in top_genes if g in marker_set]
        row[f"top{top_n}_overlap_n"] = len(overlap)
        row[f"top{top_n}_overlap_genes"] = ";".join(overlap)
        row[f"top{top_n}_jaccard"] = len(overlap) / len(set(top_genes).union(marker_set)) if top_genes or marker_set else np.nan
        row[f"top{top_n}_hypergeom_p"] = stats.hypergeom.sf(len(overlap) - 1, n_program, nh, effective_n) if nh and effective_n else np.nan

    if status == "insufficient_marker_coverage":
        return row

    rng = np.random.default_rng(pair_seed(base_seed, program, state["state_id"]))
    es, leading_idx = gsea_es_for_hit_indices(loading_values, hit_positions, args.gsea_weight)
    null_es = np.empty(args.gsea_permutations, dtype=float)
    for i in range(args.gsea_permutations):
        random_hits = rng.choice(n_program, size=nh, replace=False)
        null_es[i], _ = gsea_es_for_hit_indices(loading_values, np.sort(random_hits), args.gsea_weight)
    mean_null = float(np.nanmean(null_es))
    p_value = (1 + int(np.sum(null_es >= es))) / (1 + args.gsea_permutations)
    leading_edge = [gene for i, gene in enumerate(genes[: leading_idx + 1]) if i in set(hit_positions)]
    marker_loadings = loadings.loc[loadings["gene"].isin(hits), "loading"].to_numpy(float)
    nonmarker_loadings = loadings.loc[~loadings["gene"].isin(hits), "loading"].to_numpy(float)
    if len(marker_loadings) and len(nonmarker_loadings):
        mwu = stats.mannwhitneyu(marker_loadings, nonmarker_loadings, alternative="greater")
        auc = float(mwu.statistic / (len(marker_loadings) * len(nonmarker_loadings)))
        row["loading_auc"] = auc
        row["loading_mwu_p"] = float(mwu.pvalue)
    row["gsea_es"] = es
    row["gsea_nes"] = es / mean_null if mean_null > 0 else np.nan
    row["gsea_p"] = p_value
    row["leading_edge_genes"] = ";".join(leading_edge)
    return row


def compute_marker_enrichment(program_loadings: pd.DataFrame, states: pd.DataFrame, args: argparse.Namespace, tissue: str, cell_type: str) -> pd.DataFrame:
    top_n_list = sorted(set(parse_list(args.program_top_n_list, int) + [50, 100, 200]))
    rows = []
    for program, loadings in program_loadings.groupby("program_id", sort=False):
        for _, state in states.iterrows():
            rows.append(marker_enrichment_for_pair(program, loadings, state, args, args.random_seed, top_n_list, tissue, cell_type))
    out = pd.DataFrame(rows)
    for p_col, q_col in [("gsea_p", "gsea_q"), ("loading_mwu_p", "loading_mwu_q")]:
        out[q_col] = bh_fdr(out[p_col])
        for _, idx in out.groupby("state_type", dropna=False).groups.items():
            out.loc[idx, f"{q_col}_within_state_type"] = bh_fdr(out.loc[idx, p_col])
    return out


def score_from_state_expression(frame: pd.DataFrame, score_col: str) -> pd.Series:
    if score_col == "-log10(q_value)":
        q = pd.to_numeric(frame["q_value"], errors="coerce").clip(lower=1e-300)
        sign = np.sign(pd.to_numeric(frame.get("log2fc_weighted_vs_all_parent", 1), errors="coerce").fillna(0))
        return -np.log10(q) * sign
    if score_col not in frame.columns:
        raise SystemExit(f"--state-gene-score-col is not present in --state-expression: {score_col}")
    return pd.to_numeric(frame[score_col], errors="coerce")


def compute_expression_match(program_loadings: pd.DataFrame, state_expression: pd.DataFrame | None, states: pd.DataFrame, args: argparse.Namespace, tissue: str, cell_type: str) -> pd.DataFrame:
    columns = [
        "tissue",
        "cell_type",
        "program_id",
        "state_id",
        "state_label",
        "state_weight_type",
        "state_gene_score_col",
        "n_program_genes",
        "n_state_expression_genes",
        "n_intersection_genes",
        "spearman_r",
        "spearman_p",
        "spearman_q",
        "top_program_gene_auc",
        "top_program_gene_mwu_p",
        "top_program_gene_mwu_q",
        "top_program_gene_mean_state_score",
        "top_program_genes",
        "top_program_genes_with_state_scores",
        "match_status",
    ]
    if state_expression is None:
        return pd.DataFrame(columns=columns)
    state_expression = state_expression.copy()
    state_expression["state_id"] = state_expression["state_name"].astype(str)
    rows = []
    for weight_type in parse_list(args.state_weight_types):
        expr_weight = state_expression.loc[state_expression["state_weight_type"].astype(str).eq(weight_type)].copy()
        if expr_weight.empty:
            continue
        expr_weight["state_gene_score"] = score_from_state_expression(expr_weight, args.state_gene_score_col)
        for (state_id, state_label), expr_state in expr_weight.groupby(["state_id", "state_name"], sort=False):
            expr_state = expr_state.dropna(subset=["state_gene_score"])
            if expr_state.empty:
                continue
            score_by_gene = expr_state.groupby("gene")["state_gene_score"].mean()
            for program, loadings in program_loadings.groupby("program_id", sort=False):
                loadings = loadings.sort_values(["loading", "gene"], ascending=[False, True])
                top_genes = loadings["gene"].head(args.program_top_n).tolist()
                merged = loadings.set_index("gene")[["loading"]].join(score_by_gene.rename("state_gene_score"), how="inner")
                top_with_scores = [g for g in top_genes if g in score_by_gene.index]
                row = {
                    "tissue": tissue,
                    "cell_type": cell_type,
                    "program_id": program,
                    "state_id": state_id,
                    "state_label": display_label(state_label),
                    "state_weight_type": weight_type,
                    "state_gene_score_col": args.state_gene_score_col,
                    "n_program_genes": loadings["gene"].nunique(),
                    "n_state_expression_genes": int(score_by_gene.shape[0]),
                    "n_intersection_genes": int(merged.shape[0]),
                    "spearman_r": np.nan,
                    "spearman_p": np.nan,
                    "top_program_gene_auc": np.nan,
                    "top_program_gene_mwu_p": np.nan,
                    "top_program_gene_mean_state_score": np.nan,
                    "top_program_genes": ";".join(top_genes),
                    "top_program_genes_with_state_scores": ";".join(top_with_scores),
                    "match_status": "ok" if merged.shape[0] >= 3 else "insufficient_intersection_genes",
                }
                if merged.shape[0] >= 3:
                    corr = stats.spearmanr(merged["loading"], merged["state_gene_score"], nan_policy="omit")
                    row["spearman_r"] = float(corr.statistic)
                    row["spearman_p"] = float(corr.pvalue)
                top_scores = score_by_gene.loc[top_with_scores].dropna()
                rest_scores = score_by_gene.drop(index=top_with_scores, errors="ignore").dropna()
                if len(top_scores):
                    row["top_program_gene_mean_state_score"] = float(top_scores.mean())
                if len(top_scores) and len(rest_scores):
                    mwu = stats.mannwhitneyu(top_scores.to_numpy(float), rest_scores.to_numpy(float), alternative="greater")
                    row["top_program_gene_auc"] = float(mwu.statistic / (len(top_scores) * len(rest_scores)))
                    row["top_program_gene_mwu_p"] = float(mwu.pvalue)
                rows.append(row)
    out = pd.DataFrame(rows, columns=columns)
    if not out.empty:
        out["spearman_q"] = bh_fdr(out["spearman_p"])
        out["top_program_gene_mwu_q"] = bh_fdr(out["top_program_gene_mwu_p"])
    return out


def compute_cell_correlations(
    program_activity: pd.DataFrame | None,
    cell_state_activity: pd.DataFrame | None,
    metadata: pd.DataFrame | None,
    states: pd.DataFrame,
    args: argparse.Namespace,
    tissue: str,
    cell_type: str,
) -> pd.DataFrame:
    columns = [
        "tissue",
        "cell_type",
        "program_id",
        "state_id",
        "state_label",
        "state_type",
        "state_activity_column",
        "n_cells",
        "n_donors",
        "cell_spearman_r",
        "cell_spearman_p",
        "cell_spearman_q",
        "donor_spearman_r",
        "donor_spearman_p",
        "donor_spearman_q",
        "program_activity_mean",
        "state_activity_mean",
        "match_status",
    ]
    if program_activity is None or cell_state_activity is None:
        return pd.DataFrame(columns=columns)
    state_type = states.set_index("state_id")["state_type"].to_dict()
    state_label = states.set_index("state_id")["state_label"].to_dict()
    state_activity = cell_state_activity.copy()
    state_activity["state_id"] = state_activity["state_name"].astype(str)
    state_activity["cell_id"] = state_activity["cell_id"].astype(str)
    state_activity["cell_id_match"] = state_activity["cell_id"].map(canonical_cell_id)
    if "state_type" not in state_activity.columns:
        state_activity["state_type"] = state_activity["state_id"].map(state_type).fillna("curated_state")
    rows = []
    meta = None
    if metadata is not None and args.donor_col in metadata.columns and "cell_id" in metadata.columns:
        meta = metadata[["cell_id", args.donor_col]].dropna().copy()
        meta["cell_id"] = meta["cell_id"].astype(str)
        meta["cell_id_match"] = meta["cell_id"].map(canonical_cell_id)
    for activity_col in STATE_ACTIVITY_COLUMNS:
        if activity_col not in state_activity.columns:
            continue
        state_sub = state_activity[["cell_id", "cell_id_match", "state_id", "state_type", activity_col]].rename(columns={activity_col: "state_activity"})
        state_sub["state_activity"] = pd.to_numeric(state_sub["state_activity"], errors="coerce")
        for program, prog in program_activity.groupby("program_id", sort=False):
            prog = prog[["cell_id", "cell_id_match", "program_activity"]]
            for state_id, st in state_sub.groupby("state_id", sort=False):
                merged = prog.merge(st[["cell_id_match", "state_activity"]], on="cell_id_match", how="inner").dropna()
                row = {
                    "tissue": tissue,
                    "cell_type": cell_type,
                    "program_id": program,
                    "state_id": state_id,
                    "state_label": state_label.get(state_id, display_label(state_id)),
                    "state_type": state_type.get(state_id, st["state_type"].iloc[0] if "state_type" in st else "curated_state"),
                    "state_activity_column": activity_col,
                    "n_cells": int(len(merged)),
                    "n_donors": 0,
                    "cell_spearman_r": np.nan,
                    "cell_spearman_p": np.nan,
                    "donor_spearman_r": np.nan,
                    "donor_spearman_p": np.nan,
                    "program_activity_mean": float(merged["program_activity"].mean()) if len(merged) else np.nan,
                    "state_activity_mean": float(merged["state_activity"].mean()) if len(merged) else np.nan,
                    "match_status": "ok" if len(merged) >= 3 else "insufficient_cells",
                }
                if len(merged) >= 3 and merged["program_activity"].nunique() > 1 and merged["state_activity"].nunique() > 1:
                    corr = stats.spearmanr(merged["program_activity"], merged["state_activity"], nan_policy="omit")
                    row["cell_spearman_r"] = float(corr.statistic)
                    row["cell_spearman_p"] = float(corr.pvalue)
                if meta is not None and len(merged):
                    donor = merged.merge(meta[["cell_id_match", args.donor_col]], on="cell_id_match", how="inner")
                    donor = donor.groupby(args.donor_col)[["program_activity", "state_activity"]].mean()
                    row["n_donors"] = int(len(donor))
                    if len(donor) >= 3 and donor["program_activity"].nunique() > 1 and donor["state_activity"].nunique() > 1:
                        corr = stats.spearmanr(donor["program_activity"], donor["state_activity"], nan_policy="omit")
                        row["donor_spearman_r"] = float(corr.statistic)
                        row["donor_spearman_p"] = float(corr.pvalue)
                rows.append(row)
    out = pd.DataFrame(rows, columns=columns)
    if not out.empty:
        for col in ["cell_spearman_p", "donor_spearman_p"]:
            q_col = col.replace("_p", "_q")
            out[q_col] = bh_fdr(out[col])
            for _, idx in out.groupby(["state_type", "state_activity_column"], dropna=False).groups.items():
                out.loc[idx, f"{q_col}_within_state_type"] = bh_fdr(out.loc[idx, col])
    return out


def neglog10_q(q: Any) -> float:
    if pd.isna(q):
        return 0.0
    return min(50.0, -math.log10(max(float(q), 1e-300)))


def qc_caveat_for_state(state_id: str) -> str:
    text = norm_col(state_id)
    if any(x in text for x in ["ribosomal", "translation"]):
        return "ribosomal_or_translation"
    if any(x in text for x in ["mitochondrial", "apoptosis", "cell_death", "dying"]):
        return "mitochondrial_or_dying_cell"
    if any(x in text for x in ["offtarget", "off_target", "lineage"]):
        return "off_target_identity"
    if any(x in text for x in ["ambient", "contamination", "doublet"]):
        return "ambient_or_contamination"
    if any(x in text for x in ["heat_shock", "dissociation"]):
        return "heat_shock_or_dissociation"
    if any(x in text for x in ["immediate_early", "fos", "jun"]):
        return "immediate_early"
    return "none"


def interpretation_for(match_class: str) -> str:
    return {
        "strong_state_match": "Program loading genes are enriched for state markers and program activity tracks state activity.",
        "gene_only_state_match": "Program matches state markers but cell-level coactivity is weak or unavailable.",
        "cell_only_coactivity": "Program activity tracks state activity but marker overlap is weak.",
        "qc_dominated": "Program is QC-dominated; do not label as biological without review.",
        "mixed_state_qc": "Program has both biological state and QC/artifact evidence.",
        "insufficient_marker_coverage": "Program has insufficient marker coverage for this state.",
        "unmatched": "Program has no strong curated-state or QC match.",
    }.get(match_class, "Program has no strong curated-state or QC match.")


def build_summary(marker: pd.DataFrame, expr: pd.DataFrame, corr: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary = marker.copy()
    gradient = corr.loc[corr.get("state_activity_column", pd.Series(dtype=str)).eq("state_activity_weight_gradient")].copy() if not corr.empty else pd.DataFrame()
    if not gradient.empty:
        gradient = gradient.rename(
            columns={
                "cell_spearman_r": "cell_spearman_r_gradient",
                "cell_spearman_q": "cell_spearman_q_gradient",
                "donor_spearman_r": "donor_spearman_r_gradient",
                "donor_spearman_q": "donor_spearman_q_gradient",
            }
        )
        summary = summary.merge(
            gradient[["program_id", "state_id", "cell_spearman_r_gradient", "cell_spearman_q_gradient", "donor_spearman_r_gradient", "donor_spearman_q_gradient"]],
            on=["program_id", "state_id"],
            how="left",
        )
    expr_primary = expr.loc[expr.get("state_weight_type", pd.Series(dtype=str)).eq("gradient_percentile_squared")].copy() if not expr.empty else pd.DataFrame()
    if not expr_primary.empty:
        expr_primary = expr_primary.sort_values(["program_id", "state_id", "spearman_q", "spearman_r"], ascending=[True, True, True, False]).drop_duplicates(["program_id", "state_id"])
        expr_primary = expr_primary.rename(columns={"spearman_r": "expression_score_spearman_r", "spearman_q": "expression_score_spearman_q"})
        summary = summary.merge(expr_primary[["program_id", "state_id", "expression_score_spearman_r", "expression_score_spearman_q"]], on=["program_id", "state_id"], how="left")
    for col in ["cell_spearman_r_gradient", "cell_spearman_q_gradient", "donor_spearman_r_gradient", "donor_spearman_q_gradient", "expression_score_spearman_r", "expression_score_spearman_q"]:
        if col not in summary.columns:
            summary[col] = np.nan
    summary["best_gene_level_score"] = pd.to_numeric(summary["gsea_nes"], errors="coerce").clip(lower=0).fillna(0) * summary["gsea_q"].map(neglog10_q)
    summary["best_cell_level_score"] = pd.to_numeric(summary["cell_spearman_r_gradient"], errors="coerce").clip(lower=0).fillna(0) * summary["cell_spearman_q_gradient"].map(neglog10_q)
    summary["combined_match_score"] = summary["best_gene_level_score"] + summary["best_cell_level_score"]
    summary["qc_caveat"] = np.where(summary["state_type"].eq("qc_state"), summary["state_id"].map(qc_caveat_for_state), "none")

    best_curated = summary.loc[summary["state_type"].eq("curated_state")].sort_values("combined_match_score", ascending=False).drop_duplicates("program_id")
    best_qc = summary.loc[summary["state_type"].eq("qc_state")].sort_values("combined_match_score", ascending=False).drop_duplicates("program_id")
    best_curated_score = best_curated.set_index("program_id")["combined_match_score"].to_dict()
    best_qc_score = best_qc.set_index("program_id")["combined_match_score"].to_dict()
    best_qc_sig = best_qc.set_index("program_id")["gsea_q"].to_dict()
    best_qc_cell_q = best_qc.set_index("program_id")["cell_spearman_q_gradient"].to_dict()

    classes = []
    for _, row in summary.iterrows():
        program = row["program_id"]
        gsea_q = row.get("gsea_q")
        cell_q = row.get("cell_spearman_q_gradient")
        cell_r = row.get("cell_spearman_r_gradient")
        marker_ok = row.get("match_status") != "insufficient_marker_coverage"
        is_qc = row["state_type"] == "qc_state"
        qc_dominated = best_qc_score.get(program, -np.inf) > best_curated_score.get(program, -np.inf) and (
            (not pd.isna(best_qc_sig.get(program, np.nan)) and best_qc_sig[program] <= 0.05)
            or (not pd.isna(best_qc_cell_q.get(program, np.nan)) and best_qc_cell_q[program] <= 0.05)
        )
        strong_curated = (
            row["state_type"] == "curated_state"
            and marker_ok
            and not pd.isna(gsea_q)
            and gsea_q <= 0.05
            and row.get("gsea_nes", 0) > 0
            and not pd.isna(cell_r)
            and cell_r >= 0.20
        )
        strong_qc_for_program = qc_dominated and best_curated_score.get(program, 0) > 0 and best_qc_score.get(program, 0) > 0
        if row.get("match_status") == "insufficient_marker_coverage":
            cls = "insufficient_marker_coverage"
        elif is_qc and qc_dominated:
            cls = "qc_dominated"
        elif strong_curated and strong_qc_for_program:
            cls = "mixed_state_qc"
        elif strong_curated:
            cls = "strong_state_match"
        elif row["state_type"] == "curated_state" and marker_ok and not pd.isna(gsea_q) and gsea_q <= 0.05:
            cls = "gene_only_state_match"
        elif row["state_type"] == "curated_state" and not pd.isna(cell_q) and cell_q <= 0.05 and not pd.isna(cell_r) and cell_r >= 0.30 and (pd.isna(gsea_q) or gsea_q > 0.05 or not marker_ok):
            cls = "cell_only_coactivity"
        else:
            cls = "unmatched"
        classes.append(cls)
    summary["match_class"] = classes
    summary["interpretation"] = summary["match_class"].map(interpretation_for)

    qc_summary = build_qc_summary(summary)
    labels = build_label_suggestions(summary)
    return summary[MATCH_SUMMARY_COLUMNS], qc_summary, labels


def build_qc_summary(summary: pd.DataFrame) -> pd.DataFrame:
    qc = summary.loc[summary["state_type"].eq("qc_state")].copy()
    cols = [
        "program_id",
        "best_qc_state_id",
        "best_qc_label",
        "best_qc_gsea_q",
        "best_qc_gsea_nes",
        "best_qc_cell_spearman_r",
        "best_qc_cell_spearman_q",
        "qc_combined_match_score",
        "qc_caveat",
        "qc_recommendation",
    ]
    if qc.empty:
        return pd.DataFrame(columns=cols)
    rows = []
    for program, group in qc.groupby("program_id", sort=False):
        best = group.sort_values("combined_match_score", ascending=False).iloc[0]
        caveats = sorted(set(x for x in group.loc[group["combined_match_score"] > 0, "qc_caveat"].dropna() if x != "none"))
        caveat = "mixed_qc" if len(caveats) > 1 else caveats[0] if caveats else best["qc_caveat"]
        significant = (not pd.isna(best["gsea_q"]) and best["gsea_q"] <= 0.05) or (not pd.isna(best["cell_spearman_q_gradient"]) and best["cell_spearman_q_gradient"] <= 0.05)
        recommendation = "suppress_or_hide_by_default" if best["match_class"] == "qc_dominated" and significant else "review" if significant else "pass"
        rows.append(
            {
                "program_id": program,
                "best_qc_state_id": best["state_id"],
                "best_qc_label": best["state_label"],
                "best_qc_gsea_q": best["gsea_q"],
                "best_qc_gsea_nes": best["gsea_nes"],
                "best_qc_cell_spearman_r": best["cell_spearman_r_gradient"],
                "best_qc_cell_spearman_q": best["cell_spearman_q_gradient"],
                "qc_combined_match_score": best["combined_match_score"],
                "qc_caveat": caveat,
                "qc_recommendation": recommendation,
            }
        )
    return pd.DataFrame(rows, columns=cols)


def build_label_suggestions(summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    curated = summary.loc[summary["state_type"].eq("curated_state")]
    qc = summary.loc[summary["state_type"].eq("qc_state")]
    programs = sorted(summary["program_id"].unique())
    for program in programs:
        c = curated.loc[curated["program_id"].eq(program)].sort_values("combined_match_score", ascending=False)
        q = qc.loc[qc["program_id"].eq(program)].sort_values("combined_match_score", ascending=False)
        best_c = c.iloc[0] if not c.empty else None
        best_q = q.iloc[0] if not q.empty else None
        qc_caveat = best_q["qc_caveat"] if best_q is not None else "none"
        if best_q is not None and best_q["match_class"] == "qc_dominated":
            label = f"{qc_caveat.replace('_', ' ')}/QC program" if qc_caveat != "none" else f"{best_q['state_label']} QC program"
            quality = "qc_or_artifact"
        elif best_c is not None and best_c["match_class"] == "strong_state_match":
            label = f"{best_c['state_label']}-like program"
            quality = "high_confidence_biological"
        elif best_c is not None and best_c["match_class"] == "gene_only_state_match":
            label = f"{best_c['state_label']}-enriched program"
            quality = "exploratory_biological"
        elif best_c is not None and best_c["match_class"] == "cell_only_coactivity":
            label = f"{best_c['state_label']}-coactive program"
            quality = "exploratory_biological"
        elif best_c is not None and best_c["match_class"] == "mixed_state_qc":
            label = f"{best_c['state_label']}-like mixed QC program"
            quality = "mixed_state_qc"
        else:
            label = "unmatched data-driven program"
            quality = "unmatched"
        rows.append(
            {
                "program_id": program,
                "best_curated_state_id": best_c["state_id"] if best_c is not None else "",
                "best_curated_state_label": best_c["state_label"] if best_c is not None else "",
                "best_curated_match_class": best_c["match_class"] if best_c is not None else "",
                "best_qc_state_id": best_q["state_id"] if best_q is not None else "",
                "best_qc_label": best_q["state_label"] if best_q is not None else "",
                "qc_caveat": qc_caveat,
                "suggested_program_label": label,
                "suggested_program_quality_class": quality,
            }
        )
    return pd.DataFrame(rows)


def infer_label(frame: pd.DataFrame | None, column: str) -> str:
    if frame is not None and column in frame.columns:
        values = frame[column].dropna().astype(str).unique()
        if len(values) == 1:
            return values[0]
    return ""


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--program-loadings", type=Path, required=True)
    ap.add_argument("--state-gmt", type=Path, required=True)
    ap.add_argument("--qc-gmt", type=Path, default=None)
    ap.add_argument("--program-cell-activity", type=Path, default=None)
    ap.add_argument("--cell-state-activity", type=Path, default=None)
    ap.add_argument("--state-expression", type=Path, default=None)
    ap.add_argument("--metadata", type=Path, default=None)
    ap.add_argument("--donor-col", default="donor_id")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--tissue", default="")
    ap.add_argument("--cell-type", default="")
    ap.add_argument("--program-top-n", type=int, default=100)
    ap.add_argument("--program-top-n-list", default="50,100,200")
    ap.add_argument("--gsea-permutations", type=int, default=1000)
    ap.add_argument("--gsea-weight", type=float, default=1.0)
    ap.add_argument("--min-marker-overlap", type=int, default=3)
    ap.add_argument("--min-marker-coverage", type=float, default=0.2)
    ap.add_argument("--correlation-method", default="spearman", choices=["spearman"])
    ap.add_argument("--state-weight-types", default="gradient_percentile_squared,high_tail_percentile_90_100")
    ap.add_argument("--state-gene-score-col", default="log2fc_weighted_vs_all_parent")
    ap.add_argument("--random-seed", type=int, default=1)
    ap.add_argument("--include-qc", default=None, choices=["true", "false"])
    ap.add_argument("--api-minimal-output", action="store_true", help="Write only compact program-state outputs needed by portal APIs.")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    warnings: list[str] = []
    program_loadings = load_program_loadings(args.program_loadings)
    states = read_gmt(args.state_gmt, "curated_state")
    include_qc = parse_bool(args.include_qc, args.qc_gmt is not None)
    if include_qc and args.qc_gmt:
        states = pd.concat([states, read_gmt(args.qc_gmt, "qc_state")], ignore_index=True)
    qc_count = int(states["state_type"].eq("qc_state").sum())

    program_activity = load_program_cell_activity(args.program_cell_activity) if args.program_cell_activity else None
    cell_state_activity = read_table(args.cell_state_activity) if args.cell_state_activity else None
    state_expression = read_table(args.state_expression) if args.state_expression else None
    metadata = read_table(args.metadata) if args.metadata else None
    if program_activity is None:
        warnings.append("missing_program_cell_activity")
    if cell_state_activity is None:
        warnings.append("missing_cell_state_activity")
    if state_expression is None:
        warnings.append("missing_state_expression")

    tissue = args.tissue or infer_label(cell_state_activity, "tissue") or infer_label(state_expression, "tissue")
    cell_type = args.cell_type or infer_label(cell_state_activity, "annotated_cell_type") or infer_label(cell_state_activity, "cell_type") or infer_label(state_expression, "annotated_cell_type")

    marker = compute_marker_enrichment(program_loadings, states, args, tissue, cell_type)
    expr = compute_expression_match(program_loadings, state_expression, states, args, tissue, cell_type)
    corr = compute_cell_correlations(program_activity, cell_state_activity, metadata, states, args, tissue, cell_type)
    summary, qc_summary, labels = build_summary(marker, expr, corr)

    heat_long = summary[["program_id", "state_id", "state_type", "combined_match_score", "match_class", "qc_caveat"]].copy()
    heat_wide = heat_long.pivot_table(index="program_id", columns="state_id", values="combined_match_score", fill_value=0, aggfunc="max").reset_index()

    if marker["match_status"].eq("insufficient_marker_coverage").any():
        warnings.append("low_marker_coverage_states")
    small_programs = marker.loc[marker["n_program_genes"] < args.min_marker_overlap, "program_id"].drop_duplicates().tolist()
    if small_programs:
        warnings.append("programs_with_small_gene_universe")

    if args.api_minimal_output:
        heat_api = summary.copy()
        heat_api["correlation"] = heat_api.get("cell_spearman_r_gradient", np.nan)
        heat_cols = ["tissue", "cell_type", "state_id", "program_id", "correlation", "gsea_p", "gsea_q"]
        for col in heat_cols:
            if col not in heat_api.columns:
                heat_api[col] = "" if col in {"tissue", "cell_type", "state_id", "program_id"} else np.nan
        write_table(heat_api[heat_cols], args.out_dir / "program_state_heatmap_long.tsv.gz")
        factor_cols = ["program_id", "best_curated_state_id", "best_curated_state_label", "best_curated_match_class", "best_qc_state_id", "best_qc_label", "qc_caveat", "suggested_program_label", "suggested_program_quality_class"]
        for col in factor_cols:
            if col not in labels.columns:
                labels[col] = ""
        write_table(labels[factor_cols], args.out_dir / "program_label_suggestions.tsv.gz")
    else:
        write_table(marker, args.out_dir / "program_state_marker_enrichment.tsv.gz")
        write_table(expr, args.out_dir / "program_state_expression_score_match.tsv.gz")
        write_table(corr, args.out_dir / "program_state_cell_correlation.tsv.gz")
        write_table(summary, args.out_dir / "program_state_match_summary.tsv.gz")
        write_table(heat_wide, args.out_dir / "program_state_heatmap_matrix.tsv.gz")
        write_table(heat_long, args.out_dir / "program_state_heatmap_long.tsv.gz")
        if qc_count:
            write_table(qc_summary, args.out_dir / "program_qc_match_summary.tsv.gz")
        write_table(labels, args.out_dir / "program_label_suggestions.tsv.gz")

    run_summary = {
        "input_files": {
            "program_loadings": str(args.program_loadings),
            "state_gmt": str(args.state_gmt),
            "qc_gmt": str(args.qc_gmt) if args.qc_gmt else "",
            "program_cell_activity": str(args.program_cell_activity) if args.program_cell_activity else "",
            "cell_state_activity": str(args.cell_state_activity) if args.cell_state_activity else "",
            "state_expression": str(args.state_expression) if args.state_expression else "",
            "metadata": str(args.metadata) if args.metadata else "",
        },
        "n_programs": int(program_loadings["program_id"].nunique()),
        "n_curated_states": int(states["state_type"].eq("curated_state").sum()),
        "n_qc_states": qc_count,
        "n_program_state_tests": int(marker.shape[0]),
        "n_cell_correlation_tests": int(corr.shape[0]),
        "parameters": vars(args) | {"include_qc_resolved": include_qc},
        "software_versions": {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__, "scipy": scipy.__version__},
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "warnings": sorted(set(warnings)),
    }
    (args.out_dir / "run_summary.json").write_text(json.dumps(run_summary, indent=2, sort_keys=True, default=str) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
