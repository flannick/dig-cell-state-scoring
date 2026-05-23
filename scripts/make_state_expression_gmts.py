#!/usr/bin/env python3
"""Build state-derived GMT files from all-gene state expression summaries."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def write_table(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, compression="infer")


def safe_name(value: str) -> str:
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))
    return re.sub(r"_+", "_", value).strip("_")


def read_gmt(path: Path) -> list[list[str]]:
    rows = []
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                rows.append(parts)
    return rows


def write_gmt(rows: list[list[str]], path: Path) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        for row in rows:
            handle.write("\t".join(row) + "\n")


def filter_frame(frame: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    out = frame.copy()
    if args.weight_type:
        out = out.loc[out["state_weight_type"].eq(args.weight_type)]
    if args.exclude_low_signal_states and "expression_result_scope" in out.columns:
        out = out.loc[~out["expression_result_scope"].eq("low_signal_flagged")]
    if args.exclude_composite_required and "state_class" in out.columns:
        out = out.loc[~out["state_class"].eq("composite_required")]
    if args.min_weighted_mean_cp10k is not None:
        out = out.loc[pd.to_numeric(out["weighted_mean_cp10k"], errors="coerce").fillna(0) >= args.min_weighted_mean_cp10k]
    if args.min_weighted_pct_detected is not None:
        out = out.loc[pd.to_numeric(out["weighted_pct_detected"], errors="coerce").fillna(0) >= args.min_weighted_pct_detected]
    return out


def make_rows(frame: pd.DataFrame, method: str, args: argparse.Namespace) -> tuple[list[list[str]], list[dict[str, object]]]:
    metric_map = {
        "top_absolute_expression": ("weighted_mean_cp10k", False),
        "top_specific_fc": ("log2fc_weighted_vs_all_parent", False),
        "top_specific_logp": ("specific_logp", False),
        "top_specific_combined": ("combined_score", False),
    }
    rows: list[list[str]] = []
    membership: list[dict[str, object]] = []
    work = frame.copy()
    if method in {"top_specific_fc", "top_specific_logp", "top_specific_combined"}:
        if args.min_log2fc is not None:
            work = work.loc[pd.to_numeric(work["log2fc_weighted_vs_all_parent"], errors="coerce").fillna(-np.inf) >= args.min_log2fc]
        if args.max_q_value is not None and "q_value" in work.columns:
            work = work.loc[pd.to_numeric(work["q_value"], errors="coerce").fillna(np.inf) <= args.max_q_value]
        if args.max_p_value is not None and "p_value" in work.columns:
            work = work.loc[pd.to_numeric(work["p_value"], errors="coerce").fillna(np.inf) <= args.max_p_value]
    p_for_log = pd.to_numeric(work.get("q_value", work.get("p_value")), errors="coerce").clip(lower=1e-300)
    work["specific_logp"] = -np.log10(p_for_log)
    work["combined_score"] = pd.to_numeric(work.get("log2fc_weighted_vs_all_parent"), errors="coerce").fillna(0) * work["specific_logp"].fillna(0)
    metric, ascending = metric_map[method]
    for state, group in work.groupby("state_name", sort=False):
        group = group.sort_values(metric, ascending=ascending).head(args.top_n)
        name = f"{safe_name(state)}__{method}__top{args.top_n}"
        genes = group["gene"].astype(str).tolist()
        rows.append([name, f"{method} top {args.top_n} for {state}", *genes])
        for rank, (_, row) in enumerate(group.iterrows(), 1):
            membership.append(
                {
                    "gene_set": name,
                    "state_name": state,
                    "signature_method": method,
                    "rank": rank,
                    "gene": row["gene"],
                    "selection_metric": metric,
                    "selection_value": row.get(metric, np.nan),
                    "state_weight_type": row.get("state_weight_type", ""),
                }
            )
    return rows, membership


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--state-expression-specificity", type=Path, required=True)
    ap.add_argument("--original-state-gmt", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--top-n", type=int, default=250)
    ap.add_argument("--min-weighted-mean-cp10k", type=float, default=None)
    ap.add_argument("--min-weighted-pct-detected", type=float, default=None)
    ap.add_argument("--min-log2fc", type=float, default=None)
    ap.add_argument("--max-q-value", type=float, default=None)
    ap.add_argument("--max-p-value", type=float, default=None)
    ap.add_argument("--exclude-low-signal-states", action="store_true")
    ap.add_argument("--exclude-composite-required", action="store_true")
    ap.add_argument("--weight-type", default="gradient_percentile_squared")
    ap.add_argument("--include-combined", action="store_true")
    args = ap.parse_args()

    gmt_dir = args.out_dir / "gmt"
    gmt_dir.mkdir(parents=True, exist_ok=True)
    frame = filter_frame(read_table(args.state_expression_specificity), args)
    original = [[safe_name(row[0]), row[1], *row[2:]] for row in read_gmt(args.original_state_gmt)]
    write_gmt(original, gmt_dir / "original_markers.gmt")
    membership = [
        {"gene_set": row[0], "state_name": row[0], "signature_method": "original_markers", "rank": i, "gene": gene, "selection_metric": "curated_marker", "selection_value": np.nan, "state_weight_type": ""}
        for row in original
        for i, gene in enumerate(row[2:], 1)
    ]
    summary = [{"signature_method": "original_markers", "n_gene_sets": len(original), "path": str(gmt_dir / "original_markers.gmt")}]
    methods = ["top_absolute_expression", "top_specific_fc", "top_specific_logp"]
    if args.include_combined:
        methods.append("top_specific_combined")
    for method in methods:
        rows, member_rows = make_rows(frame, method, args)
        path = gmt_dir / f"{method}.gmt"
        write_gmt(rows, path)
        membership.extend(member_rows)
        summary.append({"signature_method": method, "n_gene_sets": len(rows), "path": str(path)})
    write_table(pd.DataFrame(membership), args.out_dir / "gmt_membership.tsv.gz")
    write_table(pd.DataFrame(summary), args.out_dir / "gmt_build_summary.tsv")


if __name__ == "__main__":
    main()
