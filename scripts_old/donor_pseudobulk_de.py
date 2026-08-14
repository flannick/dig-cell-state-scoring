#!/usr/bin/env python3
"""Run donor-level pseudobulk differential expression for a cell state."""

from __future__ import annotations

import argparse
import math

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def read_genes(path: str) -> list[str]:
    genes = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            gene = line.strip()
            if gene and not gene.startswith("#"):
                genes.append(gene)
    return list(dict.fromkeys(genes))


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series
    truthy = {"true", "t", "1", "yes", "y"}
    return series.astype(str).str.strip().str.lower().isin(truthy)


def complete_expression_grid(state_meta: pd.DataFrame, expr: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    cell_ids = state_meta["cell_id"].drop_duplicates()
    grid = pd.MultiIndex.from_product([cell_ids, genes], names=["cell_id", "gene"]).to_frame(index=False)
    expr = expr.loc[expr["gene"].isin(genes), ["cell_id", "gene", "expression"]].copy()
    expr["expression"] = pd.to_numeric(expr["expression"], errors="coerce")
    complete = grid.merge(expr, on=["cell_id", "gene"], how="left")
    complete["expression"] = complete["expression"].fillna(0.0)
    return complete


def bh_adjust(pvalues: list[float]) -> list[float]:
    pvals = np.asarray([np.nan if p is None else p for p in pvalues], dtype=float)
    qvals = np.full(len(pvals), np.nan)
    valid = np.where(~np.isnan(pvals))[0]
    if len(valid) == 0:
        return qvals.tolist()
    order = valid[np.argsort(pvals[valid])]
    ranked = pvals[order] * len(valid) / np.arange(1, len(valid) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    qvals[order] = np.minimum(ranked, 1.0)
    return qvals.tolist()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, help="Cell metadata TSV/TSV.GZ")
    parser.add_argument("--expression", required=True, help="Long expression TSV/TSV.GZ")
    parser.add_argument("--states", required=True, help="State membership TSV/TSV.GZ")
    parser.add_argument("--genes", required=True, help="Genes to test, one per line")
    parser.add_argument("--state", required=True, help="State name to test")
    parser.add_argument("--out", required=True, help="Output DE TSV/TSV.GZ")
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--donor-col", default="donor_id")
    parser.add_argument("--group-col", default="disease_group")
    parser.add_argument("--case", default="T2D")
    parser.add_argument("--control", default="control")
    parser.add_argument("--min-cells-per-donor", type=int, default=20)
    parser.add_argument("--min-donors-per-group", type=int, default=5)
    parser.add_argument("--pseudocount", type=float, default=1e-6)
    args = parser.parse_args()

    genes = read_genes(args.genes)
    if not genes:
        raise SystemExit("Gene list is empty")

    metadata = read_table(args.metadata)
    for column in [args.cell_id_col, args.donor_col, args.group_col]:
        if column not in metadata.columns:
            raise SystemExit(f"Metadata is missing required column: {column}")

    states = read_table(args.states)
    for column in ["cell_id", "state", "in_state"]:
        if column not in states.columns:
            raise SystemExit(f"State table is missing required column: {column}")

    membership = states.loc[(states["state"] == args.state) & as_bool(states["in_state"]), ["cell_id"]]
    if membership.empty:
        raise SystemExit(f"No cells assigned to state: {args.state}")

    state_meta = metadata.merge(membership, left_on=args.cell_id_col, right_on="cell_id", how="inner")
    state_meta = state_meta.loc[state_meta[args.group_col].isin([args.case, args.control])].copy()
    if state_meta.empty:
        raise SystemExit("No state cells remained after filtering to case/control groups")

    donor_counts = (
        state_meta.groupby([args.donor_col, args.group_col], as_index=False)
        .size()
        .rename(columns={"size": "n_cells"})
    )
    keep_donors = donor_counts.loc[
        donor_counts["n_cells"] >= args.min_cells_per_donor, [args.donor_col, args.group_col]
    ]
    state_meta = state_meta.merge(keep_donors, on=[args.donor_col, args.group_col], how="inner")
    if state_meta.empty:
        raise SystemExit("No donors passed the minimum cell count filter")

    expr = read_table(args.expression)
    missing_expr = {"cell_id", "gene", "expression"} - set(expr.columns)
    if missing_expr:
        raise SystemExit(f"Expression table is missing required columns: {sorted(missing_expr)}")
    state_meta = state_meta[[args.cell_id_col, args.donor_col, args.group_col]].rename(columns={args.cell_id_col: "cell_id"})
    expr = complete_expression_grid(state_meta, expr, genes)

    joined = expr.merge(
        state_meta,
        on="cell_id",
        how="inner",
    )
    if joined.empty:
        raise SystemExit("No expression rows overlapped the selected state cells")

    donor_means = (
        joined.groupby(["gene", args.donor_col, args.group_col], as_index=False)["expression"]
        .mean()
        .rename(columns={"expression": "donor_mean"})
    )

    rows = []
    for gene in genes:
        gene_means = donor_means.loc[donor_means["gene"] == gene]
        case_values = gene_means.loc[gene_means[args.group_col] == args.case, "donor_mean"].dropna().to_numpy()
        control_values = gene_means.loc[gene_means[args.group_col] == args.control, "donor_mean"].dropna().to_numpy()

        case_mean = float(np.mean(case_values)) if len(case_values) else math.nan
        control_mean = float(np.mean(control_values)) if len(control_values) else math.nan
        if len(case_values) >= args.min_donors_per_group and len(control_values) >= args.min_donors_per_group:
            pvalue = float(mannwhitneyu(case_values, control_values, alternative="two-sided").pvalue)
        else:
            pvalue = math.nan
        log2_fc = math.log2((case_mean + args.pseudocount) / (control_mean + args.pseudocount)) if (
            not math.isnan(case_mean) and not math.isnan(control_mean)
        ) else math.nan

        rows.append(
            {
                "gene": gene,
                "state": args.state,
                "case": args.case,
                "control": args.control,
                "case_donors": len(case_values),
                "control_donors": len(control_values),
                "case_mean": case_mean,
                "control_mean": control_mean,
                "log2_fc": log2_fc,
                "pvalue": pvalue,
                "test": "donor_mean_mannwhitneyu",
            }
        )

    out = pd.DataFrame(rows)
    out["qvalue"] = bh_adjust(out["pvalue"].tolist())
    out = out[
        [
            "gene",
            "state",
            "case",
            "control",
            "case_donors",
            "control_donors",
            "case_mean",
            "control_mean",
            "log2_fc",
            "pvalue",
            "qvalue",
            "test",
        ]
    ]
    out.to_csv(args.out, sep="\t", index=False, compression="infer")


if __name__ == "__main__":
    main()
