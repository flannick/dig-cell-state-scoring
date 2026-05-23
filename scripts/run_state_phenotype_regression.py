#!/usr/bin/env python3
"""Run donor-level phenotype regression for parent and state-weighted expression."""

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
    out = []
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            out.append(parts[1] if len(parts) > 1 and parts[1] else parts[0])
    return out


def resolve(directory: Path, names: list[str]) -> Path:
    for name in names:
        path = directory / name
        if path.exists():
            return path
    raise SystemExit(f"Could not find {names} in {directory}")


def load_10x(directory: Path) -> tuple[sparse.csr_matrix, list[str], list[str]]:
    genes = read_features(resolve(directory, ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"]))
    cells = read_one_column(resolve(directory, ["barcodes.tsv.gz", "barcodes.tsv"]))
    mat = mmread(resolve(directory, ["matrix.mtx.gz", "matrix.mtx"])).tocsr()
    if mat.shape == (len(genes), len(cells)):
        mat = mat.T.tocsr()
    elif mat.shape != (len(cells), len(genes)):
        raise SystemExit("Matrix shape does not match genes/cells")
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


def fit_lm(y: pd.Series, design: pd.DataFrame, phenotype: str) -> tuple[float, float, int, int, int, str]:
    frame = pd.concat([y.rename("y"), design], axis=1).dropna()
    if len(frame) < 3 or phenotype not in frame.columns or frame[phenotype].nunique() < 2:
        return np.nan, np.nan, len(frame), 0, 0, ""
    x = frame.drop(columns=["y"]).copy()
    if pd.api.types.is_numeric_dtype(x[phenotype]):
        sd = x[phenotype].std(ddof=0)
        if sd > 0:
            x[phenotype] = (x[phenotype] - x[phenotype].mean()) / sd
        coef_col = phenotype
        n_cases = n_controls = 0
    else:
        levels = sorted(x[phenotype].astype(str).unique())
        ref = levels[0]
        x[phenotype] = pd.Categorical(x[phenotype].astype(str), categories=levels)
        n_controls = int((x[phenotype].astype(str) == ref).sum())
        n_cases = int(len(x) - n_controls)
        coef_col = f"{phenotype}_{levels[1]}" if len(levels) > 1 else ""
    x = pd.get_dummies(x, drop_first=True, dtype=float)
    x.insert(0, "intercept", 1.0)
    if coef_col not in x.columns:
        matches = [c for c in x.columns if c.startswith(phenotype + "_")]
        coef_col = matches[0] if matches else phenotype
    if coef_col not in x.columns:
        return np.nan, np.nan, len(frame), n_cases, n_controls, ""
    X = x.to_numpy(float)
    Y = frame["y"].to_numpy(float)
    coef, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ coef
    dof = len(Y) - X.shape[1]
    idx = list(x.columns).index(coef_col)
    if dof <= 0:
        return float(coef[idx]), np.nan, len(frame), n_cases, n_controls, coef_col
    sigma2 = float((resid @ resid) / dof)
    cov = sigma2 * np.linalg.pinv(X.T @ X)
    se = float(np.sqrt(cov[idx, idx])) if cov[idx, idx] >= 0 else np.nan
    p = 2 * stats.t.sf(abs(coef[idx] / se), dof) if se and np.isfinite(se) and se > 0 else np.nan
    return float(coef[idx]), float(p), len(frame), n_cases, n_controls, coef_col


def donor_expression(cp10k: sparse.csr_matrix, genes: list[str], gene_idx: list[int], metadata: pd.DataFrame, cells: list[str], donor_col: str, weights: np.ndarray | None = None) -> pd.DataFrame:
    meta = metadata.set_index("cell_id").reindex(cells)
    rows = []
    for donor, labels in meta.groupby(donor_col).groups.items():
        idx = np.asarray([cells.index(label) for label in labels], dtype=int)
        w = np.ones(len(idx), dtype=float) if weights is None else weights[idx]
        denom = float(w.sum())
        values = np.full(len(gene_idx), np.nan)
        if denom > 0:
            values = np.asarray(w @ cp10k[idx, :][:, gene_idx] / denom).ravel()
        row = {"donor_id": donor, "sum_state_weight": denom}
        row.update({gene: np.log1p(value) if np.isfinite(value) else np.nan for gene, value in zip(genes, values)})
        rows.append(row)
    return pd.DataFrame(rows).set_index("donor_id")


def add_fdr(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return frame
    frame["q_global"] = bh_fdr(frame["p_value"])
    frame["q_by_trait"] = frame.groupby("phenotype", group_keys=False)["p_value"].apply(bh_fdr)
    frame["q_by_gene"] = frame.groupby("gene", group_keys=False)["p_value"].apply(bh_fdr)
    return frame


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--raw-10x-dir", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--cell-state-activity", type=Path, required=True)
    ap.add_argument("--phenotypes", type=Path, required=True)
    ap.add_argument("--genes", type=Path, default=Path(""))
    ap.add_argument("--states", default="")
    ap.add_argument("--covariates", default="")
    ap.add_argument("--donor-col", default="donor_id")
    ap.add_argument("--phenotype", default="")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    counts, cells, all_genes = load_10x(args.raw_10x_dir)
    totals = np.asarray(counts.sum(axis=1)).ravel()
    cp10k = counts.multiply(10000.0 / totals[:, None]).tocsr()
    requested = set(read_one_column(args.genes)) if args.genes and str(args.genes) else set(all_genes)
    gene_idx = [i for i, g in enumerate(all_genes) if g in requested]
    genes = [all_genes[i] for i in gene_idx]
    metadata = read_table(args.metadata).drop_duplicates("cell_id")
    phenotypes = read_table(args.phenotypes).drop_duplicates(args.donor_col).set_index(args.donor_col)
    phenotype_cols = [args.phenotype] if args.phenotype else [c for c in phenotypes.columns if c not in set(args.covariates.split(","))]
    phenotype_cols = [c for c in phenotype_cols if c]
    covariates = [c for c in args.covariates.split(",") if c]
    activity = read_table(args.cell_state_activity)
    activity = activity.loc[activity["state_type"].eq("biological") & activity["cell_id"].isin(cells)].copy()
    if args.states:
        selected = {s.strip() for s in args.states.split(",") if s.strip()}
        activity = activity.loc[activity["state_name"].isin(selected)]
    cell_pos = pd.Series(np.arange(len(cells)), index=cells)

    rows_parent = []
    parent_expr = donor_expression(cp10k, genes, gene_idx, metadata, cells, args.donor_col)
    for phenotype in phenotype_cols:
        design = phenotypes[[phenotype] + [c for c in covariates if c in phenotypes.columns]]
        for gene in genes:
            coef, p, n, n_cases, n_controls, coef_col = fit_lm(parent_expr[gene], design.reindex(parent_expr.index), phenotype)
            rows_parent.append({"gene": gene, "state_name": "", "state_weight_type": "whole_parent", "phenotype": phenotype, "coefficient": coef, "coefficient_units": "log1p_cp10k_per_standardized_phenotype_or_category", "p_value": p, "n_donors": n, "n_cases": n_cases, "n_controls": n_controls, "model_formula": f"log1p_parent_cp10k ~ {phenotype}", "expression_unit_for_model": "log1p donor mean CP10K", "interpretation_scope": "whole_parent"})

    state_rows = {"gradient_percentile_squared": [], "high_tail_percentile_90_100": []}
    contrast_rows = []
    for state, group in activity.groupby("state_name", sort=False):
        positions = cell_pos.reindex(group["cell_id"]).dropna().astype(int).to_numpy()
        for weight_type, weight_col in [("gradient_percentile_squared", "state_activity_weight_gradient"), ("high_tail_percentile_90_100", "state_activity_weight_hightail")]:
            weights = np.zeros(len(cells), dtype=float)
            weights[positions] = pd.to_numeric(group[weight_col], errors="coerce").fillna(0.0).to_numpy(float)
            expr = donor_expression(cp10k, genes, gene_idx, metadata, cells, args.donor_col, weights)
            contrast = expr[genes] - parent_expr.reindex(expr.index)[genes]
            for phenotype in phenotype_cols:
                design = phenotypes[[phenotype] + [c for c in covariates if c in phenotypes.columns]]
                for gene in genes:
                    coef, p, n, n_cases, n_controls, _ = fit_lm(expr[gene], design.reindex(expr.index), phenotype)
                    row = {"gene": gene, "state_name": state, "state_weight_type": weight_type, "phenotype": phenotype, "coefficient": coef, "coefficient_units": "log1p_cp10k_per_standardized_phenotype_or_category", "p_value": p, "n_donors": n, "n_cases": n_cases, "n_controls": n_controls, "model_formula": f"log1p_state_weighted_cp10k ~ {phenotype}", "expression_unit_for_model": "log1p donor state-weighted mean CP10K", "interpretation_scope": "state_weighted"}
                    state_rows[weight_type].append(row)
                    ccoef, cp, cn, ccases, ccontrols, _ = fit_lm(contrast[gene], design.reindex(contrast.index), phenotype)
                    contrast_rows.append({**row, "coefficient": ccoef, "p_value": cp, "n_donors": cn, "n_cases": ccases, "n_controls": ccontrols, "model_formula": f"log1p_state_minus_parent_cp10k ~ {phenotype}", "interpretation_scope": "state_specific_contrast"})

    parent = add_fdr(pd.DataFrame(rows_parent))
    grad = add_fdr(pd.DataFrame(state_rows["gradient_percentile_squared"]))
    tail = add_fdr(pd.DataFrame(state_rows["high_tail_percentile_90_100"]))
    contrast = add_fdr(pd.DataFrame(contrast_rows))
    write_table(parent, args.out_dir / "de_whole_parent.tsv.gz")
    write_table(grad, args.out_dir / "de_state_weighted_gradient.tsv.gz")
    write_table(tail, args.out_dir / "de_state_weighted_hightail.tsv.gz")
    write_table(contrast, args.out_dir / "de_state_specific_contrast.tsv.gz")
    summary = {"n_genes": len(genes), "n_states": int(activity["state_name"].nunique()), "phenotypes": phenotype_cols, "covariates": covariates, "timestamp": datetime.now().isoformat(timespec="seconds")}
    (args.out_dir / "de_run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
