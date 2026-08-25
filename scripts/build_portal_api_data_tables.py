#!/usr/bin/env python3
"""Build compact portal data tables from cell-state expression and program matching outputs."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path
import json

import numpy as np
import pandas as pd
from pandas.errors import EmptyDataError


def read_table(path: str | Path | None) -> pd.DataFrame:
    if not path or str(path) in {"", ".", "None", "none", "NULL", "null"}:
        return pd.DataFrame()
    p = Path(path)
    if p.is_dir() or not p.exists() or p.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(p, sep="	", compression="infer", low_memory=False)
    except EmptyDataError:
        return pd.DataFrame()


def write_table(frame: pd.DataFrame, path: Path, columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    out = frame.copy()
    for col in columns:
        if col not in out.columns:
            out[col] = "" if col not in {"log10_cpk", "log2fc_weighted_vs_all_parent", "p_value", "beta", "beta_uncorrected", "importance", "value", "correlation", "gsea_p", "gsea_q"} else np.nan
    out[columns].to_csv(path, sep="\t", index=False, compression="infer")


def display_factor(value: str) -> str:
    text = str(value)
    m = re.search(r"(?:factor|Factor)[_ ]*([0-9]+)$", text)
    if m:
        return f"Factor{m.group(1)}"
    m = re.search(r"liger[_-]?factor[_-]?([0-9]+)", text, re.I)
    if m:
        return f"Factor{m.group(1)}"
    return text


def ensure_expression_fields(frame: pd.DataFrame) -> pd.DataFrame:
    out = frame.copy()
    if "annotated_cell_type" in out.columns and "cell_type" not in out.columns:
        out["cell_type"] = out["annotated_cell_type"]
    if "weighted_mean_cp10k" in out.columns:
        value = pd.to_numeric(out["weighted_mean_cp10k"], errors="coerce")
    elif "weighted_mean_expression" in out.columns:
        value = pd.to_numeric(out["weighted_mean_expression"], errors="coerce")
    elif "log10_cpk" in out.columns:
        value = (10 ** pd.to_numeric(out["log10_cpk"], errors="coerce")) - 1
    else:
        value = pd.Series(np.nan, index=out.index)
    out["log10_cpk"] = np.log10(value.fillna(0).clip(lower=0) + 1.0)
    if "log2fc_weighted_vs_all_parent" not in out.columns:
        out["log2fc_weighted_vs_all_parent"] = np.nan
    if "p_value" not in out.columns:
        out["p_value"] = np.nan
    return out




def filter_default_weight(frame: pd.DataFrame, default_weight_type: str) -> pd.DataFrame:
    if frame.empty or "state_weight_type" not in frame.columns:
        return frame
    filtered = frame.loc[frame["state_weight_type"].astype(str).eq(default_weight_type)].copy()
    return filtered if not filtered.empty else frame

def program_long_loadings(path: Path, dataset: str, model: str, cell_type: str) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame(columns=["dataset", "cell_type", "model", "factor", "gene", "value"])
    df = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    first = df.columns[0]
    long = df.melt(id_vars=[first], var_name="factor", value_name="value").rename(columns={first: "gene"})
    long["dataset"] = dataset
    long["cell_type"] = cell_type
    long["model"] = model
    long["factor"] = long["factor"].map(display_factor)
    return long[["dataset", "cell_type", "model", "factor", "gene", "value"]]


def infer_cell_type_from_gene_set(gene_set: pd.Series, tissue: str, is_program: bool) -> pd.Series:
    text = gene_set.astype(str)
    if is_program:
        stripped = text.str.replace(r"^" + re.escape(tissue) + r"_", "", regex=True)
        return stripped.str.replace(r"_program_factor_.*$", "", regex=True)
    stripped = text.str.replace(r"^" + re.escape(tissue) + r"_", "", regex=True)
    return stripped.str.rsplit("_", n=1).str[0]


def normalize_pigean(frame: pd.DataFrame, id_col: str, tissue: str, dataset: str, model: str, is_program: bool) -> pd.DataFrame:
    if frame.empty:
        cols = ["dataset", "cell_type", "model", "factor", "trait", "beta", "beta_uncorrected"] if is_program else ["tissue", "cell_type", "state_name", "trait", "beta", "beta_uncorrected"]
        return pd.DataFrame(columns=cols)
    out = frame.copy()
    lower = {str(c).lower(): c for c in out.columns}
    if id_col not in out.columns:
        for candidate in ["gene_set", "gene_set_name", "state_name", "name", "set"]:
            col = lower.get(candidate)
            if col is not None:
                out[id_col] = out[col]
                break
    if id_col not in out.columns:
        cols = ["dataset", "cell_type", "model", "factor", "trait", "beta", "beta_uncorrected"] if is_program else ["tissue", "cell_type", "state_name", "trait", "beta", "beta_uncorrected"]
        return pd.DataFrame(columns=cols)
    if "trait" not in out.columns:
        for candidate in ["trait", "phenotype", "trait_internal", "y", "trait_name"]:
            col = lower.get(candidate)
            if col is not None:
                out["trait"] = out[col]
                break
    if "trait" not in out.columns:
        out["trait"] = ""
    for target, candidates in {"beta": ["beta"], "beta_uncorrected": ["beta_uncorrected", "beta_uncorrected_orig"]}.items():
        if target not in out.columns:
            for candidate in candidates:
                col = lower.get(candidate)
                if col is not None:
                    out[target] = out[col]
                    break
        if target not in out.columns:
            out[target] = np.nan
        out[target] = pd.to_numeric(out[target], errors="coerce")
    if is_program:
        out["dataset"] = dataset
        out["model"] = model
        if "cell_type" not in out.columns:
            out["cell_type"] = infer_cell_type_from_gene_set(out[id_col], tissue, True)
        out["factor"] = out[id_col].map(display_factor)
        return out[["dataset", "cell_type", "model", "factor", "trait", "beta", "beta_uncorrected"]]
    out["tissue"] = out.get("tissue", tissue)
    if "cell_type" not in out.columns:
        out["cell_type"] = infer_cell_type_from_gene_set(out[id_col], tissue, False)
    out["state_name"] = out[id_col]
    return out[["tissue", "cell_type", "state_name", "trait", "beta", "beta_uncorrected"]]

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--dataset", default="")
    ap.add_argument("--model", default="")
    ap.add_argument("--cell-state-expression", type=Path, default=None)
    ap.add_argument("--program-expression", type=Path, default=None)
    ap.add_argument("--cell-type-expression", type=Path, default=None)
    ap.add_argument("--program-loadings-manifest", type=Path, default=None, help="TSV with cell_type and program_loadings columns")
    ap.add_argument("--program-match-dir", type=Path, default=None)
    ap.add_argument("--cell-state-pigean", type=Path, default=None)
    ap.add_argument("--program-pigean", type=Path, default=None)
    ap.add_argument("--program-factors", type=Path, default=None)
    ap.add_argument("--default-state-weight-type", default="gradient_percentile_squared")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    state_expr = ensure_expression_fields(filter_default_weight(read_table(args.cell_state_expression), args.default_state_weight_type))
    if not state_expr.empty:
        state_expr["tissue"] = state_expr.get("tissue", args.tissue)
    write_table(state_expr, args.out_dir / "cell_state_expression.tsv.gz", ["gene", "tissue", "cell_type", "state_name", "log10_cpk", "log2fc_weighted_vs_all_parent", "p_value"])

    cell_type_expr = ensure_expression_fields(read_table(args.cell_type_expression))
    if not cell_type_expr.empty:
        cell_type_expr["tissue"] = cell_type_expr.get("tissue", args.tissue)
    write_table(cell_type_expr, args.out_dir / "cell_type_expression.tsv.gz", ["gene", "tissue", "cell_type", "log10_cpk", "log2fc_weighted_vs_all_parent", "p_value"])

    program_expr = ensure_expression_fields(filter_default_weight(read_table(args.program_expression), args.default_state_weight_type))
    if not program_expr.empty:
        program_expr["dataset"] = args.dataset
        program_expr["model"] = args.model
        if "factor" not in program_expr.columns:
            program_expr["factor"] = program_expr.get("state_name", "").map(display_factor)
        else:
            program_expr["factor"] = program_expr["factor"].map(display_factor)
    write_table(program_expr, args.out_dir / "program_expression.tsv.gz", ["dataset", "cell_type", "model", "factor", "gene", "log10_cpk", "log2fc_weighted_vs_all_parent", "p_value"])

    loading_frames = []
    manifest = read_table(args.program_loadings_manifest)
    if not manifest.empty:
        for row in manifest.itertuples(index=False):
            path = Path(getattr(row, "program_loadings"))
            cell_type = str(getattr(row, "cell_type"))
            loading_frames.append(program_long_loadings(path, args.dataset, args.model, cell_type))
    loadings = pd.concat(loading_frames, ignore_index=True) if loading_frames else pd.DataFrame()
    write_table(loadings, args.out_dir / "program_gene_loadings.tsv.gz", ["dataset", "cell_type", "model", "factor", "gene", "value"])
    top_gene_lookup = {}
    if not loadings.empty:
        ranked = loadings.sort_values(["cell_type", "factor", "value"], ascending=[True, True, False])
        top_gene_lookup = ranked.groupby(["cell_type", "factor"])["gene"].apply(lambda x: ";".join(x.head(20).astype(str))).to_dict()

    heat_frames = []
    label_frames = []
    if args.program_match_dir and args.program_match_dir.exists():
        for p in sorted(args.program_match_dir.glob("*/program_state_heatmap_long.tsv.gz")):
            heat_frames.append(read_table(p))
        for p in sorted(args.program_match_dir.glob("*/program_label_suggestions.tsv.gz")):
            labels = read_table(p)
            cell_type = p.parent.name
            if not labels.empty:
                labels["cell_type"] = cell_type
            label_frames.append(labels)
    factor_labels = {}
    for p in sorted(args.program_match_dir.glob("*/mouse_msigdb/factors.json")):
        with open(p, 'r') as f:
            for line in f:
                dict_line = json.loads(line.strip())
                factor_labels[display_factor(dict_line['factor'])] = dict_line
    heat = pd.concat(heat_frames, ignore_index=True) if heat_frames else pd.DataFrame()
    if not heat.empty:
        if "state_name" not in heat.columns and "state_id" in heat.columns:
            heat["state_name"] = heat["state_id"]
        heat["tissue"] = heat.get("tissue", args.tissue)
    write_table(heat, args.out_dir / "program_state_heatmap.tsv.gz", ["tissue", "cell_type", "state_name", "program_id", "correlation", "gsea_p", "gsea_q"])

    labels = pd.concat(label_frames, ignore_index=True) if label_frames else pd.DataFrame()
    if labels.empty and not loadings.empty:
        labels = loadings[["cell_type", "factor"]].drop_duplicates().copy()
        labels["program_id"] = labels["factor"]
        labels["suggested_program_label"] = labels["factor"]
        labels["best_curated_state_id"] = ""
        labels["best_qc_state_id"] = ""
    if not labels.empty:
        labels["dataset"] = args.dataset
        labels["model"] = args.model
        labels["factor"] = labels["program_id"].map(display_factor)
        labels["label"] = labels["factor"].map(lambda x: factor_labels.get(x['label'], ""))
        labels["quality"] = labels.get("suggested_program_quality_class", "")
        labels["importance"] = np.nan
        labels["top_cells"] = ""
        labels["top_gene_sets"] = ""
        labels["top_genes"] = [top_gene_lookup.get((str(ct), str(factor)), "") for ct, factor in zip(labels["cell_type"], labels["factor"])]
        labels["top_traits"] = ""
        labels["significant_cell_states"] = labels.get("best_curated_state_id", "")
        labels["qc_cell_states"] = labels.get("best_qc_state_id", "")
    write_table(labels, args.out_dir / "program_factor_metadata.tsv.gz", ["dataset", "cell_type", "model", "factor", "importance", "label", "top_cells", "top_gene_sets", "top_genes", "top_traits", "significant_cell_states", "qc_cell_states"])

    state_pigean = normalize_pigean(read_table(args.cell_state_pigean), "state_name", args.tissue, args.dataset, args.model, False)
    write_table(state_pigean, args.out_dir / "cell_state_pigean_trait_results.tsv.gz", ["tissue", "cell_type", "state_name", "trait", "beta", "beta_uncorrected"])
    program_pigean = normalize_pigean(read_table(args.program_pigean), "factor", args.tissue, args.dataset, args.model, True)
    write_table(program_pigean, args.out_dir / "program_pigean_trait_results.tsv.gz", ["dataset", "cell_type", "model", "factor", "trait", "beta", "beta_uncorrected"])


if __name__ == "__main__":
    main()
