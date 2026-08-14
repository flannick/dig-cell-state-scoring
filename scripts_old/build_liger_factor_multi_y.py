#!/usr/bin/env python3
"""Build a PIGEAN multi-y input from LIGER factor gene loadings.

Each tissue/cell-type/factor becomes one multi-y trait. The top positive-loading
factor genes are assigned 0-1 weights by dividing by the maximum loading within
that selected factor gene list. The weight is written to Direct, Combined, and
Indirect so PIGEAN can use the standard multi-y gene-stat columns.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


def snake_id(value: object) -> str:
    text = str(value).strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return re.sub(r"_+", "_", text).strip("_")


def display_factor(value: str) -> str:
    text = str(value)
    match = re.search(r"(?:program[_-]?)?factor[_-]?([0-9]+)$", text, re.I)
    if match:
        return f"Factor{match.group(1)}"
    match = re.search(r"(?:factor|Factor)[_ ]*([0-9]+)$", text)
    if match:
        return f"Factor{match.group(1)}"
    return text


def read_loadings_long(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", compression="infer")
    if frame.empty:
        return pd.DataFrame(columns=["gene", "factor", "loading"])
    columns_lower = {str(c).lower(): c for c in frame.columns}
    program_col = next((columns_lower[c] for c in ["program_id", "program", "factor", "module", "component"] if c in columns_lower), None)
    gene_col = next((columns_lower[c] for c in ["gene", "gene_symbol", "gene_name"] if c in columns_lower), None)
    loading_col = next((columns_lower[c] for c in ["loading", "weight", "value", "loading_score"] if c in columns_lower), None)
    if program_col and gene_col and loading_col:
        out = frame[[gene_col, program_col, loading_col]].copy()
        out.columns = ["gene", "factor", "loading"]
        out["loading"] = pd.to_numeric(out["loading"], errors="coerce")
        return out.dropna(subset=["gene", "factor", "loading"])

    with path.open("r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    header = lines[0].split("\t") if lines else []
    first_body = lines[1].split("\t") if len(lines) > 1 else []
    if header and first_body and len(first_body) == len(header) + 1:
        rows = []
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) >= len(header) + 1:
                rows.append(parts[: len(header) + 1])
        wide = pd.DataFrame(rows, columns=["gene"] + header)
    else:
        gene_col = frame.columns[0]
        wide = frame.rename(columns={gene_col: "gene"})
    long = wide.melt(id_vars=["gene"], var_name="factor", value_name="loading")
    long["loading"] = pd.to_numeric(long["loading"], errors="coerce")
    return long.dropna(subset=["gene", "factor", "loading"])


def trait_id(tissue: str, cell_type_dir: str, factor: str) -> str:
    return f"{snake_id(tissue)}__{snake_id(cell_type_dir)}__{snake_id(factor)}"


def factor_rows(group: pd.DataFrame, *, top_n: int, min_loading: float, weight_floor: float, weight_ceiling: float) -> pd.DataFrame:
    sub = group.copy()
    sub["loading"] = pd.to_numeric(sub["loading"], errors="coerce")
    sub = sub.dropna(subset=["gene", "loading"])
    sub = sub[sub["loading"].gt(min_loading)]
    sub = sub.sort_values(["loading", "gene"], ascending=[False, True])
    if top_n > 0:
        sub = sub.head(top_n)
    sub = sub.drop_duplicates(subset=["gene"], keep="first")
    if sub.empty:
        return pd.DataFrame(columns=["gene", "loading", "weight", "rank"])
    max_loading = float(sub["loading"].max())
    if max_loading <= 0:
        return pd.DataFrame(columns=["gene", "loading", "weight", "rank"])
    sub = sub.copy()
    sub["weight"] = (sub["loading"] / max_loading).clip(weight_floor, weight_ceiling)
    sub["rank"] = range(1, len(sub) + 1)
    return sub[["gene", "loading", "weight", "rank"]]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--liger-root", type=Path, default=Path("data/external/liger"))
    ap.add_argument("--out", type=Path, required=True, help="Output PIGEAN multi-y TSV/TSV.GZ")
    ap.add_argument("--manifest-out", type=Path, required=True)
    ap.add_argument("--top-n", type=int, default=500, help="Top positive-loading genes per factor; use 0 for all positive loadings.")
    ap.add_argument("--min-loading", type=float, default=0.0)
    ap.add_argument("--weight-floor", type=float, default=0.0)
    ap.add_argument("--weight-ceiling", type=float, default=1.0)
    ap.add_argument("--min-genes", type=int, default=5)
    args = ap.parse_args()

    loading_files = sorted(args.liger_root.glob("*/*/gene_loadings.tsv"))
    if not loading_files:
        raise SystemExit(f"No */*/gene_loadings.tsv files found under {args.liger_root}")

    y_rows = []
    manifest_rows = []
    for loading_path in loading_files:
        tissue = loading_path.parent.parent.name
        cell_type_dir = loading_path.parent.name
        long = read_loadings_long(loading_path)
        for factor_raw, group in long.groupby("factor", sort=True):
            factor_label = display_factor(str(factor_raw))
            rows = factor_rows(
                group,
                top_n=args.top_n,
                min_loading=args.min_loading,
                weight_floor=args.weight_floor,
                weight_ceiling=args.weight_ceiling,
            )
            tid = trait_id(tissue, cell_type_dir, factor_label)
            status = "ok" if len(rows) >= args.min_genes else "skipped_too_few_genes"
            if status == "ok":
                for row in rows.itertuples(index=False):
                    y_rows.append({
                        "Gene": str(row.gene),
                        "Trait_Internal": tid,
                        "Direct": float(row.weight),
                        "Combined": float(row.weight),
                        "Indirect": float(row.weight),
                    })
            top = rows.head(20)
            manifest_rows.append({
                "trait": tid,
                "tissue": tissue,
                "cell_type_dir": cell_type_dir,
                "cell_type_id": snake_id(cell_type_dir),
                "factor": factor_label,
                "factor_raw": str(factor_raw),
                "n_genes": len(rows),
                "top_n_requested": args.top_n,
                "min_loading": args.min_loading,
                "status": status,
                "gene_loadings": str(loading_path),
                "top_genes": ";".join(top["gene"].astype(str).tolist()),
                "top_gene_weights": ";".join(f"{x:.6g}" for x in top["weight"].astype(float).tolist()),
                "top_gene_loadings": ";".join(f"{x:.6g}" for x in top["loading"].astype(float).tolist()),
            })

    y = pd.DataFrame(y_rows, columns=["Gene", "Trait_Internal", "Direct", "Combined", "Indirect"])
    manifest = pd.DataFrame(manifest_rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.manifest_out.parent.mkdir(parents=True, exist_ok=True)
    y.to_csv(args.out, sep="\t", index=False, compression="infer")
    manifest.to_csv(args.manifest_out, sep="\t", index=False, compression="infer")
    n_traits = int(manifest[manifest["status"].eq("ok")]["trait"].nunique()) if not manifest.empty else 0
    print(f"Found {len(loading_files)} loading files")
    print(f"Wrote {len(y)} multi-y rows for {n_traits} traits to {args.out}")
    print(f"Wrote {len(manifest)} manifest rows to {args.manifest_out}")


if __name__ == "__main__":
    main()
