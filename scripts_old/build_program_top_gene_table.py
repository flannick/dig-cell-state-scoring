#!/usr/bin/env python3
"""Build a cell_type / factor / top genes table from program loading files."""

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


def read_program_map(path: Path | None) -> dict[str, str]:
    if not path or not path.exists() or path.stat().st_size == 0:
        return {}
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "program_dir" not in frame.columns or "cell_type" not in frame.columns:
        raise SystemExit(f"Program map must contain program_dir and cell_type columns: {path}")
    return dict(zip(frame["program_dir"], frame["cell_type"]))


def long_from_loadings(path: Path) -> pd.DataFrame:
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

    # LIGER exports often have a header with only factor names and no explicit
    # gene-column label. In that case data rows have one extra field. Parse this
    # path line-by-line because pandas treats the ragged header as an error.
    with path.open("r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    header = lines[0].split("\t") if lines else []
    first_body = lines[1].split("\t") if len(lines) > 1 else []
    if header and first_body and len(first_body) == len(header) + 1:
        rows = []
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) < 2:
                continue
            rows.append(parts[: len(header) + 1])
        wide = pd.DataFrame(rows, columns=["gene"] + header)
    else:
        gene_col = frame.columns[0]
        wide = frame.rename(columns={gene_col: "gene"})
    long = wide.melt(id_vars=["gene"], var_name="factor", value_name="loading")
    long["loading"] = pd.to_numeric(long["loading"], errors="coerce")
    return long.dropna(subset=["gene", "factor", "loading"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--program-root", type=Path, required=True, help="Directory containing per-cell-type program subdirectories.")
    ap.add_argument("--program-cell-type-map", type=Path, default=None, help="Optional TSV with program_dir and cell_type columns.")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--min-loading", type=float, default=None, help="Optional minimum loading for top-gene inclusion.")
    args = ap.parse_args()

    program_map = read_program_map(args.program_cell_type_map)
    rows = []
    loading_files = sorted(args.program_root.glob("*/gene_loadings.tsv"))
    if not loading_files:
        raise SystemExit(f"No */gene_loadings.tsv files found under {args.program_root}")

    for loadings_path in loading_files:
        program_dir = loadings_path.parent.name
        cell_type = program_map.get(program_dir, snake_id(program_dir))
        long = long_from_loadings(loadings_path)
        if args.min_loading is not None:
            long = long[long["loading"].ge(args.min_loading)]
        for factor, group in long.groupby("factor", sort=True):
            top = group.sort_values(["loading", "gene"], ascending=[False, True]).head(args.top_n).copy()
            if top.empty:
                continue
            genes = top["gene"].astype(str).tolist()
            weights = [f"{x:.6g}" for x in top["loading"].astype(float).tolist()]
            rows.append({
                "cell_type": cell_type,
                "program_dir": program_dir,
                "factor": display_factor(str(factor)),
                "factor_raw": str(factor),
                "top_n": len(top),
                "top_genes": ";".join(genes),
                "top_gene_loadings": ";".join(weights),
                "top_genes_with_loadings": ";".join(f"{g}:{w}" for g, w in zip(genes, weights)),
                "program_loadings": str(loadings_path),
            })

    out = pd.DataFrame(rows).sort_values(["cell_type", "factor_raw"])
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False, compression="infer")
    print(f"Wrote {len(out)} program-factor rows to {args.out}")


if __name__ == "__main__":
    main()
