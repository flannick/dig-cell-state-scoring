#!/usr/bin/env python3
"""Split a tissue/cell-type/state GMT into tissue and cell-type GMT files.

State names are expected to begin with a tissue token followed by a parent
cell-type token, for example ``tissue_a_parent_cell_type_state_label``. The
splitter uses metadata tissue and cell-type labels as the vocabulary and
chooses the longest normalized cell-type prefix that matches each state name
after the tissue prefix.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def norm_name(value: str) -> str:
    value = str(value).strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return re.sub(r"_+", "_", value).strip("_")


def read_gmt(path: Path) -> list[list[str]]:
    rows: list[list[str]] = []
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                rows.append(parts)
    if not rows:
        raise SystemExit(f"No GMT rows found in {path}")
    return rows


def infer_scope(state_name: str, tissues: list[str], cell_types: list[str]) -> tuple[str, str]:
    norm_state = norm_name(state_name)
    tissue_match = ""
    remainder = norm_state
    for tissue in sorted(tissues, key=len, reverse=True):
        if norm_state.startswith(tissue + "_"):
            tissue_match = tissue
            remainder = norm_state[len(tissue) + 1 :]
            break
    if not tissue_match:
        parts = norm_state.split("_", 1)
        tissue_match = parts[0]
        remainder = parts[1] if len(parts) > 1 else ""
    cell_match = ""
    for cell_type in sorted(cell_types, key=len, reverse=True):
        if remainder.startswith(cell_type + "_"):
            cell_match = cell_type
            break
    return tissue_match, cell_match


def write_gmt(rows: list[list[str]], path: Path) -> None:
    with path.open("w", encoding="utf-8") as out:
        for row in rows:
            out.write("\t".join(row) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gmt", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--tissue-col", default="tissue")
    ap.add_argument("--cell-type-col", default="annotated_cell_type")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--unmatched-out", type=Path, default=Path(""))
    args = ap.parse_args()

    metadata = pd.read_csv(args.metadata, sep="\t", compression="infer", low_memory=False)
    if args.tissue_col not in metadata.columns:
        raise SystemExit(f"Metadata is missing tissue column: {args.tissue_col}")
    if args.cell_type_col not in metadata.columns:
        raise SystemExit(f"Metadata is missing cell type column: {args.cell_type_col}")
    tissues = sorted({norm_name(x) for x in metadata[args.tissue_col].dropna().unique()})
    cell_types = sorted({norm_name(x) for x in metadata[args.cell_type_col].dropna().unique()})
    rows = read_gmt(args.gmt)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    grouped: dict[tuple[str, str], list[list[str]]] = {}
    summary = []
    unmatched = []
    for row in rows:
        tissue, cell_type = infer_scope(row[0], tissues, cell_types)
        summary.append({"state_name": row[0], "tissue": tissue, "cell_type": cell_type, "n_genes": max(len(row) - 2, 0)})
        if not cell_type:
            unmatched.append(row)
            continue
        grouped.setdefault((tissue, cell_type), []).append(row)

    out_rows = []
    for (tissue, cell_type), group_rows in sorted(grouped.items()):
        tissue_dir = args.out_dir / tissue
        tissue_dir.mkdir(parents=True, exist_ok=True)
        path = tissue_dir / f"{cell_type}.gmt"
        write_gmt(group_rows, path)
        out_rows.append({"tissue": tissue, "cell_type": cell_type, "n_states": len(group_rows), "gmt": str(path)})
    if unmatched:
        unmatched_path = args.unmatched_out or (args.out_dir / "unmatched.gmt")
        write_gmt(unmatched, unmatched_path)
    pd.DataFrame(summary).to_csv(args.out_dir / "state_scope_summary.tsv", sep="\t", index=False)
    pd.DataFrame(out_rows).to_csv(args.out_dir / "split_gmt_manifest.tsv", sep="\t", index=False)
    print(f"Wrote {len(out_rows)} split GMT files to {args.out_dir}")
    if unmatched:
        print(f"Warning: {len(unmatched)} GMT rows did not match a metadata cell type")


if __name__ == "__main__":
    main()
