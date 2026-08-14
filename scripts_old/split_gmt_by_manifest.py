#!/usr/bin/env python3
"""Split a state GMT into tissue/cell-type GMTs using a state manifest."""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import pandas as pd


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def safe_label(value: object) -> str:
    value = str(value).strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return re.sub(r"_+", "_", value).strip("_") or "unknown"


def read_gmt(path: Path) -> dict[str, list[str]]:
    rows: dict[str, list[str]] = {}
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                rows[parts[0]] = parts
    if not rows:
        raise SystemExit(f"No GMT rows found in {path}")
    return rows


def write_gmt(rows: list[list[str]], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as out:
        for row in rows:
            out.write("\t".join(row) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gmt", type=Path, required=True)
    ap.add_argument("--state-manifest", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--tissue-col", default="tissue")
    ap.add_argument("--cell-type-col", default="cell_type")
    ap.add_argument("--require-all-states", action="store_true")
    args = ap.parse_args()

    gmt = read_gmt(args.gmt)
    manifest = pd.read_csv(args.state_manifest, sep="\t", compression="infer", low_memory=False)
    required = {"state_name", args.tissue_col, args.cell_type_col}
    missing = required - set(manifest.columns)
    if missing:
        raise SystemExit(f"State manifest is missing required column(s): {', '.join(sorted(missing))}")
    manifest = manifest.loc[manifest["state_name"].isin(gmt)].copy()
    missing_states = sorted(set(gmt) - set(manifest["state_name"]))
    if missing_states and args.require_all_states:
        preview = ", ".join(missing_states[:10])
        suffix = "..." if len(missing_states) > 10 else ""
        raise SystemExit(f"Manifest missing {len(missing_states)} GMT state(s): {preview}{suffix}")

    grouped: dict[tuple[str, str], list[list[str]]] = {}
    for _, row in manifest.iterrows():
        key = (safe_label(row[args.tissue_col]), safe_label(row[args.cell_type_col]))
        grouped.setdefault(key, []).append(gmt[str(row["state_name"])])

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_rows = []
    for (tissue, cell_type), rows in sorted(grouped.items()):
        path = args.out_dir / tissue / f"{cell_type}.gmt"
        write_gmt(rows, path)
        out_rows.append({"tissue": tissue, "cell_type": cell_type, "n_states": len(rows), "gmt": str(path)})
    pd.DataFrame(out_rows).to_csv(args.out_dir / "split_gmt_manifest.tsv", sep="\t", index=False)
    if missing_states:
        pd.DataFrame({"state_name": missing_states}).to_csv(args.out_dir / "manifest_missing_states.tsv", sep="\t", index=False)
    print(f"Wrote {len(out_rows)} manifest-split GMT files to {args.out_dir}")


if __name__ == "__main__":
    main()
