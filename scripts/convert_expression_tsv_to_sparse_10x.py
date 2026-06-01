#!/usr/bin/env python3
"""Convert a dense/wide expression TSV into a sparse 10x-like Matrix Market directory."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
from datetime import datetime
from pathlib import Path

import numpy as np

from matrix_value_types import infer_matrix_value_type, resolve_value_type


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8" if "t" in mode else None)
    return open(path, mode, encoding="utf-8" if "t" in mode else None)


def write_gzip_lines(path: Path, lines) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        for line in lines:
            handle.write(line)


def infer_orientation(first_col: str) -> str:
    token = first_col.strip().lower()
    if token in {"gene", "genes", "gene_id", "gene_symbol", "symbol", "feature", "features"}:
        return "gene_by_cell"
    if token in {"cell", "cell_id", "barcode", "barcodes", "id"}:
        return "cell_by_gene"
    return "gene_by_cell"


def parse_float(value: str) -> float:
    if value == "" or value.upper() == "NA" or value.lower() == "nan":
        return 0.0
    return float(value)


def read_include_ids(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    ids: set[str] = set()
    with open_text(path, "rt") as handle:
        for line in handle:
            value = line.strip().split("	")[0]
            if value:
                ids.add(value)
    return ids


def choose_indices(labels: list[str], include_ids: set[str] | None, sample_fraction: float, sample_seed: int) -> list[int]:
    if include_ids is not None:
        indices = [i for i, label in enumerate(labels) if label in include_ids]
    else:
        indices = list(range(len(labels)))
    if sample_fraction < 1.0:
        if sample_fraction <= 0:
            raise SystemExit("--cell-sample-fraction must be > 0")
        rng = np.random.default_rng(sample_seed)
        keep = rng.random(len(indices)) < sample_fraction
        indices = [idx for idx, flag in zip(indices, keep) if bool(flag)]
        if not indices and labels:
            indices = [int(rng.integers(0, len(labels)))]
    return indices


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--matrix-tsv", type=Path, required=True, help="Dense/wide TSV or TSV.GZ expression matrix")
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--orientation", choices=["auto", "gene_by_cell", "cell_by_gene"], default="auto")
    ap.add_argument("--value-type", choices=["auto", "raw_counts", "linear_cp10k", "log1p_cp10k", "linear_normalized", "log1p_normalized", "scaled"], default="auto")
    ap.add_argument("--zero-epsilon", type=float, default=0.0, help="Treat absolute values at or below this threshold as zero")
    ap.add_argument("--max-infer-values", type=int, default=200000)
    ap.add_argument("--progress-every-rows", type=int, default=1000)
    ap.add_argument("--cell-include", type=Path, default=None, help="Optional one-column list of cells to keep during conversion")
    ap.add_argument("--cell-sample-fraction", type=float, default=1.0, help="Random fraction of cells to keep before writing sparse output")
    ap.add_argument("--cell-sample-seed", type=int, default=1)
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix="matrix_entries_", dir=str(args.out_dir)))
    entries_path = tmp / "entries.tsv"
    sample_values: list[float] = []
    n_rows = 0
    n_cols = 0
    nnz = 0
    row_names: list[str] = []
    col_names: list[str] = []
    original_n_cells = 0
    original_n_genes = 0
    kept_label_indices: list[int] = []

    try:
        with open_text(args.matrix_tsv, "rt") as handle, open(entries_path, "w", encoding="utf-8") as entries:
            header = handle.readline().rstrip("\n").split("\t")
            if len(header) < 2:
                raise SystemExit("Expression TSV must have at least one ID column and one value column")
            orientation = infer_orientation(header[0]) if args.orientation == "auto" else args.orientation
            labels = header[1:]
            include_ids = read_include_ids(args.cell_include)
            if orientation == "gene_by_cell":
                original_n_cells = len(labels)
                kept_label_indices = choose_indices(labels, include_ids, args.cell_sample_fraction, args.cell_sample_seed)
                col_names = [labels[i] for i in kept_label_indices]
                cell_index_map = {old_i: new_i + 1 for new_i, old_i in enumerate(kept_label_indices)}
            else:
                original_n_genes = len(labels)
                row_names = labels
                kept_row_counter = 0
            for line in handle:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                feature = parts[0]
                values = parts[1:]
                if len(values) != len(labels):
                    raise SystemExit(f"Row {n_rows + 2} has {len(values)} values but expected {len(labels)}")
                n_rows += 1
                if orientation == "gene_by_cell":
                    row_names.append(feature)
                    gene_idx = n_rows
                    for old_cell_zero in kept_label_indices:
                        value = parse_float(values[old_cell_zero])
                        if abs(value) <= args.zero_epsilon:
                            continue
                        nnz += 1
                        if len(sample_values) < args.max_infer_values:
                            sample_values.append(value)
                        entries.write(f"{gene_idx}\t{cell_index_map[old_cell_zero]}\t{value:.12g}\n")
                else:
                    keep_cell = include_ids is None or feature in include_ids
                    if keep_cell and args.cell_sample_fraction < 1.0:
                        row_rng = np.random.default_rng(args.cell_sample_seed + n_rows)
                        keep_cell = bool(row_rng.random() < args.cell_sample_fraction)
                    if not keep_cell:
                        continue
                    kept_row_counter += 1
                    col_names.append(feature)
                    cell_idx = kept_row_counter
                    for gene_zero, raw_value in enumerate(values):
                        value = parse_float(raw_value)
                        if abs(value) <= args.zero_epsilon:
                            continue
                        nnz += 1
                        if len(sample_values) < args.max_infer_values:
                            sample_values.append(value)
                        entries.write(f"{gene_zero + 1}\t{cell_idx}\t{value:.12g}\n")
                if args.progress_every_rows and n_rows % args.progress_every_rows == 0:
                    print(f"[convert-expression] processed {n_rows} rows; nnz={nnz}", flush=True)
        if orientation == "gene_by_cell":
            genes = row_names
            cells = col_names
            original_n_genes = len(row_names)
        else:
            genes = row_names
            cells = col_names
            original_n_cells = n_rows
        n_genes = len(genes)
        n_cells = len(cells)
        if n_genes == 0 or n_cells == 0:
            raise SystemExit("No genes or cells found in expression TSV")
        sample_matrix = np.asarray(sample_values, dtype=float)
        inferred = infer_matrix_value_type(sample_matrix, max_values=args.max_infer_values)
        if args.value_type == "auto":
            value_type = inferred["inferred_value_type"]
            if inferred["confidence"] == "low":
                raise SystemExit(f"Could not confidently infer matrix value type: {inferred}")
            if value_type == "scaled":
                raise SystemExit("Matrix appears to contain scaled/negative values; sparse rank inputs must be nonnegative")
        else:
            value_type = args.value_type
            # Validate explicit type against sampled values using shared logic.
            value_type, inferred = resolve_value_type(args.value_type, sample_matrix, context="matrix")
        matrix_path = args.out_dir / "matrix.mtx.gz"
        with gzip.open(matrix_path, "wt", encoding="utf-8") as out, open(entries_path, "r", encoding="utf-8") as entries:
            out.write("%%MatrixMarket matrix coordinate real general\n")
            out.write("% generated by convert_expression_tsv_to_sparse_10x.py\n")
            out.write(f"{n_genes} {n_cells} {nnz}\n")
            shutil.copyfileobj(entries, out)
        write_gzip_lines(args.out_dir / "features.tsv.gz", (f"{gene}\t{gene}\tGene Expression\n" for gene in genes))
        write_gzip_lines(args.out_dir / "barcodes.tsv.gz", (f"{cell}\n" for cell in cells))
        report = {
            "input_matrix": str(args.matrix_tsv),
            "orientation": orientation,
            "matrix_value_type": value_type,
            "inference": inferred,
            "n_genes": n_genes,
            "n_cells": n_cells,
            "original_n_genes": original_n_genes,
            "original_n_cells": original_n_cells,
            "cell_sample_fraction": args.cell_sample_fraction,
            "cell_sample_seed": args.cell_sample_seed,
            "cell_include": str(args.cell_include) if args.cell_include else None,
            "nnz": nnz,
            "zero_epsilon": args.zero_epsilon,
            "timestamp": datetime.now().isoformat(timespec="seconds"),
        }
        (args.out_dir / "matrix_value_type.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        (args.out_dir / "value_type_inference_report.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
        print(f"Wrote {args.out_dir} with {n_genes} genes, {n_cells} cells, {nnz} nonzero values, value_type={value_type}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
