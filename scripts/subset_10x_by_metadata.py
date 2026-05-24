#!/usr/bin/env python3
"""Subset a 10x-style sparse Matrix Market directory by metadata filters.

This streams the Matrix Market coordinate file and writes a new 10x directory
without materializing the full sparse matrix in memory. It is intended for
large rank-universe matrices where downstream scoring should be run
incrementally by tissue, cell type, or another metadata-defined group.
"""

from __future__ import annotations

import argparse
import gzip
import re
import shutil
import tempfile
from pathlib import Path

import pandas as pd


def open_text(path: Path, mode: str = "rt"):
    return gzip.open(path, mode, encoding="utf-8") if path.suffix == ".gz" else open(path, mode, encoding="utf-8")


def resolve_10x_path(directory: Path, names: list[str]) -> Path:
    for name in names:
        path = directory / name
        if path.exists():
            return path
    raise SystemExit(f"Could not find any of {', '.join(names)} in {directory}")


def read_one_column(path: Path) -> list[str]:
    values: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            value = line.rstrip("\n").split("\t")[0]
            if value:
                values.append(value)
    return values


def parse_filter(value: str) -> tuple[str, set[str]]:
    if "=" not in value:
        raise SystemExit("--metadata-filter must have the form column=value1,value2")
    column, raw_values = value.split("=", 1)
    values = {item.strip() for item in raw_values.split(",") if item.strip()}
    if not column or not values:
        raise SystemExit("--metadata-filter must include a column and at least one value")
    return column, values


def safe_label(value: str) -> str:
    value = value.strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return re.sub(r"_+", "_", value).strip("_") or "subset"


def read_selected_cells(metadata_path: Path, cell_col: str, filters: list[str]) -> list[str]:
    metadata = read_filtered_metadata(metadata_path, cell_col, filters)
    cells = metadata[cell_col].drop_duplicates().astype(str).tolist()
    if not cells:
        raise SystemExit("Metadata filters selected zero cells")
    return cells


def read_filtered_metadata(metadata_path: Path, cell_col: str, filters: list[str]) -> pd.DataFrame:
    metadata = pd.read_csv(metadata_path, sep="\t", compression="infer", low_memory=False, dtype=str).fillna("")
    if cell_col not in metadata.columns:
        raise SystemExit(f"Metadata is missing cell ID column: {cell_col}")
    keep = pd.Series(True, index=metadata.index)
    for filter_text in filters:
        column, values = parse_filter(filter_text)
        if column not in metadata.columns:
            raise SystemExit(f"Metadata is missing filter column: {column}")
        keep &= metadata[column].isin(values)
    out = metadata.loc[keep].copy()
    if out.empty:
        raise SystemExit("Metadata filters selected zero cells")
    return out


def copy_text_gz(src: Path, dest: Path) -> None:
    with open_text(src) as inp, gzip.open(dest, "wt", encoding="utf-8") as out:
        shutil.copyfileobj(inp, out)


def subset_barcodes(src: Path, dest: Path, selected_positions: set[int]) -> dict[int, int]:
    old_to_new: dict[int, int] = {}
    new_index = 0
    with open_text(src) as inp, gzip.open(dest, "wt", encoding="utf-8") as out:
        for old_index, line in enumerate(inp, start=1):
            if old_index in selected_positions:
                new_index += 1
                old_to_new[old_index] = new_index
                out.write(line)
    return old_to_new


def subset_cell_totals(src: Path, dest: Path, selected_cells: set[str]) -> None:
    with open_text(src) as inp, gzip.open(dest, "wt", encoding="utf-8") as out:
        header = inp.readline()
        if header:
            out.write(header)
        for line in inp:
            cell = line.rstrip("\n").split("\t", 1)[0]
            if cell in selected_cells:
                out.write(line)


def write_subset_matrix(
    src: Path,
    dest: Path,
    n_genes: int,
    n_cells: int,
    old_to_new_cell: dict[int, int],
) -> int:
    nnz = 0
    with tempfile.NamedTemporaryFile(prefix="subset_10x_body_", suffix=".gz", delete=False) as tmp:
        tmp_body = Path(tmp.name)
    try:
        with open_text(src) as inp, gzip.open(tmp_body, "wt", encoding="utf-8") as body:
            banner = inp.readline()
            if not banner.startswith("%%MatrixMarket"):
                raise SystemExit(f"{src} is not a Matrix Market file")
            output_rows = n_genes
            output_cols = n_cells
            comments: list[str] = []
            for line in inp:
                if line.startswith("%"):
                    comments.append(line)
                    continue
                n_rows, n_cols, _ = [int(x) for x in line.split()[:3]]
                genes_by_cells = n_rows == n_genes and n_cols >= n_cells
                cells_by_genes = n_cols == n_genes and n_rows >= n_cells
                if not genes_by_cells and not cells_by_genes:
                    raise SystemExit(
                        f"Matrix dimensions {n_rows} x {n_cols} are incompatible with "
                        f"{n_genes} features and {n_cells} selected cells"
                    )
                if cells_by_genes:
                    output_rows = n_cells
                    output_cols = n_genes
                break
            else:
                raise SystemExit(f"{src} ended before dimensions line")

            for line in inp:
                parts = line.split()
                if len(parts) < 3:
                    continue
                row, col = int(parts[0]), int(parts[1])
                value = parts[2]
                if genes_by_cells:
                    new_col = old_to_new_cell.get(col)
                    if new_col is None:
                        continue
                    body.write(f"{row} {new_col} {value}\n")
                else:
                    new_row = old_to_new_cell.get(row)
                    if new_row is None:
                        continue
                    body.write(f"{new_row} {col} {value}\n")
                nnz += 1

        with gzip.open(dest, "wt", encoding="utf-8") as out:
            out.write("%%MatrixMarket matrix coordinate real general\n")
            out.write(f"{output_rows} {output_cols} {nnz}\n")
        with open(dest, "ab") as out_bin, open(tmp_body, "rb") as body_bin:
            shutil.copyfileobj(body_bin, out_bin)
    finally:
        tmp_body.unlink(missing_ok=True)
    return nnz


def group_output_dir(base: Path, key: tuple[str, ...]) -> Path:
    path = base
    for value in key:
        path = path / safe_label(value)
    return path


def batch_groups(metadata: pd.DataFrame, cell_col: str, split_by: list[str]) -> dict[tuple[str, ...], list[str]]:
    missing = [col for col in split_by if col not in metadata.columns]
    if missing:
        raise SystemExit(f"Metadata is missing split column(s): {', '.join(missing)}")
    cell_key_counts = metadata[[cell_col] + split_by].drop_duplicates().groupby(cell_col).size()
    ambiguous = cell_key_counts[cell_key_counts > 1]
    if not ambiguous.empty:
        raise SystemExit(f"{len(ambiguous)} cells map to multiple --split-by groups")
    grouped: dict[tuple[str, ...], list[str]] = {}
    metadata = metadata.drop_duplicates(cell_col)
    for key_values, group in metadata.groupby(split_by, sort=True, dropna=False):
        key = key_values if isinstance(key_values, tuple) else (key_values,)
        key = tuple(str(value) for value in key)
        cells = group[cell_col].drop_duplicates().astype(str).tolist()
        if cells:
            grouped[key] = cells
    if not grouped:
        raise SystemExit("--split-by produced zero non-empty groups")
    return grouped


def copy_features_to_groups(features_in: Path, group_dirs: dict[tuple[str, ...], Path]) -> None:
    for out_dir in group_dirs.values():
        out_dir.mkdir(parents=True, exist_ok=True)
        copy_text_gz(features_in, out_dir / "features.tsv.gz")


def write_batch_barcodes(
    barcodes_in: Path,
    group_to_cells: dict[tuple[str, ...], set[str]],
    group_dirs: dict[tuple[str, ...], Path],
) -> tuple[dict[int, tuple[tuple[str, ...], int]], dict[tuple[str, ...], int], set[str]]:
    handles = {key: gzip.open(group_dirs[key] / "barcodes.tsv.gz", "wt", encoding="utf-8") for key in group_dirs}
    old_to_group_new: dict[int, tuple[tuple[str, ...], int]] = {}
    counts = {key: 0 for key in group_dirs}
    selected_cells = set().union(*group_to_cells.values()) if group_to_cells else set()
    seen: set[str] = set()
    try:
        with open_text(barcodes_in) as inp:
            for old_index, line in enumerate(inp, start=1):
                cell = line.rstrip("\n").split("\t", 1)[0]
                for key, cells in group_to_cells.items():
                    if cell not in cells:
                        continue
                    counts[key] += 1
                    old_to_group_new[old_index] = (key, counts[key])
                    handles[key].write(line)
                    seen.add(cell)
                    break
    finally:
        for handle in handles.values():
            handle.close()
    return old_to_group_new, counts, selected_cells - seen


def write_batch_cell_totals(
    totals_in: Path,
    group_to_cells: dict[tuple[str, ...], set[str]],
    group_dirs: dict[tuple[str, ...], Path],
) -> None:
    cell_to_key = {cell: key for key, cells in group_to_cells.items() for cell in cells}
    handles = {key: gzip.open(group_dirs[key] / "cell_total_counts.tsv.gz", "wt", encoding="utf-8") for key in group_dirs}
    try:
        with open_text(totals_in) as inp:
            header = inp.readline()
            if header:
                for handle in handles.values():
                    handle.write(header)
            for line in inp:
                cell = line.rstrip("\n").split("\t", 1)[0]
                key = cell_to_key.get(cell)
                if key is not None:
                    handles[key].write(line)
    finally:
        for handle in handles.values():
            handle.close()


def write_batch_matrices(
    src: Path,
    group_dirs: dict[tuple[str, ...], Path],
    n_genes: int,
    group_cell_counts: dict[tuple[str, ...], int],
    old_to_group_new: dict[int, tuple[tuple[str, ...], int]],
) -> dict[tuple[str, ...], int]:
    tmp_paths: dict[tuple[str, ...], Path] = {}
    body_handles = {}
    nnz = {key: 0 for key in group_dirs}
    output_dims: dict[tuple[str, ...], tuple[int, int]] = {}
    try:
        for key in group_dirs:
            tmp = tempfile.NamedTemporaryFile(prefix="subset_10x_batch_body_", suffix=".gz", delete=False)
            tmp_paths[key] = Path(tmp.name)
            tmp.close()
            body_handles[key] = gzip.open(tmp_paths[key], "wt", encoding="utf-8")
        with open_text(src) as inp:
            banner = inp.readline()
            if not banner.startswith("%%MatrixMarket"):
                raise SystemExit(f"{src} is not a Matrix Market file")
            for line in inp:
                if line.startswith("%"):
                    continue
                n_rows, n_cols, _ = [int(x) for x in line.split()[:3]]
                genes_by_cells = n_rows == n_genes
                cells_by_genes = n_cols == n_genes
                if not genes_by_cells and not cells_by_genes:
                    raise SystemExit(f"Matrix dimensions {n_rows} x {n_cols} are incompatible with {n_genes} features")
                for key, n_cells in group_cell_counts.items():
                    output_dims[key] = (n_genes, n_cells) if genes_by_cells else (n_cells, n_genes)
                break
            else:
                raise SystemExit(f"{src} ended before dimensions line")

            for line in inp:
                parts = line.split()
                if len(parts) < 3:
                    continue
                row, col = int(parts[0]), int(parts[1])
                value = parts[2]
                old_cell = col if genes_by_cells else row
                mapped = old_to_group_new.get(old_cell)
                if mapped is None:
                    continue
                key, new_cell = mapped
                if genes_by_cells:
                    body_handles[key].write(f"{row} {new_cell} {value}\n")
                else:
                    body_handles[key].write(f"{new_cell} {col} {value}\n")
                nnz[key] += 1
    finally:
        for handle in body_handles.values():
            handle.close()

    for key, out_dir in group_dirs.items():
        rows, cols = output_dims[key]
        dest = out_dir / "matrix.mtx.gz"
        with gzip.open(dest, "wt", encoding="utf-8") as out:
            out.write("%%MatrixMarket matrix coordinate real general\n")
            out.write(f"{rows} {cols} {nnz[key]}\n")
        with open(dest, "ab") as out_bin, open(tmp_paths[key], "rb") as body_bin:
            shutil.copyfileobj(body_bin, out_bin)
        tmp_paths[key].unlink(missing_ok=True)
    return nnz


def write_batch_split(
    matrix_in: Path,
    features_in: Path,
    barcodes_in: Path,
    totals_in: Path | None,
    metadata: pd.DataFrame,
    cell_col: str,
    split_by: list[str],
    out_dir: Path,
    filters: list[str],
) -> None:
    group_cells = batch_groups(metadata, cell_col, split_by)
    group_to_cells = {key: set(cells) for key, cells in group_cells.items()}
    group_dirs = {key: group_output_dir(out_dir, key) for key in group_cells}
    out_dir.mkdir(parents=True, exist_ok=True)
    copy_features_to_groups(features_in, group_dirs)
    old_to_group_new, group_cell_counts, missing = write_batch_barcodes(barcodes_in, group_to_cells, group_dirs)
    if missing:
        raise SystemExit(f"Input 10x barcodes are missing {len(missing)} selected cells")
    n_genes = len(read_one_column(features_in))
    nnz = write_batch_matrices(matrix_in, group_dirs, n_genes, group_cell_counts, old_to_group_new)
    if totals_in is not None:
        write_batch_cell_totals(totals_in, group_to_cells, group_dirs)
    rows = []
    for key in sorted(group_dirs):
        row = {
            "group_label": "/".join(safe_label(value) for value in key),
            "out_dir": str(group_dirs[key]),
            "n_genes": n_genes,
            "n_cells": group_cell_counts[key],
            "nnz": nnz[key],
            "metadata_filters": ";".join(filters),
        }
        for col, value in zip(split_by, key):
            row[col] = value
        rows.append(row)
        pd.DataFrame(
            [
                {"metric": "n_genes", "value": n_genes},
                {"metric": "n_cells", "value": group_cell_counts[key]},
                {"metric": "nnz", "value": nnz[key]},
                {"metric": "metadata_filters", "value": ";".join(filters)},
                {"metric": "split_by", "value": ",".join(split_by)},
                {"metric": "group_label", "value": row["group_label"]},
            ]
        ).to_csv(group_dirs[key] / "subset_summary.tsv", sep="\t", index=False)
    pd.DataFrame(rows).to_csv(out_dir / "split_summary.tsv", sep="\t", index=False)
    print(f"Wrote {len(rows)} split 10x directories to {out_dir}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input-10x-dir", type=Path, required=True)
    ap.add_argument("--metadata", type=Path, required=True)
    ap.add_argument("--metadata-filter", action="append", default=[])
    ap.add_argument("--split-by", default="", help="Comma-separated metadata columns for amortized batch splitting")
    ap.add_argument("--cell-id-col", default="cell_id")
    ap.add_argument("--out-dir", type=Path, required=True)
    args = ap.parse_args()

    matrix_in = resolve_10x_path(args.input_10x_dir, ["matrix.mtx.gz", "matrix.mtx"])
    features_in = resolve_10x_path(args.input_10x_dir, ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"])
    barcodes_in = resolve_10x_path(args.input_10x_dir, ["barcodes.tsv.gz", "barcodes.tsv"])
    totals_in = None
    for name in ["cell_total_counts.tsv.gz", "cell_total_counts.tsv"]:
        path = args.input_10x_dir / name
        if path.exists():
            totals_in = path
            break

    split_by = [col.strip() for col in args.split_by.split(",") if col.strip()]
    if split_by:
        metadata = read_filtered_metadata(args.metadata, args.cell_id_col, args.metadata_filter)
        write_batch_split(
            matrix_in,
            features_in,
            barcodes_in,
            totals_in,
            metadata,
            args.cell_id_col,
            split_by,
            args.out_dir,
            args.metadata_filter,
        )
        return

    selected_cells = read_selected_cells(args.metadata, args.cell_id_col, args.metadata_filter)
    selected_set = set(selected_cells)
    all_cells = read_one_column(barcodes_in)
    selected_positions = {i for i, cell in enumerate(all_cells, start=1) if cell in selected_set}
    missing = selected_set - set(all_cells)
    if missing:
        raise SystemExit(f"Input 10x barcodes are missing {len(missing)} selected cells")

    args.out_dir.mkdir(parents=True, exist_ok=True)
    features_out = args.out_dir / "features.tsv.gz"
    barcodes_out = args.out_dir / "barcodes.tsv.gz"
    matrix_out = args.out_dir / "matrix.mtx.gz"
    copy_text_gz(features_in, features_out)
    old_to_new = subset_barcodes(barcodes_in, barcodes_out, selected_positions)
    n_genes = len(read_one_column(features_in))
    nnz = write_subset_matrix(matrix_in, matrix_out, n_genes, len(old_to_new), old_to_new)
    if totals_in is not None:
        subset_cell_totals(totals_in, args.out_dir / "cell_total_counts.tsv.gz", selected_set)
    summary = pd.DataFrame(
        [
            {"metric": "n_genes", "value": n_genes},
            {"metric": "n_cells", "value": len(old_to_new)},
            {"metric": "nnz", "value": nnz},
            {"metric": "metadata_filters", "value": ";".join(args.metadata_filter)},
        ]
    )
    summary.to_csv(args.out_dir / "subset_summary.tsv", sep="\t", index=False)
    print(f"Wrote {len(old_to_new)} cells and {nnz} nonzero entries to {args.out_dir}")


if __name__ == "__main__":
    main()
