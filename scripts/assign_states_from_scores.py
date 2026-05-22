#!/usr/bin/env python3
"""Assign cells to states from a cell x state score table."""

from __future__ import annotations

import argparse
import json

import pandas as pd


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series
    return series.astype(str).str.strip().str.lower().isin({"true", "t", "1", "yes", "y"})


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", required=True, help="Long score TSV/TSV.GZ with cell_id, state, score")
    parser.add_argument("--metadata", required=True, help="Cell metadata TSV/TSV.GZ")
    parser.add_argument("--out", required=True, help="Output membership TSV/TSV.GZ")
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--cell-type-col", default="")
    parser.add_argument("--state-cell-type-map", default="", help="Optional TSV with state and cell_type columns")
    parser.add_argument("--method", choices=["quantile", "top"], default="quantile")
    parser.add_argument("--quantile", type=float, default=0.75)
    parser.add_argument("--within", choices=["all", "cell_type", "donor", "cell_type_donor"], default="cell_type")
    parser.add_argument("--donor-col", default="donor_id")
    parser.add_argument("--allow-cross-cell-type", action="store_true")
    args = parser.parse_args()

    metadata = read_table(args.metadata)
    if args.cell_id_col not in metadata.columns:
        raise SystemExit(f"Metadata is missing cell ID column: {args.cell_id_col}")
    scores = read_table(args.scores)
    for column in ["cell_id", "state", "score"]:
        if column not in scores.columns:
            raise SystemExit(f"Scores are missing required column: {column}")
    scores["score"] = pd.to_numeric(scores["score"], errors="coerce")

    meta_cols = [args.cell_id_col]
    if args.cell_type_col:
        if args.cell_type_col not in metadata.columns:
            raise SystemExit(f"Metadata is missing cell type column: {args.cell_type_col}")
        meta_cols.append(args.cell_type_col)
    if args.within in {"donor", "cell_type_donor"}:
        if args.donor_col not in metadata.columns:
            raise SystemExit(f"Metadata is missing donor column: {args.donor_col}")
        meta_cols.append(args.donor_col)
    frame = scores.merge(metadata[meta_cols].rename(columns={args.cell_id_col: "cell_id"}), on="cell_id", how="inner")

    if args.state_cell_type_map:
        state_map = read_table(args.state_cell_type_map)
        for column in ["state", "cell_type"]:
            if column not in state_map.columns:
                raise SystemExit(f"State cell type map is missing required column: {column}")
        frame = frame.merge(state_map[["state", "cell_type"]], on="state", how="left")
        if args.cell_type_col and not args.allow_cross_cell_type:
            frame = frame.loc[frame[args.cell_type_col].astype(str) == frame["cell_type"].astype(str)].copy()
    elif args.cell_type_col and not args.allow_cross_cell_type:
        print("Warning: no --state-cell-type-map supplied; state scores are not filtered by cell type")

    if frame.empty:
        raise SystemExit("No score rows remained after metadata and cell-type filtering")

    if args.method == "top":
        frame["threshold"] = frame.groupby("cell_id")["score"].transform("max")
        frame["in_state"] = frame["score"] == frame["threshold"]
    else:
        group_cols = ["state"]
        if args.within in {"cell_type", "cell_type_donor"}:
            if not args.cell_type_col:
                raise SystemExit("--cell-type-col is required when --within uses cell_type")
            group_cols.append(args.cell_type_col)
        if args.within in {"donor", "cell_type_donor"}:
            group_cols.append(args.donor_col)
        frame["threshold"] = frame.groupby(group_cols)["score"].transform(lambda s: s.quantile(args.quantile))
        frame["in_state"] = frame["score"] >= frame["threshold"]

    frame["rule"] = json.dumps(
        {
            "method": args.method,
            "quantile": args.quantile,
            "within": args.within,
            "allow_cross_cell_type": args.allow_cross_cell_type,
        },
        sort_keys=True,
        separators=(",", ":"),
    )
    out_cols = ["cell_id", "state", "in_state", "score", "threshold", "rule"]
    if args.cell_type_col and args.cell_type_col in frame.columns:
        frame["cell_type"] = frame[args.cell_type_col]
        out_cols.append("cell_type")
    frame[out_cols].to_csv(args.out, sep="\t", index=False, compression="infer")


if __name__ == "__main__":
    main()
