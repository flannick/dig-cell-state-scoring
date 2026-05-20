#!/usr/bin/env python3
"""Assign genes to states from GMT markers and optional state association DE."""

from __future__ import annotations

import argparse

import pandas as pd


def read_gmt(path: str) -> pd.DataFrame:
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            state, library, *genes = parts
            for gene in genes:
                if gene:
                    rows.append({"gene": gene, "state": state, "library": library, "is_marker": True})
    return pd.DataFrame(rows).drop_duplicates()


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gmt", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--state-association-de", default="")
    parser.add_argument("--state-regex", default="")
    parser.add_argument("--qvalue-threshold", type=float, default=0.05)
    parser.add_argument("--min-log-fc", type=float, default=0.0)
    args = parser.parse_args()

    markers = read_gmt(args.gmt)
    if args.state_regex:
        markers = markers.loc[markers["state"].str.contains(args.state_regex, regex=True)].copy()
    markers["assignment_source"] = "gmt_marker"
    for col in ["log_fc", "pvalue", "qvalue"]:
        markers[col] = pd.NA
    markers["assignment_pass"] = True

    outputs = [markers[["gene", "state", "library", "assignment_source", "is_marker", "log_fc", "pvalue", "qvalue", "assignment_pass"]]]
    if args.state_association_de:
        de = read_table(args.state_association_de)
        rename = {"log2_fc": "log_fc", "logFC": "log_fc", "FDR": "qvalue", "PValue": "pvalue"}
        de = de.rename(columns={k: v for k, v in rename.items() if k in de.columns})
        required = {"gene", "state", "log_fc", "pvalue", "qvalue"}
        missing = required - set(de.columns)
        if missing:
            raise SystemExit(f"State association DE table is missing required columns: {sorted(missing)}")
        de = de[["gene", "state", "log_fc", "pvalue", "qvalue"]].copy()
        if args.state_regex:
            de = de.loc[de["state"].str.contains(args.state_regex, regex=True)].copy()
        marker_pairs = set(zip(markers["gene"], markers["state"]))
        de["library"] = pd.NA
        de["assignment_source"] = "state_association_de"
        de["is_marker"] = [pair in marker_pairs for pair in zip(de["gene"], de["state"])]
        de["assignment_pass"] = (pd.to_numeric(de["qvalue"], errors="coerce") <= args.qvalue_threshold) & (
            pd.to_numeric(de["log_fc"], errors="coerce") > args.min_log_fc
        )
        outputs.append(de[["gene", "state", "library", "assignment_source", "is_marker", "log_fc", "pvalue", "qvalue", "assignment_pass"]])

    out = pd.concat(outputs, ignore_index=True).drop_duplicates()
    out.to_csv(args.out, sep="\t", index=False, compression="infer")


if __name__ == "__main__":
    main()
