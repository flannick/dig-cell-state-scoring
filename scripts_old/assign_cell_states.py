#!/usr/bin/env python3
"""Assign single cells to marker-defined or metadata-defined states."""

from __future__ import annotations

import argparse
import json

import pandas as pd


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def write_table(df: pd.DataFrame, path: str) -> None:
    df.to_csv(path, sep="\t", index=False, compression="infer")


def load_spec(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def score_module(expr: pd.DataFrame, genes: list[str], cell_ids: pd.Series) -> pd.DataFrame:
    grid = pd.MultiIndex.from_product([cell_ids, genes], names=["cell_id", "gene"]).to_frame(index=False)
    sub = expr.loc[expr["gene"].isin(genes), ["cell_id", "gene", "expression"]].copy()
    sub["expression"] = pd.to_numeric(sub["expression"], errors="coerce")
    merged = grid.merge(sub, on=["cell_id", "gene"], how="left")
    merged["expression"] = merged["expression"].fillna(0.0)
    return merged.groupby("cell_id", as_index=False)["expression"].mean().rename(columns={"expression": "score"})


def assign_module_quantile(
    metadata: pd.DataFrame,
    expr: pd.DataFrame,
    state: dict,
    cell_col: str,
    donor_col: str,
) -> pd.DataFrame:
    genes = state["genes"]
    direction = state.get("direction", "high")
    quantile = float(state.get("quantile", 0.75))
    within = state.get("within", "all")

    base_cols = [cell_col]
    if within == "donor":
        base_cols.append(donor_col)
    scores = score_module(expr, genes, metadata[cell_col])
    frame = metadata[base_cols].rename(columns={cell_col: "cell_id"}).merge(scores, on="cell_id", how="left")
    if within == "donor":
        frame["threshold"] = frame.groupby(donor_col)["score"].transform(lambda s: s.quantile(quantile))
    elif within == "all":
        frame["threshold"] = frame["score"].quantile(quantile)
    else:
        raise ValueError(f"Unsupported module_quantile within value: {within}")

    if direction == "high":
        frame["in_state"] = frame["score"] >= frame["threshold"]
    elif direction == "low":
        frame["in_state"] = frame["score"] <= frame["threshold"]
    else:
        raise ValueError(f"Unsupported module_quantile direction: {direction}")

    frame["state"] = state["name"]
    frame["rule"] = json.dumps(state, sort_keys=True, separators=(",", ":"))
    return frame[["cell_id", "state", "in_state", "score", "threshold", "rule"]]


def assign_metadata_rule(metadata: pd.DataFrame, state: dict, cell_col: str) -> pd.DataFrame:
    column = state["column"]
    if column not in metadata.columns:
        raise ValueError(f"Metadata column not found for state {state['name']}: {column}")

    rule_type = state["type"]
    values = state.get("values")
    if rule_type == "metadata_equals":
        if values is not None:
            raise ValueError("metadata_equals expects 'value', not 'values'")
        in_state = metadata[column].astype(str) == str(state["value"])
    elif rule_type == "metadata_in":
        if not isinstance(values, list):
            raise ValueError("metadata_in expects a list-valued 'values' field")
        in_state = metadata[column].astype(str).isin([str(v) for v in values])
    else:
        raise ValueError(f"Unsupported metadata rule type: {rule_type}")

    out = pd.DataFrame(
        {
            "cell_id": metadata[cell_col],
            "state": state["name"],
            "in_state": in_state,
            "score": pd.NA,
            "threshold": pd.NA,
            "rule": json.dumps(state, sort_keys=True, separators=(",", ":")),
        }
    )
    return out


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True, help="Cell metadata TSV/TSV.GZ")
    parser.add_argument("--expression", required=True, help="Long expression TSV/TSV.GZ")
    parser.add_argument("--state-spec", required=True, help="JSON state specification")
    parser.add_argument("--out", required=True, help="Output membership TSV/TSV.GZ")
    args = parser.parse_args()

    spec = load_spec(args.state_spec)
    cell_col = spec.get("cell_id_col", "cell_id")
    donor_col = spec.get("donor_col", "donor_id")

    metadata = read_table(args.metadata)
    if cell_col not in metadata.columns:
        raise SystemExit(f"Metadata is missing required cell id column: {cell_col}")

    expr = read_table(args.expression)
    required_expr = {"cell_id", "gene", "expression"}
    missing_expr = required_expr - set(expr.columns)
    if missing_expr:
        raise SystemExit(f"Expression table is missing required columns: {sorted(missing_expr)}")

    outputs = []
    for state in spec["states"]:
        state_type = state["type"]
        if state_type == "module_quantile":
            if state.get("within", "all") == "donor" and donor_col not in metadata.columns:
                raise SystemExit(
                    f"Metadata is missing donor column '{donor_col}' required by donor-threshold state: {state['name']}"
                )
            outputs.append(assign_module_quantile(metadata, expr, state, cell_col, donor_col))
        elif state_type in {"metadata_equals", "metadata_in"}:
            outputs.append(assign_metadata_rule(metadata, state, cell_col))
        else:
            raise SystemExit(f"Unsupported state type for {state.get('name', '<unnamed>')}: {state_type}")

    out = pd.concat(outputs, ignore_index=True)
    out["in_state"] = out["in_state"].fillna(False).astype(bool)
    write_table(out, args.out)


if __name__ == "__main__":
    main()
