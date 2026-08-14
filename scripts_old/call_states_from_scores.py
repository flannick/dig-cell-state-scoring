#!/usr/bin/env python3
"""Call multi-label cell states from continuous scores and calibrated thresholds."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


ACTIVE = "active"
INACTIVE = "inactive"
AMBIGUOUS = "ambiguous"
INSUFFICIENT = "insufficient_coverage"
QC_FLAGGED = "qc_flagged"


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def write_table(frame: pd.DataFrame, path: str) -> None:
    frame.to_csv(path, sep="\t", index=False, compression="infer")


def load_rules(path: str) -> dict[str, dict[str, Any]]:
    if not path:
        return {}
    with open(path, "r", encoding="utf-8") as handle:
        obj = json.load(handle)
    states = obj.get("states", obj if isinstance(obj, list) else [])
    rules: dict[str, dict[str, Any]] = {}
    for state in states:
        if "state" in state:
            name = state["state"]
        else:
            name = state["name"]
        rules[name] = state
    return rules


def default_kind(state: str) -> str:
    if state.startswith("qc_"):
        return "qc_flag"
    if any(token in state for token in ["stress", "interferon", "proliferative", "cell_cycle", "senescence"]):
        return "process_flag" if state.startswith("qc_") else "biological"
    return "biological"


def score_column(scores: pd.DataFrame) -> str:
    if "ucell_score" in scores.columns and scores["ucell_score"].notna().any():
        return "ucell_score"
    return "score"


def normalize_thresholds(thresholds: pd.DataFrame) -> pd.DataFrame:
    required = {"state", "threshold_value"}
    missing = sorted(required - set(thresholds.columns))
    if missing:
        raise SystemExit(f"Thresholds are missing required column(s): {', '.join(missing)}")
    cols = ["state", "threshold_value"]
    optional = ["threshold_method", "null_percentile", "mixture_boundary", "n_cells_used"]
    cols.extend([c for c in optional if c in thresholds.columns])
    out = thresholds[cols].drop_duplicates("state").copy()
    out["threshold_value"] = pd.to_numeric(out["threshold_value"], errors="coerce")
    return out


def add_rule_columns(frame: pd.DataFrame, rules: dict[str, dict[str, Any]]) -> pd.DataFrame:
    kinds = []
    required_active = []
    required_inactive = []
    excluded_active = []
    for state in frame["state"]:
        rule = rules.get(state, {})
        kinds.append(rule.get("kind", default_kind(state)))
        required_active.append(";".join(rule.get("required_active", [])))
        required_inactive.append(";".join(rule.get("required_inactive", [])))
        excluded_active.append(";".join(rule.get("excluded_active", [])))
    frame["state_kind"] = kinds
    frame["required_active"] = required_active
    frame["required_inactive"] = required_inactive
    frame["excluded_active"] = excluded_active
    return frame


def initial_calls(frame: pd.DataFrame, min_markers_present: int, min_coverage_fraction: float) -> pd.DataFrame:
    coverage_col = "marker_coverage_fraction" if "marker_coverage_fraction" in frame.columns else None
    present_col = "marker_genes_present" if "marker_genes_present" in frame.columns else "n_markers_found"
    if present_col not in frame.columns:
        frame[present_col] = np.nan
    frame[present_col] = pd.to_numeric(frame[present_col], errors="coerce")
    if coverage_col is None:
        frame["marker_coverage_fraction"] = np.nan
        coverage_col = "marker_coverage_fraction"
    frame[coverage_col] = pd.to_numeric(frame[coverage_col], errors="coerce")

    insufficient = (frame[present_col] < min_markers_present) | (frame[coverage_col] < min_coverage_fraction)
    missing_threshold = frame["threshold_value"].isna()
    missing_score = frame["score_for_call"].isna()
    active = frame["score_for_call"] >= frame["threshold_value"]
    frame["confidence"] = frame["score_for_call"] - frame["threshold_value"]
    frame["call"] = np.where(active, ACTIVE, INACTIVE)
    frame.loc[missing_threshold | missing_score, "call"] = AMBIGUOUS
    frame.loc[insufficient, "call"] = INSUFFICIENT
    frame["reason"] = np.where(active, "score_above_null99", "score_below_threshold")
    frame.loc[missing_threshold, "reason"] = "missing_threshold"
    frame.loc[missing_score, "reason"] = "missing_score"
    frame.loc[insufficient, "reason"] = "insufficient_marker_coverage"
    return frame


def apply_rules(frame: pd.DataFrame) -> pd.DataFrame:
    active_pairs = frame.loc[frame["call"].eq(ACTIVE), ["cell_id", "state"]]
    active_by_cell = active_pairs.groupby("cell_id")["state"].agg(set).to_dict()

    calls = frame["call"].tolist()
    reasons = frame["reason"].tolist()
    for idx, row in frame.iterrows():
        if row["call"] != ACTIVE:
            continue
        active_states = active_by_cell.get(row["cell_id"], set())
        missing_required = [s for s in str(row["required_active"]).split(";") if s and s not in active_states]
        forbidden = [s for s in str(row["excluded_active"]).split(";") if s and s in active_states]
        should_be_inactive = [s for s in str(row["required_inactive"]).split(";") if s and s in active_states]
        if missing_required:
            calls[idx] = AMBIGUOUS
            reasons[idx] = "required_state_not_active:" + ",".join(missing_required)
        elif should_be_inactive:
            calls[idx] = AMBIGUOUS
            reasons[idx] = "required_inactive_state_active:" + ",".join(should_be_inactive)
        elif forbidden:
            calls[idx] = QC_FLAGGED
            reasons[idx] = "excluded_qc_flag_active:" + ",".join(forbidden)
    frame["call"] = calls
    frame["reason"] = reasons
    return frame


def build_annotations(calls: pd.DataFrame, metadata: pd.DataFrame | None, cell_id_col: str, parent_cell_type_col: str) -> pd.DataFrame:
    active = calls.loc[calls["call"].isin([ACTIVE, QC_FLAGGED])].copy()
    active_bio = active.loc[(active["call"].eq(ACTIVE)) & (active["state_kind"].eq("biological"))]
    process = active.loc[(active["call"].eq(ACTIVE)) & (active["state_kind"].eq("process_flag"))]
    qc = active.loc[active["state_kind"].eq("qc_flag")]

    cells = pd.DataFrame({"cell_id": sorted(calls["cell_id"].unique())})
    for name, subset in [
        ("active_biological_states", active_bio),
        ("active_process_flags", process),
        ("qc_flags", qc),
    ]:
        grouped = subset.groupby("cell_id")["state"].agg(lambda x: ";".join(sorted(set(x)))).rename(name)
        cells = cells.merge(grouped, how="left", left_on="cell_id", right_index=True)
        cells[name] = cells[name].fillna("none")

    if metadata is not None and parent_cell_type_col and parent_cell_type_col in metadata.columns:
        meta = metadata[[cell_id_col, parent_cell_type_col]].rename(
            columns={cell_id_col: "cell_id", parent_cell_type_col: "parent_cell_type"}
        )
        cells = cells.merge(meta.drop_duplicates("cell_id"), on="cell_id", how="left")
    else:
        cells["parent_cell_type"] = ""

    def interpretation(row: pd.Series) -> str:
        states = row["active_biological_states"]
        flags = row["active_process_flags"]
        qc_flags = row["qc_flags"]
        parts = []
        if states != "none":
            parts.append(states.replace(";", "; "))
        if flags != "none":
            parts.append("process flags: " + flags.replace(";", "; "))
        if qc_flags != "none":
            parts.append("QC flags: " + qc_flags.replace(";", "; "))
        return " | ".join(parts) if parts else "no active state calls"

    cells["primary_interpretation"] = cells.apply(interpretation, axis=1)
    return cells[
        [
            "cell_id",
            "parent_cell_type",
            "active_biological_states",
            "active_process_flags",
            "qc_flags",
            "primary_interpretation",
        ]
    ]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", required=True, help="Long score TSV/TSV.GZ")
    parser.add_argument("--thresholds", required=True, help="Threshold TSV/TSV.GZ")
    parser.add_argument("--out", required=True, help="Long state-call TSV/TSV.GZ")
    parser.add_argument("--annotation-out", default="", help="Optional multi-label cell annotation TSV/TSV.GZ")
    parser.add_argument("--metadata", default="", help="Optional metadata TSV/TSV.GZ for parent cell type")
    parser.add_argument("--cell-id-col", default="cell_id")
    parser.add_argument("--parent-cell-type-col", default="")
    parser.add_argument("--rules", default="", help="Optional JSON rule config")
    parser.add_argument("--min-markers-present", type=int, default=5)
    parser.add_argument("--min-coverage-fraction", type=float, default=0.5)
    args = parser.parse_args()

    scores = read_table(args.scores)
    for column in ["cell_id", "state"]:
        if column not in scores.columns:
            raise SystemExit(f"Scores are missing required column: {column}")
    call_score = score_column(scores)
    scores["score_for_call"] = pd.to_numeric(scores[call_score], errors="coerce")
    thresholds = normalize_thresholds(read_table(args.thresholds))
    rules = load_rules(args.rules)

    frame = scores.merge(thresholds, on="state", how="left")
    frame = add_rule_columns(frame, rules)
    frame = initial_calls(frame, args.min_markers_present, args.min_coverage_fraction)
    frame = apply_rules(frame)
    frame["in_state"] = frame["call"].eq(ACTIVE)
    frame["call_score_column"] = call_score

    preferred = [
        "cell_id",
        "state",
        "state_kind",
        "score_for_call",
        "threshold_value",
        "call",
        "confidence",
        "reason",
        "in_state",
        "call_score_column",
    ]
    remaining = [c for c in frame.columns if c not in preferred]
    write_table(frame[preferred + remaining], args.out)

    if args.annotation_out:
        metadata = read_table(args.metadata) if args.metadata else None
        annotations = build_annotations(frame, metadata, args.cell_id_col, args.parent_cell_type_col)
        write_table(annotations, args.annotation_out)


if __name__ == "__main__":
    main()
