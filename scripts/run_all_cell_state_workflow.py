#!/usr/bin/env python3
"""Config-driven multi-group cell-state workflow runner."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml


SCRIPT_DIR = Path(__file__).resolve().parent


def safe_label(value: object) -> str:
    value = str(value).strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return re.sub(r"_+", "_", value).strip("_") or "unknown"


def read_config(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise SystemExit("Workflow config must be a YAML mapping")
    return data


def require(config: dict, key: str) -> str:
    value = config.get(key)
    if not value:
        raise SystemExit(f"Config is missing required key: {key}")
    return str(value)


def run_command(command: list[str], dry_run: bool, log_path: Path) -> tuple[str, float]:
    start = time.time()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(" ".join(command) + "\n", encoding="utf-8")
    if dry_run:
        return "planned", 0.0
    result = subprocess.run(command, text=True, capture_output=True)
    (log_path.with_suffix(".stdout.txt")).write_text(result.stdout, encoding="utf-8")
    (log_path.with_suffix(".stderr.txt")).write_text(result.stderr, encoding="utf-8")
    elapsed = time.time() - start
    if result.returncode != 0:
        return f"failed:{result.returncode}", elapsed
    return "completed", elapsed


def metadata_filter(split_by: list[str], group: pd.Series) -> str:
    return ";".join(f"{col}={group[col]}" for col in split_by)


def parent_filter(split_by: list[str], group: pd.Series) -> str:
    if len(split_by) == 1:
        return f"{split_by[0]}={group[split_by[0]]}"
    return ";".join(f"{col}={group[col]}" for col in split_by)


def split_dir(base: Path, split_by: list[str], group: pd.Series) -> Path:
    path = base
    for col in split_by:
        path = path / safe_label(group[col])
    return path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    config = read_config(args.config)
    out_dir = Path(require(config, "out_dir"))
    out_dir.mkdir(parents=True, exist_ok=True)
    dry_run = bool(args.dry_run or config.get("dry_run", False))
    split_by = config.get("split_by", ["tissue", "cell_type"])
    if isinstance(split_by, str):
        split_by = [x.strip() for x in split_by.split(",") if x.strip()]
    split_by = list(split_by)
    metadata_path = Path(require(config, "metadata"))
    metadata = pd.read_csv(metadata_path, sep="\t", compression="infer", low_memory=False, dtype=str).fillna("")
    missing_split = [col for col in split_by if col not in metadata.columns]
    if missing_split:
        raise SystemExit(f"Metadata is missing split_by column(s): {', '.join(missing_split)}")
    groups = metadata[split_by].drop_duplicates().sort_values(split_by).reset_index(drop=True)

    commands = []
    timings = []
    failures = []
    paths = {
        "rank_subsets": out_dir / "rank_10x_by_group",
        "raw_subsets": out_dir / "raw_10x_by_group",
        "split_gmts": out_dir / "state_gmts_by_group",
        "logs": out_dir / "command_logs",
        "workflow": out_dir / "groups",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)

    split_cmd = [
        sys.executable,
        str(SCRIPT_DIR / "subset_10x_by_metadata.py"),
        "--input-10x-dir",
        require(config, "rank_10x_dir"),
        "--metadata",
        str(metadata_path),
        "--split-by",
        ",".join(split_by),
        "--out-dir",
        str(paths["rank_subsets"]),
    ]
    if config.get("metadata_filter"):
        split_cmd.extend(["--metadata-filter", str(config["metadata_filter"])])
    status, elapsed = run_command(split_cmd, dry_run, paths["logs"] / "split_rank_10x.sh")
    timings.append({"step": "split_rank_10x", "status": status, "seconds": round(elapsed, 3)})
    commands.append({"step": "split_rank_10x", "command": " ".join(split_cmd)})
    if status.startswith("failed"):
        failures.append({"step": "split_rank_10x", "status": status})

    raw_10x_dir = str(config.get("raw_10x_dir", config.get("rank_10x_dir", "")))
    if raw_10x_dir and raw_10x_dir != str(config.get("rank_10x_dir")):
        raw_cmd = split_cmd.copy()
        raw_cmd[raw_cmd.index(require(config, "rank_10x_dir"))] = raw_10x_dir
        raw_cmd[raw_cmd.index(str(paths["rank_subsets"]))] = str(paths["raw_subsets"])
        status, elapsed = run_command(raw_cmd, dry_run, paths["logs"] / "split_raw_10x.sh")
        timings.append({"step": "split_raw_10x", "status": status, "seconds": round(elapsed, 3)})
        commands.append({"step": "split_raw_10x", "command": " ".join(raw_cmd)})
    else:
        paths["raw_subsets"] = paths["rank_subsets"]

    gmt_cmd = [
        sys.executable,
        str(SCRIPT_DIR / "split_gmt_by_manifest.py"),
        "--gmt",
        require(config, "states_gmt"),
        "--state-manifest",
        require(config, "state_manifest"),
        "--out-dir",
        str(paths["split_gmts"]),
        "--require-all-states",
    ]
    status, elapsed = run_command(gmt_cmd, dry_run, paths["logs"] / "split_state_gmt.sh")
    timings.append({"step": "split_state_gmt", "status": status, "seconds": round(elapsed, 3)})
    commands.append({"step": "split_state_gmt", "command": " ".join(gmt_cmd)})

    for _, group in groups.iterrows():
        label = "__".join(safe_label(group[col]) for col in split_by)
        group_dir = paths["workflow"] / label
        rank_dir = split_dir(paths["rank_subsets"], split_by, group)
        raw_dir = split_dir(paths["raw_subsets"], split_by, group)
        gmt_path = paths["split_gmts"] / safe_label(group[split_by[0]]) / f"{safe_label(group[split_by[-1]])}.gmt"
        if len(split_by) == 1:
            gmt_path = paths["split_gmts"] / "all" / f"{safe_label(group[split_by[0]])}.gmt"
        scoring_dir = group_dir / "workflow"
        expr_dir = group_dir / "expression"
        gmt_out = group_dir / "state_expression_gmts"
        de_dir = group_dir / "de"
        command_log = paths["logs"] / f"{label}.commands.txt"
        group_commands = []

        scoring_cmd = [
            sys.executable,
            str(SCRIPT_DIR / "run_cmdkp_state_scoring.py"),
            "--expression-matrix",
            require(config, "expression_matrix"),
            "--cell-metadata",
            str(metadata_path),
            "--states-gmt",
            str(gmt_path),
            "--state-manifest",
            require(config, "state_manifest"),
            "--require-state-manifest",
            "--rank-10x-dir",
            str(rank_dir),
            "--qc-raw-10x-dir",
            str(raw_dir),
            "--parent-cell-filter",
            parent_filter(split_by, group),
            "--out-dir",
            str(scoring_dir),
        ]
        if config.get("qc_gmt"):
            scoring_cmd.extend(["--qc-states-gmt", str(config["qc_gmt"])])
        if config.get("allow_small_rank_universe"):
            scoring_cmd.append("--allow-small-rank-universe")
        group_commands.append(scoring_cmd)

        expr_cmd = [
            sys.executable,
            str(SCRIPT_DIR / "summarize_state_expression.py"),
            "--raw-10x-dir",
            str(raw_dir),
            "--metadata",
            str(metadata_path),
            "--cell-state-activity",
            str(scoring_dir / "cell_state_activity.tsv.gz"),
            "--states-gmt",
            str(gmt_path),
            "--out-dir",
            str(expr_dir),
        ]
        if config.get("output_format"):
            expr_cmd.extend(["--output-format", str(config["output_format"])])
        group_commands.append(expr_cmd)

        make_gmt_cmd = [
            sys.executable,
            str(SCRIPT_DIR / "make_state_expression_gmts.py"),
            "--state-expression-specificity",
            str(expr_dir / "all_gene_state_expression_specificity_cp10k.tsv.gz"),
            "--original-state-gmt",
            str(gmt_path),
            "--out-dir",
            str(gmt_out),
        ]
        group_commands.append(make_gmt_cmd)

        if config.get("phenotypes"):
            de_cmd = [
                sys.executable,
                str(SCRIPT_DIR / "run_state_phenotype_regression.py"),
                "--raw-10x-dir",
                str(raw_dir),
                "--metadata",
                str(metadata_path),
                "--cell-state-activity",
                str(scoring_dir / "cell_state_activity.tsv.gz"),
                "--donor-state-expression",
                str(expr_dir / "donor_state_weighted_expression_cp10k.tsv.gz"),
                "--phenotypes",
                str(config["phenotypes"]),
                "--out-dir",
                str(de_dir),
            ]
            if config.get("phenotype_config"):
                de_cmd.extend(["--phenotype-config", str(config["phenotype_config"])])
            group_commands.append(de_cmd)

        for idx, command in enumerate(group_commands, start=1):
            step = f"{label}_{idx}_{Path(command[1]).stem}"
            status, elapsed = run_command(command, dry_run, paths["logs"] / f"{step}.sh")
            timings.append({"step": step, "status": status, "seconds": round(elapsed, 3)})
            commands.append({"step": step, "command": " ".join(command)})
            with command_log.open("a", encoding="utf-8") as handle:
                handle.write(" ".join(command) + "\n")
            if status.startswith("failed"):
                failures.append({"step": step, "status": status})
                break

    pd.DataFrame(commands).to_csv(out_dir / "command_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(timings).to_csv(out_dir / "timing_summary.tsv", sep="\t", index=False)
    pd.DataFrame(failures).to_csv(out_dir / "failed_groups.tsv", sep="\t", index=False)
    groups.to_csv(out_dir / "group_manifest.tsv", sep="\t", index=False)
    summary = {
        "config": str(args.config),
        "dry_run": dry_run,
        "n_groups": int(len(groups)),
        "n_failures": int(len(failures)),
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    (out_dir / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if failures:
        raise SystemExit(f"Workflow completed with {len(failures)} failed step(s); see failed_groups.tsv")


if __name__ == "__main__":
    main()
