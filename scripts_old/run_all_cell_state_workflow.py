#!/usr/bin/env python3
"""Config-driven multi-group cell-state workflow runner."""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import Any

import pandas as pd
import yaml


SCRIPT_DIR = Path(__file__).resolve().parent


def safe_label(value: object) -> str:
    value = str(value).strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return re.sub(r"_+", "_", value).strip("_") or "unknown"


def read_config(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise SystemExit("Workflow config must be a YAML mapping")
    return data


def require(config: dict[str, Any], key: str) -> str:
    value = config.get(key)
    if not value:
        raise SystemExit(f"Config is missing required key: {key}")
    return str(value)


def config_list(config: dict[str, Any], key: str, default: list[str] | None = None) -> list[str]:
    value = config.get(key, default or [])
    if isinstance(value, str):
        return [x for x in value.split() if x]
    if isinstance(value, list):
        return [str(x) for x in value]
    if isinstance(value, dict):
        out: list[str] = []
        for name, option_value in value.items():
            flag = f"--{str(name).replace('_', '-')}"
            if isinstance(option_value, bool):
                if option_value:
                    out.append(flag)
            elif option_value is not None:
                out.extend([flag, str(option_value)])
        return out
    raise SystemExit(f"Config key {key} must be a string, list, or mapping")


def run_command(command: list[str], dry_run: bool, log_path: Path, resume_output: Path | None = None, resume: bool = False) -> dict[str, Any]:
    start = time.time()
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(" ".join(command) + "\n", encoding="utf-8")
    if resume and resume_output is not None and resume_output.exists():
        return {"status": "resumed", "seconds": 0.0, "command": " ".join(command)}
    if dry_run:
        return {"status": "planned", "seconds": 0.0, "command": " ".join(command)}
    result = subprocess.run(command, text=True, capture_output=True)
    (log_path.with_suffix(".stdout.txt")).write_text(result.stdout, encoding="utf-8")
    (log_path.with_suffix(".stderr.txt")).write_text(result.stderr, encoding="utf-8")
    elapsed = round(time.time() - start, 3)
    status = "completed" if result.returncode == 0 else f"failed:{result.returncode}"
    return {"status": status, "seconds": elapsed, "command": " ".join(command)}


def parent_filter(split_by: list[str], group: pd.Series) -> str:
    return ";".join(f"{col}={group[col]}" for col in split_by)


def split_dir(base: Path, split_by: list[str], group: pd.Series) -> Path:
    path = base
    for col in split_by:
        path = path / safe_label(group[col])
    return path


def group_label(split_by: list[str], group: pd.Series) -> str:
    return "__".join(safe_label(group[col]) for col in split_by)


def gmt_for_group(base: Path, split_by: list[str], group: pd.Series) -> Path:
    if len(split_by) == 1:
        return base / "all" / f"{safe_label(group[split_by[0]])}.gmt"
    return base / safe_label(group[split_by[0]]) / f"{safe_label(group[split_by[-1]])}.gmt"


def run_group(group: pd.Series, context: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    config = context["config"]
    split_by = context["split_by"]
    paths = context["paths"]
    dry_run = context["dry_run"]
    resume = context["resume"]
    metadata_path = context["metadata_path"]
    label = group_label(split_by, group)
    group_dir = paths["workflow"] / label
    rank_dir = split_dir(paths["rank_subsets"], split_by, group)
    raw_dir = split_dir(paths["raw_subsets"], split_by, group)
    gmt_path = gmt_for_group(paths["split_gmts"], split_by, group)
    command_log = paths["logs"] / f"{label}.commands.txt"
    records: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []

    if not gmt_path.exists() and not dry_run:
        message = {"group_label": label, "reason": "missing_group_gmt", "gmt": str(gmt_path)}
        if context["skip_groups_without_gmt"]:
            skipped.append(message)
            return {"commands": records, "timings": [], "failures": failures, "skipped": skipped}
        failures.append({"step": f"{label}_missing_gmt", "status": "failed:missing_group_gmt", **message})
        return {"commands": records, "timings": [], "failures": failures, "skipped": skipped}

    scoring_dir = group_dir / "workflow"
    expr_dir = group_dir / "expression"
    gmt_out = group_dir / "state_expression_gmts"
    de_dir = group_dir / "de"
    pigean_dir = group_dir / "pigean"

    group_commands: list[tuple[str, list[str], Path | None]] = []
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
    if config.get("rank_value_type"):
        scoring_cmd.extend(["--rank-value-type", str(config["rank_value_type"])])
    if config.get("expression_kind"):
        scoring_cmd.extend(["--expression-kind", str(config["expression_kind"])])
    if config.get("allow_small_rank_universe"):
        scoring_cmd.append("--allow-small-rank-universe")
    scoring_cmd.extend(config_list(config, "scoring_extra_args"))
    group_commands.append(("run_cmdkp_state_scoring", scoring_cmd, scoring_dir / "cell_state_activity.tsv.gz"))

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
    if config.get("expression_value_type"):
        expr_cmd.extend(["--expression-value-type", str(config["expression_value_type"])])
    if config.get("donor_expression_genes"):
        expr_cmd.extend(["--donor-expression-genes", str(config["donor_expression_genes"])])
    if "write_donor_state_expression" in config:
        expr_cmd.append("--write-donor-state-expression" if config["write_donor_state_expression"] else "--no-write-donor-state-expression")
    if config.get("output_format"):
        expr_cmd.extend(["--output-format", str(config["output_format"])])
    expr_cmd.extend(config_list(config, "expression_extra_args"))
    group_commands.append(("summarize_state_expression", expr_cmd, expr_dir / "all_gene_state_expression_specificity_cp10k.tsv.gz"))

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
    make_gmt_cmd.extend(config_list(config, "gmt_extra_args"))
    group_commands.append(("make_state_expression_gmts", make_gmt_cmd, gmt_out / "gmt_build_summary.tsv"))

    if config.get("pigean"):
        pigean = config["pigean"] or {}
        pigean_cmd = [
            sys.executable,
            str(SCRIPT_DIR / "run_pigean_multi_y_for_state_gmts.py"),
            "--gmt-dir",
            str(gmt_out / "gmt"),
            "--out-dir",
            str(pigean_dir),
            "--multi-y-input",
            str(pigean.get("multi_y_input", pigean.get("traits", ""))),
        ]
        if pigean.get("pigean_command"):
            pigean_cmd.extend(["--pigean-command", str(pigean["pigean_command"])])
        if pigean.get("pigean_bin"):
            pigean_cmd.extend(["--pigean-bin", str(pigean["pigean_bin"])])
        if pigean.get("pigean_command_template"):
            pigean_cmd.extend(["--pigean-command-template", str(pigean["pigean_command_template"])])
        if pigean.get("methods"):
            methods = pigean["methods"]
            pigean_cmd.extend(["--methods", ",".join(methods) if isinstance(methods, list) else str(methods)])
        if pigean.get("extra_args"):
            pigean_cmd.extend(["--extra-args", str(pigean["extra_args"])])
        if pigean.get("dry_run", dry_run):
            pigean_cmd.append("--dry-run")
        group_commands.append(("run_pigean_multi_y_for_state_gmts", pigean_cmd, pigean_dir / "pigean_run_manifest.tsv"))

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
            "--phenotypes",
            str(config["phenotypes"]),
            "--out-dir",
            str(de_dir),
        ]
        donor_expression = expr_dir / "donor_state_weighted_expression_cp10k.tsv.gz"
        if config.get("write_donor_state_expression", True):
            de_cmd.extend(["--donor-state-expression", str(donor_expression)])
        if config.get("phenotype_config"):
            de_cmd.extend(["--phenotype-config", str(config["phenotype_config"])])
        de_cmd.extend(config_list(config, "de_extra_args"))
        group_commands.append(("run_state_phenotype_regression", de_cmd, de_dir / "de_run_summary.json"))

    for idx, (name, command, resume_output) in enumerate(group_commands, start=1):
        step = f"{label}_{idx}_{name}"
        result = run_command(command, dry_run, paths["logs"] / f"{step}.sh", resume_output, resume)
        record = {"step": step, "group_label": label, "status": result["status"], "seconds": result["seconds"], "command": result["command"]}
        records.append(record)
        with command_log.open("a", encoding="utf-8") as handle:
            handle.write(result["command"] + "\n")
        if str(result["status"]).startswith("failed"):
            failures.append(record)
            break

    return {"commands": records, "timings": records, "failures": failures, "skipped": skipped}


def combine_group_pigean(groups: pd.DataFrame, split_by: list[str], workflow_dir: Path, out_dir: Path) -> None:
    combined = []
    for _, group in groups.iterrows():
        label = group_label(split_by, group)
        path = workflow_dir / label / "pigean" / "combined_pigean_state_trait_results.tsv.gz"
        if not path.exists():
            continue
        frame = pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
        frame.insert(0, "group_label", label)
        for col in reversed(split_by):
            frame.insert(1, col, group[col])
        combined.append(frame)
    if combined:
        pd.concat(combined, ignore_index=True).to_csv(out_dir / "combined_pigean_state_trait_results.tsv.gz", sep="\t", index=False, compression="infer")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, required=True)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--skip-groups-without-gmt", action="store_true")
    args = ap.parse_args()

    config = read_config(args.config)
    out_dir = Path(require(config, "out_dir"))
    out_dir.mkdir(parents=True, exist_ok=True)
    dry_run = bool(args.dry_run or config.get("dry_run", False))
    resume = bool(args.resume or config.get("resume", False))
    skip_groups_without_gmt = bool(args.skip_groups_without_gmt or config.get("skip_groups_without_gmt", False))
    jobs = int(config.get("jobs", args.jobs))
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

    paths = {
        "rank_subsets": out_dir / "rank_10x_by_group",
        "raw_subsets": out_dir / "raw_10x_by_group",
        "split_gmts": out_dir / "state_gmts_by_group",
        "logs": out_dir / "command_logs",
        "workflow": out_dir / "groups",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)

    commands: list[dict[str, Any]] = []
    timings: list[dict[str, Any]] = []
    failures: list[dict[str, Any]] = []
    skipped: list[dict[str, Any]] = []

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
    result = run_command(split_cmd, dry_run, paths["logs"] / "split_rank_10x.sh", paths["rank_subsets"] / "split_summary.tsv", resume)
    commands.append({"step": "split_rank_10x", **result})
    timings.append({"step": "split_rank_10x", "status": result["status"], "seconds": result["seconds"]})
    if str(result["status"]).startswith("failed"):
        failures.append({"step": "split_rank_10x", "status": result["status"]})

    raw_10x_dir = str(config.get("raw_10x_dir", config.get("rank_10x_dir", "")))
    if raw_10x_dir and raw_10x_dir != str(config.get("rank_10x_dir")):
        raw_cmd = split_cmd.copy()
        raw_cmd[raw_cmd.index(require(config, "rank_10x_dir"))] = raw_10x_dir
        raw_cmd[raw_cmd.index(str(paths["rank_subsets"]))] = str(paths["raw_subsets"])
        result = run_command(raw_cmd, dry_run, paths["logs"] / "split_raw_10x.sh", paths["raw_subsets"] / "split_summary.tsv", resume)
        commands.append({"step": "split_raw_10x", **result})
        timings.append({"step": "split_raw_10x", "status": result["status"], "seconds": result["seconds"]})
        if str(result["status"]).startswith("failed"):
            failures.append({"step": "split_raw_10x", "status": result["status"]})
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
    ]
    if not skip_groups_without_gmt:
        gmt_cmd.append("--require-all-states")
    result = run_command(gmt_cmd, dry_run, paths["logs"] / "split_state_gmt.sh", paths["split_gmts"] / "split_gmt_manifest.tsv", resume)
    commands.append({"step": "split_state_gmt", **result})
    timings.append({"step": "split_state_gmt", "status": result["status"], "seconds": result["seconds"]})
    if str(result["status"]).startswith("failed"):
        failures.append({"step": "split_state_gmt", "status": result["status"]})

    context = {
        "config": config,
        "split_by": split_by,
        "paths": paths,
        "dry_run": dry_run,
        "resume": resume,
        "metadata_path": metadata_path,
        "skip_groups_without_gmt": skip_groups_without_gmt,
    }
    if jobs <= 1:
        results = [run_group(group, context) for _, group in groups.iterrows()]
    else:
        with ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = [pool.submit(run_group, group, context) for _, group in groups.iterrows()]
            results = [future.result() for future in as_completed(futures)]
    for result_group in results:
        commands.extend(result_group["commands"])
        timings.extend(result_group["timings"])
        failures.extend(result_group["failures"])
        skipped.extend(result_group["skipped"])

    if config.get("pigean") and not dry_run:
        combine_group_pigean(groups, split_by, paths["workflow"], out_dir)

    pd.DataFrame(commands).to_csv(out_dir / "command_manifest.tsv", sep="\t", index=False)
    pd.DataFrame(timings).to_csv(out_dir / "timing_summary.tsv", sep="\t", index=False)
    pd.DataFrame(failures).to_csv(out_dir / "failed_groups.tsv", sep="\t", index=False)
    pd.DataFrame(skipped).to_csv(out_dir / "skipped_groups.tsv", sep="\t", index=False)
    groups.to_csv(out_dir / "group_manifest.tsv", sep="\t", index=False)
    summary = {
        "config": str(args.config),
        "dry_run": dry_run,
        "resume": resume,
        "jobs": jobs,
        "skip_groups_without_gmt": skip_groups_without_gmt,
        "n_groups": int(len(groups)),
        "n_failures": int(len(failures)),
        "n_skipped_groups": int(len(skipped)),
        "pigean_scope": "per_group_per_signature_method",
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    (out_dir / "run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    if failures:
        raise SystemExit(f"Workflow completed with {len(failures)} failed step(s); see failed_groups.tsv")


if __name__ == "__main__":
    main()
