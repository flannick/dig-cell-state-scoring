#!/usr/bin/env python3
"""Run PIGEAN betas for API-minimal tissue pipeline GMTs and combine results."""

from __future__ import annotations

import argparse
import gzip
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd


def read_gmt_names(path: Path) -> list[str]:
    names = []
    with path.open() as handle:
        for line in handle:
            if line.strip():
                names.append(line.split("\t", 1)[0])
    return names


def write_pigean_x_from_gmt(gmt_path: Path, out_path: Path) -> Path:
    """Write PIGEAN --X-in format from GMT by dropping the description column."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with gmt_path.open() as src, out_path.open("w", encoding="utf-8") as dst:
        for line in src:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            genes = [g for g in parts[2:] if g]
            if genes:
                dst.write("\t".join([parts[0]] + genes) + "\n")
    return out_path


def display_factor(value: str) -> str:
    text = str(value)
    m = re.search(r"(?:program[_-]?)?factor[_-]?([0-9]+)$", text, re.I)
    if m:
        return f"Factor{m.group(1)}"
    m = re.search(r"(?:factor|Factor)[_ ]*([0-9]+)$", text)
    if m:
        return f"Factor{m.group(1)}"
    return text


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists() or path.stat().st_size == 0:
        return pd.DataFrame()
    try:
        return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)
    except pd.errors.EmptyDataError:
        return pd.DataFrame()


def pick_col(frame: pd.DataFrame, names: list[str]) -> str | None:
    by_lower = {str(c).lower(): c for c in frame.columns}
    for name in names:
        if name.lower() in by_lower:
            return by_lower[name.lower()]
    return None


def generate_hp_exomes_blacklist(multi_y_in: Path, out_path: Path, pheno_col: str = "Trait_Internal") -> Path:
    opener = gzip.open if str(multi_y_in).endswith(".gz") else open
    traits: set[str] = set()
    with opener(multi_y_in, "rt") as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if pheno_col not in header:
            raise SystemExit(f"Cannot generate trait blacklist: {pheno_col} not found in {multi_y_in}")
        idx = header.index(pheno_col)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if idx >= len(parts):
                continue
            trait = parts[idx]
            if trait.startswith("HP_") or trait.startswith("exomes_"):
                traits.add(trait)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("".join(f"{trait}\n" for trait in sorted(traits)), encoding="utf-8")
    return out_path


def resolve_trait_blacklist(path_text: str, multi_y_in: Path, out_dir: Path) -> Path:
    if path_text in {"", "auto", "AUTO", "hp_exomes", "HP_EXOMES"}:
        return generate_hp_exomes_blacklist(multi_y_in, out_dir / "trait_blacklist_hp_exomes.txt")
    blacklist = Path(path_text)
    if blacklist.exists():
        return blacklist
    fallback = Path("../blanc_screen/results/pigean_all_cell_states_full_universe_default_no_hpo_exomes/trait_blacklist_exomes_hp.txt")
    if fallback.exists():
        return fallback
    return generate_hp_exomes_blacklist(multi_y_in, out_dir / "trait_blacklist_hp_exomes.txt")


def normalize_one(frame: pd.DataFrame, *, gmt: Path, cell_type: str, tissue: str, dataset: str, model: str, kind: str) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame()
    gene_set_col = pick_col(frame, ["Gene_Set", "gene_set", "gene_set_name", "state_name", "factor", "set", "name"])
    trait_col = pick_col(frame, ["trait", "Trait", "phenotype", "Trait_Internal", "y"])
    beta_col = pick_col(frame, ["beta", "Beta"])
    beta_unc_col = pick_col(frame, ["beta_uncorrected", "Beta_uncorrected", "beta_uncorrected_orig"])
    if gene_set_col is None:
        return pd.DataFrame()
    out = pd.DataFrame()
    out["gene_set_id"] = frame[gene_set_col].astype(str)
    out["trait"] = frame[trait_col].astype(str) if trait_col else ""
    out["beta"] = pd.to_numeric(frame[beta_col], errors="coerce") if beta_col else pd.NA
    out["beta_uncorrected"] = pd.to_numeric(frame[beta_unc_col], errors="coerce") if beta_unc_col else pd.NA
    if kind == "program":
        out["dataset"] = dataset
        out["cell_type"] = cell_type
        out["model"] = model
        out["factor"] = out["gene_set_id"].map(display_factor)
        return out[["dataset", "cell_type", "model", "factor", "trait", "beta", "beta_uncorrected"]]
    out["tissue"] = tissue
    out["cell_type"] = cell_type
    out["state_name"] = out["gene_set_id"]
    return out[["tissue", "cell_type", "state_name", "trait", "beta", "beta_uncorrected"]]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gmt-dir", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--combined-out", type=Path, required=True)
    ap.add_argument("--kind", choices=["curated", "program"], required=True)
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--dataset", default="")
    ap.add_argument("--model", default="")
    ap.add_argument("--python", default=sys.executable)
    ap.add_argument("--pythonpath", default="/Users/flannick/codex-workspace/analysis/pigean_optimize/pigean/src")
    ap.add_argument("--multi-y-in", default="../resources/pigean/data/large/all.gene_stats.large.gt1.out.gz")
    ap.add_argument("--trait-blacklist-in", default="auto", help="Trait blacklist path, or auto to generate HP_/exomes_ blacklist from --multi-y-in")
    ap.add_argument("--gene-universe-in", default="../resources/pigean/data/reference/NCBI37.3.plink.gene.loc")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.combined_out.parent.mkdir(parents=True, exist_ok=True)
    multi_y_in = Path(args.multi_y_in)
    gene_universe_in = Path(args.gene_universe_in)
    missing = [str(x) for x in [multi_y_in, gene_universe_in] if not x.exists()]
    if missing:
        raise SystemExit("Missing PIGEAN input files: " + ", ".join(missing))
    blacklist = resolve_trait_blacklist(args.trait_blacklist_in, multi_y_in, args.combined_out.parent)
    if not blacklist.exists():
        raise SystemExit(f"Missing PIGEAN trait blacklist: {blacklist}")
    print(f"Using PIGEAN trait blacklist {blacklist}", flush=True)

    frames = []
    for gmt in sorted(args.gmt_dir.glob("*.gmt")):
        if not gmt.stat().st_size or not read_gmt_names(gmt):
            continue
        cell_type = gmt.stem
        run_dir = args.out_dir / cell_type
        stats_out = run_dir / "gene_set_stats.debug.out.gz"
        run_dir.mkdir(parents=True, exist_ok=True)
        pigean_x = write_pigean_x_from_gmt(gmt, run_dir / "input.x.tsv")
        cmd = [
                args.python, "-m", "pigean", "betas",
                "--X-in", str(pigean_x),
                "--multi-y-in", str(multi_y_in),
                "--multi-y-id-col", "Gene",
                "--multi-y-pheno-col", "Trait_Internal",
                "--multi-y-log-bf-col", "Direct",
                "--multi-y-combined-col", "Combined",
                "--multi-y-prior-col", "Indirect",
                "--multi-y-trait-blacklist-in", str(blacklist),
                "--gene-universe-in", str(gene_universe_in),
                "--gene-universe-id-col", "6",
                "--gene-universe-no-header",
                "--gene-set-stats-out", str(stats_out),
                "--params-out", str(run_dir / "params.out.gz"),
                "--log-file", str(run_dir / "run.log"),
                "--warnings-file", str(run_dir / "warnings.log"),
                "--output-detail", "debug",
                "--deterministic",
                "--hide-progress",
                "--min-gene-set-size", "1",
                "--filter-gene-set-p", "1",
                "--max-gene-set-read-p", "1",
                "--no-filter-negative",
                "--prune-gene-sets", "1.1",
                "--weighted-prune-gene-sets", "1.1",
        ]
        command_text = "env PYTHONPATH=%s %s\n" % (args.pythonpath, " ".join(cmd))
        command_path = run_dir / "run_command.txt"
        previous_command = command_path.read_text(encoding="utf-8") if command_path.exists() else ""
        needs_run = args.force or not stats_out.exists() or stats_out.stat().st_size == 0 or previous_command != command_text
        if needs_run:
            for stale in [stats_out, run_dir / "params.out.gz"]:
                if stale.exists():
                    stale.unlink()
            command_path.write_text(command_text, encoding="utf-8")
            env = os.environ.copy()
            env["PYTHONPATH"] = args.pythonpath
            with (run_dir / "stdout.log").open("w") as stdout, (run_dir / "stderr.log").open("w") as stderr:
                subprocess.run(cmd, env=env, stdout=stdout, stderr=stderr, check=True)
        frame = read_table(stats_out)
        norm = normalize_one(frame, gmt=gmt, cell_type=cell_type, tissue=args.tissue, dataset=args.dataset, model=args.model, kind=args.kind)
        if not norm.empty:
            frames.append(norm)
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    if combined.empty:
        cols = ["dataset", "cell_type", "model", "factor", "trait", "beta", "beta_uncorrected"] if args.kind == "program" else ["tissue", "cell_type", "state_name", "trait", "beta", "beta_uncorrected"]
        combined = pd.DataFrame(columns=cols)
    combined.to_csv(args.combined_out, sep="\t", index=False, compression="infer")


if __name__ == "__main__":
    main()
