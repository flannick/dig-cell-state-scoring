#!/usr/bin/env python3
"""Run small-model PIGEAN for each LIGER factor as a weighted gene list.

The script reads per-cell-type ``gene_loadings.tsv`` files, rescales each
factor's positive loadings to 0-1 weights, writes one weighted input gene list
per factor, and runs ``pigean betas`` in independent/beta_uncorrected mode
against the bundled small MSigDB and mouse gene-set matrices.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from pathlib import Path

import pandas as pd


def snake_id(value: object) -> str:
    text = str(value).strip().lower()
    text = re.sub(r"[^a-z0-9]+", "_", text)
    return re.sub(r"_+", "_", text).strip("_")


def display_factor(value: str) -> str:
    text = str(value)
    match = re.search(r"(?:program[_-]?)?factor[_-]?([0-9]+)$", text, re.I)
    if match:
        return f"Factor{match.group(1)}"
    match = re.search(r"(?:factor|Factor)[_ ]*([0-9]+)$", text)
    if match:
        return f"Factor{match.group(1)}"
    return text


def read_program_map(path: Path | None) -> dict[str, str]:
    if not path or not path.exists() or path.stat().st_size == 0:
        return {}
    frame = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    if "program_dir" not in frame.columns or "cell_type" not in frame.columns:
        raise SystemExit(f"Program map must contain program_dir and cell_type columns: {path}")
    return dict(zip(frame["program_dir"], frame["cell_type"]))


def read_loadings_long(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", compression="infer")
    if frame.empty:
        return pd.DataFrame(columns=["gene", "factor", "loading"])
    columns_lower = {str(c).lower(): c for c in frame.columns}
    program_col = next((columns_lower[c] for c in ["program_id", "program", "factor", "module", "component"] if c in columns_lower), None)
    gene_col = next((columns_lower[c] for c in ["gene", "gene_symbol", "gene_name"] if c in columns_lower), None)
    loading_col = next((columns_lower[c] for c in ["loading", "weight", "value", "loading_score"] if c in columns_lower), None)
    if program_col and gene_col and loading_col:
        out = frame[[gene_col, program_col, loading_col]].copy()
        out.columns = ["gene", "factor", "loading"]
        out["loading"] = pd.to_numeric(out["loading"], errors="coerce")
        return out.dropna(subset=["gene", "factor", "loading"])

    # LIGER exports often have factor names in the header and an unlabeled gene
    # column in data rows. Pandas cannot read those ragged tables directly.
    with path.open("r", encoding="utf-8") as handle:
        lines = [line.rstrip("\n") for line in handle if line.strip()]
    header = lines[0].split("\t") if lines else []
    first_body = lines[1].split("\t") if len(lines) > 1 else []
    if header and first_body and len(first_body) == len(header) + 1:
        rows = []
        for line in lines[1:]:
            parts = line.split("\t")
            if len(parts) >= len(header) + 1:
                rows.append(parts[: len(header) + 1])
        wide = pd.DataFrame(rows, columns=["gene"] + header)
    else:
        gene_col = frame.columns[0]
        wide = frame.rename(columns={gene_col: "gene"})
    long = wide.melt(id_vars=["gene"], var_name="factor", value_name="loading")
    long["loading"] = pd.to_numeric(long["loading"], errors="coerce")
    return long.dropna(subset=["gene", "factor", "loading"])


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


def normalize_stats(frame: pd.DataFrame, *, tissue: str, cell_type: str, program_dir: str, factor_raw: str, factor: str) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame()
    gene_set_col = pick_col(frame, ["Gene_Set", "gene_set", "gene_set_name", "set", "name"])
    beta_col = pick_col(frame, ["beta", "Beta"])
    beta_unc_col = pick_col(frame, ["beta_uncorrected", "Beta_uncorrected", "beta_uncorrected_orig"])
    p_col = pick_col(frame, ["p", "p_value", "P", "P_Value"])
    label_col = pick_col(frame, ["label", "Label", "source", "model"])
    if gene_set_col is None:
        return pd.DataFrame()
    out = pd.DataFrame(index=frame.index)
    out["tissue"] = tissue
    out["cell_type"] = cell_type
    out["program_dir"] = program_dir
    out["factor"] = factor
    out["factor_raw"] = factor_raw
    out["gene_set_id"] = frame[gene_set_col].astype(str)
    out["gene_set_source"] = frame[label_col].astype(str) if label_col else ""
    out["beta"] = pd.to_numeric(frame[beta_col], errors="coerce") if beta_col else pd.NA
    out["beta_uncorrected"] = pd.to_numeric(frame[beta_unc_col], errors="coerce") if beta_unc_col else pd.NA
    out["p_value"] = pd.to_numeric(frame[p_col], errors="coerce") if p_col else pd.NA
    return out


def factor_weight_table(group: pd.DataFrame, *, top_n: int, min_loading: float, weight_floor: float, weight_ceiling: float) -> pd.DataFrame:
    sub = group.copy()
    sub["loading"] = pd.to_numeric(sub["loading"], errors="coerce")
    sub = sub.dropna(subset=["gene", "loading"])
    sub = sub[sub["loading"].gt(min_loading)]
    sub = sub.sort_values(["loading", "gene"], ascending=[False, True])
    if top_n > 0:
        sub = sub.head(top_n)
    if sub.empty:
        return pd.DataFrame(columns=["Gene", "Weight", "Loading", "Rank"])
    max_loading = float(sub["loading"].max())
    if max_loading <= 0:
        return pd.DataFrame(columns=["Gene", "Weight", "Loading", "Rank"])
    sub = sub.drop_duplicates(subset=["gene"], keep="first").copy()
    sub["Weight"] = (sub["loading"] / max_loading).clip(weight_floor, weight_ceiling)
    sub["Rank"] = range(1, len(sub) + 1)
    return sub.rename(columns={"gene": "Gene", "loading": "Loading"})[["Gene", "Weight", "Loading", "Rank"]]


def write_weighted_gene_list(table: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(path, sep="\t", index=False)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--program-root", type=Path, required=True, help="Directory containing per-cell-type LIGER program subdirectories.")
    ap.add_argument("--program-cell-type-map", type=Path, default=None, help="Optional TSV with program_dir and cell_type columns.")
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--combined-out", type=Path, default=None, help="Default: <out-dir>/factor_small_pigean_results.tsv.gz")
    ap.add_argument("--python", default="../.venv/bin/python", help="Python executable used to run PIGEAN.")
    ap.add_argument("--pythonpath", default="/Users/flannick/codex-workspace/analysis/pigean_optimize/pigean/src")
    ap.add_argument("--msigdb-x-in", default="../resources/pigean/data/small/gene_set_list_msigdb_nohp.txt")
    ap.add_argument("--mouse-x-in", default="../resources/pigean/data/small/gene_set_list_mouse_2024.txt")
    ap.add_argument("--gene-universe-in", default="../resources/pigean/data/reference/NCBI37.3.plink.gene.loc")
    ap.add_argument("--top-n", type=int, default=500, help="Top positive-loading genes per factor; use 0 for all positive loadings.")
    ap.add_argument("--min-loading", type=float, default=0.0)
    ap.add_argument("--weight-floor", type=float, default=0.001)
    ap.add_argument("--weight-ceiling", type=float, default=0.999)
    ap.add_argument("--min-genes", type=int, default=5)
    ap.add_argument("--limit-factors", type=int, default=0, help="Debug limit after sorting factors; 0 means no limit.")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    x_inputs = [Path(args.msigdb_x_in), Path(args.mouse_x_in)]
    required = x_inputs + [Path(args.gene_universe_in), Path(args.python)]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise SystemExit("Missing required PIGEAN inputs: " + ", ".join(missing))

    combined_out = args.combined_out or args.out_dir / "factor_small_pigean_results.tsv.gz"
    program_map = read_program_map(args.program_cell_type_map)
    loading_files = sorted(args.program_root.glob("*/gene_loadings.tsv"))
    if not loading_files:
        raise SystemExit(f"No */gene_loadings.tsv files found under {args.program_root}")

    jobs = []
    for loadings_path in loading_files:
        program_dir = loadings_path.parent.name
        cell_type = program_map.get(program_dir, snake_id(program_dir))
        long = read_loadings_long(loadings_path)
        for factor_raw, group in long.groupby("factor", sort=True):
            factor = display_factor(str(factor_raw))
            weighted = factor_weight_table(
                group,
                top_n=args.top_n,
                min_loading=args.min_loading,
                weight_floor=args.weight_floor,
                weight_ceiling=args.weight_ceiling,
            )
            jobs.append((cell_type, program_dir, str(factor_raw), factor, weighted))
    jobs = sorted(jobs, key=lambda x: (x[0], x[3], x[2]))
    if args.limit_factors > 0:
        jobs = jobs[: args.limit_factors]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    frames = []
    manifest_rows = []
    total = len(jobs)
    print(f"Prepared {total} factor PIGEAN jobs from {len(loading_files)} loading files", flush=True)
    for index, (cell_type, program_dir, factor_raw, factor, weighted) in enumerate(jobs, start=1):
        safe_cell = snake_id(cell_type)
        safe_factor = snake_id(factor)
        run_dir = args.out_dir / "pigean" / safe_cell / safe_factor
        gene_list = args.out_dir / "weighted_gene_lists" / safe_cell / f"{safe_factor}.tsv"
        stats_out = run_dir / "gene_set_stats.beta_uncorrected.out.gz"
        run_dir.mkdir(parents=True, exist_ok=True)
        if len(weighted) < args.min_genes:
            print(f"[{index}/{total}] skipping {cell_type} {factor}: only {len(weighted)} genes after filters", flush=True)
            manifest_rows.append({
                "tissue": args.tissue,
                "cell_type": cell_type,
                "program_dir": program_dir,
                "factor": factor,
                "factor_raw": factor_raw,
                "n_genes": len(weighted),
                "weighted_gene_list": str(gene_list),
                "pigean_run_dir": str(run_dir),
                "status": "skipped_too_few_genes",
            })
            continue
        write_weighted_gene_list(weighted, gene_list)
        cmd = [
            args.python,
            "-m",
            "pigean",
            "betas",
            "--gene-list-in",
            str(gene_list),
            "--gene-list-id-col",
            "Gene",
            "--gene-list-prob-col",
            "Weight",
            "--X-in",
            str(x_inputs[0]),
            "--X-in",
            str(x_inputs[1]),
            "--gene-universe-in",
            str(args.gene_universe_in),
            "--gene-universe-id-col",
            "6",
            "--gene-universe-no-header",
            "--gene-set-stats-out",
            str(stats_out),
            "--params-out",
            str(run_dir / "params.out.gz"),
            "--log-file",
            str(run_dir / "run.log"),
            "--warnings-file",
            str(run_dir / "warnings.log"),
            "--output-detail",
            "debug",
            "--deterministic",
            "--hide-progress",
            "--independent-betas-only",
            "--max-no-write-gene-set-beta-uncorrected",
            "0",
            "--min-gene-set-size",
            "1",
            "--filter-gene-set-p",
            "1",
            "--max-gene-set-read-p",
            "1",
            "--no-filter-negative",
            "--prune-gene-sets",
            "1.1",
            "--weighted-prune-gene-sets",
            "1.1",
        ]
        command_text = "env PYTHONPATH=%s %s\n" % (args.pythonpath, " ".join(cmd))
        command_path = run_dir / "run_command.txt"
        previous_command = command_path.read_text(encoding="utf-8") if command_path.exists() else ""
        needs_run = args.force or not stats_out.exists() or stats_out.stat().st_size == 0 or previous_command != command_text
        if args.dry_run:
            status = "dry_run_would_run" if needs_run else "dry_run_reuse"
            print(f"[{index}/{total}] {status}: {cell_type} {factor} genes={len(weighted)}", flush=True)
        elif needs_run:
            print(f"[{index}/{total}] running PIGEAN: {cell_type} {factor} genes={len(weighted)}", flush=True)
            for stale in [stats_out, run_dir / "params.out.gz"]:
                if stale.exists():
                    stale.unlink()
            command_path.write_text(command_text, encoding="utf-8")
            env = os.environ.copy()
            env["PYTHONPATH"] = args.pythonpath
            with (run_dir / "stdout.log").open("w", encoding="utf-8") as stdout, (run_dir / "stderr.log").open("w", encoding="utf-8") as stderr:
                subprocess.run(cmd, env=env, stdout=stdout, stderr=stderr, check=True)
            status = "ran"
        else:
            print(f"[{index}/{total}] reusing PIGEAN: {cell_type} {factor} genes={len(weighted)}", flush=True)
            status = "reused"
        manifest_rows.append({
            "tissue": args.tissue,
            "cell_type": cell_type,
            "program_dir": program_dir,
            "factor": factor,
            "factor_raw": factor_raw,
            "n_genes": len(weighted),
            "weighted_gene_list": str(gene_list),
            "pigean_run_dir": str(run_dir),
            "gene_set_stats": str(stats_out),
            "status": status,
        })
        if not args.dry_run and stats_out.exists() and stats_out.stat().st_size > 0:
            frame = read_table(stats_out)
            norm = normalize_stats(frame, tissue=args.tissue, cell_type=cell_type, program_dir=program_dir, factor_raw=factor_raw, factor=factor)
            if not norm.empty:
                frames.append(norm)

    manifest = pd.DataFrame(manifest_rows)
    manifest_out = args.out_dir / "factor_small_pigean_manifest.tsv.gz"
    manifest.to_csv(manifest_out, sep="\t", index=False, compression="infer")
    combined_out.parent.mkdir(parents=True, exist_ok=True)
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(
        columns=["tissue", "cell_type", "program_dir", "factor", "factor_raw", "gene_set_id", "gene_set_source", "beta", "beta_uncorrected", "p_value"]
    )
    combined.to_csv(combined_out, sep="\t", index=False, compression="infer")
    print(f"Wrote manifest to {manifest_out}", flush=True)
    print(f"Wrote combined PIGEAN results to {combined_out}", flush=True)


if __name__ == "__main__":
    main()
