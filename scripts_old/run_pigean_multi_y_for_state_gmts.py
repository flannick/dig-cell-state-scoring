#!/usr/bin/env python3
"""Run or stage PIGEAN multi-y analyses for separate state GMT methods."""

from __future__ import annotations

import argparse
import gzip
import json
import shlex
import shutil
import subprocess
from datetime import datetime
from pathlib import Path

import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def write_table(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, compression="infer")


def find_output(directory: Path) -> Path | None:
    candidates = list(directory.glob("*gene_set*stats*.gz")) + list(directory.glob("*.out.gz")) + list(directory.glob("*.tsv.gz"))
    return candidates[0] if candidates else None


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--pigean-bin", default="")
    ap.add_argument("--pigean-command", default="")
    ap.add_argument(
        "--pigean-command-template",
        default="",
        help="Optional template using {pigean}, {gmt}, {multi_y}, {out_dir}, {out}, {method}, and {extra_args}",
    )
    ap.add_argument("--multi-y-input", default="")
    ap.add_argument("--traits", default="", help="Alias for --multi-y-input when the local PIGEAN command expects trait input")
    ap.add_argument("--gmt-dir", type=Path, default=Path(""))
    ap.add_argument("--gmt", action="append", default=[])
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--methods", default="original_markers,top_absolute_expression,top_specific_fc,top_specific_logp")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--extra-args", default="")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    command = args.pigean_command or args.pigean_bin or "pigean"
    multi_y = args.multi_y_input or args.traits
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    gmt_by_method: dict[str, Path] = {}
    for gmt in args.gmt:
        path = Path(gmt)
        gmt_by_method[path.stem] = path
    if args.gmt_dir:
        for method in methods:
            path = args.gmt_dir / f"{method}.gmt"
            if path.exists():
                gmt_by_method[method] = path
    if not gmt_by_method:
        raise SystemExit("No GMT files found")
    command_available = bool(shutil.which(shlex.split(command)[0]))
    manifest = []
    combined = []
    for method in methods:
        if method not in gmt_by_method:
            continue
        method_dir = args.out_dir / method
        method_dir.mkdir(parents=True, exist_ok=True)
        out_prefix = method_dir / "gene_set_stats.out.gz"
        if args.pigean_command_template:
            rendered = args.pigean_command_template.format(
                pigean=command,
                gmt=str(gmt_by_method[method]),
                multi_y=multi_y,
                out_dir=str(method_dir),
                out=str(out_prefix),
                method=method,
                extra_args=args.extra_args,
            )
            cmd = shlex.split(rendered)
        else:
            cmd = shlex.split(command) + ["betas", "--gmt-in", str(gmt_by_method[method]), "--multi-y-in", multi_y, "--out", str(out_prefix)]
            if args.extra_args:
                cmd.extend(shlex.split(args.extra_args))
        cmd_text = " ".join(shlex.quote(x) for x in cmd)
        (method_dir / "run_command.txt").write_text(cmd_text + "\n", encoding="utf-8")
        run = {
            "method": method,
            "gmt": str(gmt_by_method[method]),
            "command": cmd_text,
            "ran": False,
            "returncode": None,
            "timestamp": datetime.now().isoformat(timespec="seconds"),
        }
        if command_available and not args.dry_run and multi_y:
            proc = subprocess.run(cmd, text=True, capture_output=True)
            (method_dir / "stdout.txt").write_text(proc.stdout, encoding="utf-8")
            (method_dir / "stderr.txt").write_text(proc.stderr, encoding="utf-8")
            run["ran"] = True
            run["returncode"] = proc.returncode
            if proc.returncode != 0:
                run["error"] = "pigean_command_failed"
            out = find_output(method_dir)
            if out is not None:
                try:
                    frame = read_table(out)
                    frame["signature_method"] = method
                    combined.append(frame)
                    run["output"] = str(out)
                except Exception as exc:  # pragma: no cover - summary only
                    run["output_read_error"] = str(exc)
        else:
            shell = method_dir / "run_pigean.sh"
            shell.write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + cmd_text + "\n", encoding="utf-8")
            shell.chmod(0o755)
            run["runnable_script"] = str(shell)
            run["error"] = "pigean_not_run_command_unavailable_or_dry_run_or_missing_multi_y"
        (method_dir / "run_summary.json").write_text(json.dumps(run, indent=2), encoding="utf-8")
        manifest.append(run)
    pd.DataFrame(manifest).to_csv(args.out_dir / "pigean_run_manifest.tsv", sep="\t", index=False)
    (args.out_dir / "pigean_run_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    if combined:
        write_table(pd.concat(combined, ignore_index=True), args.out_dir / "combined_pigean_state_trait_results.tsv.gz")


if __name__ == "__main__":
    main()
