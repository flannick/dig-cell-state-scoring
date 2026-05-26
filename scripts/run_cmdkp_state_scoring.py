#!/usr/bin/env python3
"""Run the general CMDKP cell-state scoring and calling workflow.

The runner consumes an exported cell x gene expression table plus cell metadata.
It scores biological state GMTs and the auxiliary bad-cell QC GMT with local
rank-based UCell- and AUCell-style statistics and writes the standard CMDKP
output bundle.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sys
import time
from datetime import datetime
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from scipy.io import mmread
from scipy import stats
from scipy.ndimage import gaussian_filter1d


PROCESS_TOKENS = {
    "upr",
    "er_stress",
    "interferon",
    "ifn",
    "mhc",
    "oxidative",
    "hypoxia",
    "proliferation",
    "proliferative",
    "cell_cycle",
    "cycling",
    "fibrosis",
    "remodeling",
}
COMPOSITE_TOKENS = {
    "dedifferentiation",
    "dedifferentiated",
    "low_identity",
    "disallowed",
    "senescence",
    "senescence_like",
    "doublet",
    "contamination",
}
HEMOGLOBIN_PREFIXES = ("HBA", "HBB", "HBD", "HBG", "HBM", "HBQ")
IDENTITY_TOKENS = {"identity", "marker", "canonical", "lineage"}
FUNCTION_TOKENS = {"function", "secretory", "metabolic", "mature"}
RARE_PROCESS_TOKENS = {"apoptosis", "cell_death", "necrosis", "rare"}
QC_TOKENS = {"qc", "ambient", "contamination", "doublet", "mitochondrial", "ribosomal", "hemoglobin", "platelet"}


@dataclass
class GeneSet:
    name: str
    description: str
    genes: list[str]
    meta: dict[str, str]


@dataclass
class SparseRankUniverse:
    matrix: sparse.csr_matrix
    cells: list[str]
    genes: list[str]


def norm_name(value: str) -> str:
    value = str(value).strip().lower()
    value = re.sub(r"[^a-z0-9]+", "_", value)
    return re.sub(r"_+", "_", value).strip("_")


def read_table(path: str) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", compression="infer", low_memory=False)


def open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8") if str(path).endswith(".gz") else open(path, "r", encoding="utf-8")


def write_table(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, compression="infer")


def empty_table(columns: list[str]) -> pd.DataFrame:
    return pd.DataFrame({column: pd.Series(dtype="object") for column in columns})


def output_path(args: argparse.Namespace, attr: str, default_name: str) -> Path:
    explicit = getattr(args, attr)
    if explicit:
        return Path(explicit)
    return Path(args.out_dir) / default_name


def read_expression_input(path: str) -> pd.DataFrame:
    frame = read_table(path)
    if {"cell_id", "gene", "expression"}.issubset(frame.columns):
        return frame[["cell_id", "gene", "expression"]].copy()
    if "cell_id" not in frame.columns:
        first = frame.columns[0]
        frame = frame.rename(columns={first: "cell_id"})
    long = frame.melt(id_vars=["cell_id"], var_name="gene", value_name="expression")
    long = long.loc[pd.to_numeric(long["expression"], errors="coerce").fillna(0) != 0].copy()
    return long


def load_threshold_yaml(path: str) -> dict[str, float]:
    if not path:
        return {}
    with open_text(path) as handle:
        obj = yaml.safe_load(handle) or {}
    if isinstance(obj, dict) and "states" in obj:
        obj = obj["states"]
    out: dict[str, float] = {}
    if isinstance(obj, dict):
        for key, value in obj.items():
            if isinstance(value, dict):
                value = value.get("threshold", value.get("aucell_score"))
            if value is not None:
                out[str(key)] = float(value)
    return out


def load_state_class_config(path: str) -> dict[str, str]:
    if not path:
        return {}
    with open_text(path) as handle:
        obj = yaml.safe_load(handle) or {}
    if isinstance(obj, dict) and "states" in obj:
        obj = obj["states"]
    out: dict[str, str] = {}
    if isinstance(obj, dict):
        for key, value in obj.items():
            if isinstance(value, dict):
                value = value.get("state_class", value.get("class"))
            if value:
                out[str(key)] = str(value)
    return out


def read_gene_list(path: str) -> list[str]:
    if not path:
        return []
    genes: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            value = line.strip()
            if not value or value.startswith("#"):
                continue
            genes.append(value.split()[0])
    return list(dict.fromkeys(genes))


def read_one_column(path: str) -> list[str]:
    values: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            value = line.strip()
            if value:
                values.append(value.split("\t")[0])
    return values


def read_10x_features(path: str) -> list[str]:
    genes: list[str] = []
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2 and parts[1]:
                genes.append(parts[1])
            elif parts and parts[0]:
                genes.append(parts[0])
    return genes


def resolve_10x_path(directory: str, names: list[str]) -> str:
    base = Path(directory)
    for name in names:
        path = base / name
        if path.exists():
            return str(path)
    raise SystemExit(f"Could not find any of {', '.join(names)} in {directory}")


def load_sparse_rank_universe(args: argparse.Namespace, metadata: pd.DataFrame) -> SparseRankUniverse | None:
    matrix_path = args.rank_matrix_mtx
    genes_path = args.rank_genes
    cells_path = args.rank_cells
    if args.rank_10x_dir:
        matrix_path = matrix_path or resolve_10x_path(args.rank_10x_dir, ["matrix.mtx.gz", "matrix.mtx"])
        genes_path = genes_path or resolve_10x_path(args.rank_10x_dir, ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"])
        cells_path = cells_path or resolve_10x_path(args.rank_10x_dir, ["barcodes.tsv.gz", "barcodes.tsv"])
    if not matrix_path:
        return None
    if not genes_path or not cells_path:
        raise SystemExit("--rank-matrix-mtx requires --rank-genes and --rank-cells, or use --rank-10x-dir")
    mat = mmread(matrix_path).tocsr()
    genes = read_10x_features(genes_path)
    cells = read_one_column(cells_path)
    if mat.shape == (len(genes), len(cells)):
        mat = mat.T.tocsr()
    elif mat.shape != (len(cells), len(genes)):
        raise SystemExit(
            f"Rank matrix shape {mat.shape} does not match cells x genes "
            f"({len(cells)}, {len(genes)}) or genes x cells ({len(genes)}, {len(cells)})"
        )
    cell_to_pos = {cell: i for i, cell in enumerate(cells)}
    missing = [cell for cell in metadata["cell_id"] if cell not in cell_to_pos]
    if missing:
        raise SystemExit(f"Rank universe is missing {len(missing)} metadata cell IDs")
    order = [cell_to_pos[cell] for cell in metadata["cell_id"]]
    mat = mat[order, :].tocsr()
    return SparseRankUniverse(matrix=mat, cells=metadata["cell_id"].astype(str).tolist(), genes=genes)


def load_sparse_10x_dir(directory: str, metadata: pd.DataFrame) -> SparseRankUniverse | None:
    if not directory:
        return None
    sparse_args = argparse.Namespace(rank_10x_dir=directory, rank_matrix_mtx="", rank_genes="", rank_cells="")
    return load_sparse_rank_universe(sparse_args, metadata)


def resolve_query_genes(args: argparse.Namespace, matrix: pd.DataFrame) -> list[str]:
    requested = read_gene_list(args.query_genes) if args.query_genes else []
    if args.query_gene:
        requested.extend(g.strip() for g in args.query_gene.split(",") if g.strip())
    if not requested:
        return list(matrix.columns)
    requested = list(dict.fromkeys(requested))
    present = [gene for gene in requested if gene in matrix.columns]
    missing = [gene for gene in requested if gene not in matrix.columns]
    if missing:
        print(f"Warning: {len(missing)} query genes were not found in the expression matrix", file=sys.stderr)
    if not present:
        raise SystemExit("No query genes were found in the expression matrix")
    return present


def apply_metadata_filter(metadata: pd.DataFrame, filter_expr: str) -> pd.DataFrame:
    if not filter_expr:
        return metadata
    keep = pd.Series(True, index=metadata.index)
    for term in filter_expr.split(";"):
        if not term:
            continue
        if "=" not in term:
            raise SystemExit("--parent-cell-filter must have form column=value1,value2 or semicolon-separated filters")
        column, values = term.split("=", 1)
        column = column.strip()
        allowed = {value.strip() for value in values.split(",") if value.strip()}
        if column not in metadata.columns:
            raise SystemExit(f"--parent-cell-filter column not found in metadata: {column}")
        if not allowed:
            raise SystemExit("--parent-cell-filter did not include any values")
        keep &= metadata[column].astype(str).isin(allowed)
    out = metadata.loc[keep].copy()
    if out.empty:
        raise SystemExit("--parent-cell-filter removed all metadata rows")
    return out


def bh_fdr(pvalues: pd.Series) -> pd.Series:
    p = pd.to_numeric(pvalues, errors="coerce")
    out = pd.Series(np.nan, index=p.index)
    valid = p.notna()
    if valid.sum() == 0:
        return out
    order = p[valid].sort_values().index
    ranked = p.loc[order].to_numpy()
    n = len(ranked)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out.loc[order] = np.clip(q, 0, 1)
    return out


def parse_description(desc: str) -> dict[str, str]:
    meta: dict[str, str] = {}
    for part in str(desc).split(";"):
        if "=" in part:
            key, value = part.split("=", 1)
            meta[key.strip()] = value.strip()
    return meta


def read_gmt(path: str, regex: str = "") -> list[GeneSet]:
    pattern = re.compile(regex) if regex else None
    sets: list[GeneSet] = []
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            name = parts[0]
            if pattern and not pattern.search(name):
                continue
            genes = list(dict.fromkeys(g for g in parts[2:] if g))
            sets.append(GeneSet(name=name, description=parts[1], genes=genes, meta=parse_description(parts[1])))
    if not sets:
        raise SystemExit(f"No GMT rows found in {path}")
    return sets


def load_state_manifest(path: str) -> dict[str, dict[str, object]]:
    if not path:
        return {}
    frame = read_table(path)
    if "state_name" not in frame.columns:
        if "name" in frame.columns:
            frame = frame.rename(columns={"name": "state_name"})
        else:
            raise SystemExit("--state-manifest must include a state_name column")
    manifest: dict[str, dict[str, object]] = {}
    bool_cols = {"is_composite_required", "is_qc", "allow_hard_call"}
    for _, row in frame.dropna(subset=["state_name"]).iterrows():
        values = {k: v for k, v in row.to_dict().items() if pd.notna(v)}
        for col in bool_cols & set(values):
            values[col] = str(values[col]).strip().lower() in {"1", "true", "yes", "y"}
        manifest[str(values["state_name"])] = values
    return manifest


def validate_state_manifest(
    manifest: dict[str, dict[str, object]],
    gene_sets: list[GeneSet],
    require: bool = False,
) -> None:
    if not require:
        return
    required_cols = {"state_name", "tissue", "cell_type", "state_class", "is_composite_required"}
    if not manifest:
        raise SystemExit("--require-state-manifest was set but no --state-manifest was supplied")
    missing_states = [gene_set.name for gene_set in gene_sets if gene_set.name not in manifest]
    if missing_states:
        preview = ", ".join(missing_states[:10])
        suffix = "..." if len(missing_states) > 10 else ""
        raise SystemExit(f"State manifest is missing {len(missing_states)} GMT state(s): {preview}{suffix}")
    bad_rows = []
    for state_name, row in manifest.items():
        missing = [col for col in required_cols if col == "state_name" and not state_name or col != "state_name" and col not in row]
        if missing:
            bad_rows.append(f"{state_name}: {','.join(missing)}")
    if bad_rows:
        preview = "; ".join(bad_rows[:10])
        suffix = "..." if len(bad_rows) > 10 else ""
        raise SystemExit(f"State manifest rows are missing required production columns: {preview}{suffix}")


def apply_state_manifest(gene_sets: list[GeneSet], manifest: dict[str, dict[str, object]]) -> list[GeneSet]:
    if not manifest:
        return gene_sets
    for gene_set in gene_sets:
        if gene_set.name in manifest:
            gene_set.meta.update({k: str(v) for k, v in manifest[gene_set.name].items() if k != "state_name"})
            gene_set.meta["manifest_present"] = "true"
        else:
            gene_set.meta["manifest_present"] = "false"
    return gene_sets


class StepTimer:
    def __init__(self) -> None:
        self.started = time.time()
        self.last = self.started
        self.rows: list[dict[str, object]] = []

    def mark(self, step: str, **extra: object) -> None:
        now = time.time()
        row = {
            "step": step,
            "seconds_since_previous": round(now - self.last, 3),
            "seconds_since_start": round(now - self.started, 3),
        }
        row.update(extra)
        self.rows.append(row)
        self.last = now


def harmonize_expression(
    expression: pd.DataFrame,
    expression_kind: str,
    gene_map_path: str,
    duplicate_method: str,
) -> tuple[pd.DataFrame, dict[str, str]]:
    required = {"cell_id", "gene", "expression"}
    missing = sorted(required - set(expression.columns))
    if missing:
        raise SystemExit(f"Expression table is missing required column(s): {', '.join(missing)}")
    expr = expression[["cell_id", "gene", "expression"]].copy()
    expr["expression"] = pd.to_numeric(expr["expression"], errors="coerce").fillna(0.0)
    mapping_info = {"gene_id_type": "assumed_hgnc_symbol", "duplicate_collapse_method": "none"}

    if gene_map_path:
        gene_map = read_table(gene_map_path)
        for col in ["gene_id", "gene_symbol"]:
            if col not in gene_map.columns:
                raise SystemExit(f"Gene map is missing required column: {col}")
        gene_map = gene_map[["gene_id", "gene_symbol"]].dropna().drop_duplicates("gene_id")
        expr = expr.merge(gene_map, left_on="gene", right_on="gene_id", how="left")
        expr["gene"] = expr["gene_symbol"].fillna(expr["gene"])
        expr = expr.drop(columns=["gene_id", "gene_symbol"])
        mapping_info["gene_id_type"] = "mapped_to_hgnc_symbol"

    expr["gene"] = expr["gene"].astype(str)
    if duplicate_method == "auto":
        duplicate_method = "sum" if expression_kind == "raw_counts" else "max"
    mapping_info["duplicate_collapse_method"] = duplicate_method
    if duplicate_method == "sum":
        expr = expr.groupby(["cell_id", "gene"], as_index=False)["expression"].sum()
    elif duplicate_method == "mean":
        expr = expr.groupby(["cell_id", "gene"], as_index=False)["expression"].mean()
    elif duplicate_method == "max":
        expr = expr.groupby(["cell_id", "gene"], as_index=False)["expression"].max()
    else:
        raise SystemExit("--duplicate-collapse must be auto, sum, mean, or max")
    return expr, mapping_info


def expression_matrix(expression: pd.DataFrame, cells: pd.Series) -> pd.DataFrame:
    matrix = expression.pivot(index="cell_id", columns="gene", values="expression")
    matrix = matrix.reindex(cells).fillna(0.0)
    return matrix.astype(float)


def rank_matrix_for_ucell(matrix: pd.DataFrame, max_rank: int) -> pd.DataFrame:
    ranks = matrix.rank(axis=1, ascending=False, method="average")
    ranks = ranks.where(ranks <= max_rank, max_rank + 1)
    return ranks


def rank_matrix_for_aucell(matrix: pd.DataFrame) -> pd.DataFrame:
    return matrix.rank(axis=1, ascending=False, method="average")


def ucell_score_from_ranks(ranks: pd.DataFrame, genes: list[str], max_rank: int) -> pd.Series:
    present = [g for g in genes if g in ranks.columns]
    n = len(present)
    if n == 0:
        return pd.Series(np.nan, index=ranks.index)
    max_u = n * max_rank - (n * (n + 1)) / 2
    if max_u <= 0:
        return pd.Series(np.nan, index=ranks.index)
    rank_sum = ranks[present].sum(axis=1)
    u_stat = rank_sum - (n * (n + 1)) / 2
    return (1 - (u_stat / max_u)).clip(0, 1)


def aucell_score_from_ranks(ranks: pd.DataFrame, genes: list[str], auc_max_rank: int) -> pd.Series:
    present = [g for g in genes if g in ranks.columns]
    if not present or auc_max_rank < 1:
        return pd.Series(np.nan, index=ranks.index)
    top = ranks[present].where(ranks[present] <= auc_max_rank)
    recovery = (auc_max_rank - top + 1).clip(lower=0).sum(axis=1)
    max_recovery = len(present) * auc_max_rank
    if max_recovery <= 0:
        return pd.Series(np.nan, index=ranks.index)
    return (recovery / max_recovery).clip(0, 1)


def sparse_rank_scores_for_gene_set(
    universe: SparseRankUniverse,
    genes: list[str],
    auc_max_rank: int,
    max_rank: int,
) -> tuple[pd.Series, pd.Series, list[str]]:
    gene_to_pos = {gene: i for i, gene in enumerate(universe.genes)}
    present = [gene for gene in genes if gene in gene_to_pos]
    marker_positions = {gene_to_pos[gene] for gene in present}
    auc_scores = np.full(len(universe.cells), np.nan, dtype=float)
    ucell_scores = np.full(len(universe.cells), np.nan, dtype=float)
    if not present:
        return pd.Series(auc_scores, index=universe.cells), pd.Series(ucell_scores, index=universe.cells), present
    max_u = len(present) * max_rank - (len(present) * (len(present) + 1)) / 2
    for cell_idx in range(universe.matrix.shape[0]):
        start, end = universe.matrix.indptr[cell_idx], universe.matrix.indptr[cell_idx + 1]
        cols = universe.matrix.indices[start:end]
        vals = universe.matrix.data[start:end]
        if len(cols) == 0:
            auc_scores[cell_idx] = 0.0
            ucell_scores[cell_idx] = 0.0 if max_u > 0 else np.nan
            continue
        order = np.lexsort((cols, -vals))
        ranked_cols = cols[order]
        rank_by_col: dict[int, int] = {}
        auc_recovery = 0.0
        for rank, col in enumerate(ranked_cols[: max(max_rank, auc_max_rank)], start=1):
            if col in marker_positions:
                if rank <= max_rank:
                    rank_by_col[col] = rank
                if rank <= auc_max_rank:
                    auc_recovery += auc_max_rank - rank + 1
        auc_scores[cell_idx] = auc_recovery / (len(present) * auc_max_rank) if auc_max_rank > 0 else np.nan
        if max_u > 0:
            rank_sum = sum(rank_by_col.get(col, max_rank + 1) for col in marker_positions)
            u_stat = rank_sum - (len(present) * (len(present) + 1)) / 2
            ucell_scores[cell_idx] = np.clip(1 - (u_stat / max_u), 0, 1)
    return (
        pd.Series(auc_scores, index=universe.cells),
        pd.Series(ucell_scores, index=universe.cells),
        present,
    )


def sparse_rank_scores_for_gene_sets(
    universe: SparseRankUniverse,
    gene_sets: list[GeneSet],
    auc_max_rank: int,
    max_rank: int,
) -> dict[str, tuple[pd.Series, pd.Series, list[str]]]:
    """Score many gene sets from one sparse per-cell ranking pass."""
    gene_to_pos = {gene: i for i, gene in enumerate(universe.genes)}
    set_positions: dict[str, set[int]] = {}
    pos_to_sets: dict[int, list[str]] = {}
    present_by_set: dict[str, list[str]] = {}
    for gene_set in gene_sets:
        present = [gene for gene in gene_set.genes if gene in gene_to_pos]
        positions = {gene_to_pos[gene] for gene in present}
        present_by_set[gene_set.name] = present
        set_positions[gene_set.name] = positions
        for position in positions:
            pos_to_sets.setdefault(position, []).append(gene_set.name)

    auc_values = {
        gene_set.name: np.full(len(universe.cells), np.nan if not set_positions[gene_set.name] else 0.0, dtype=float)
        for gene_set in gene_sets
    }
    rank_sums = {
        gene_set.name: np.zeros(len(universe.cells), dtype=float)
        for gene_set in gene_sets
        if set_positions[gene_set.name]
    }
    hit_counts = {
        gene_set.name: np.zeros(len(universe.cells), dtype=np.int32)
        for gene_set in gene_sets
        if set_positions[gene_set.name]
    }
    top_rank = max(max_rank, auc_max_rank)
    for cell_idx in range(universe.matrix.shape[0]):
        start, end = universe.matrix.indptr[cell_idx], universe.matrix.indptr[cell_idx + 1]
        cols = universe.matrix.indices[start:end]
        vals = universe.matrix.data[start:end]
        if len(cols) == 0:
            continue
        if len(cols) > top_rank:
            top_unsorted = np.argpartition(-vals, top_rank - 1)[:top_rank]
            order = top_unsorted[np.lexsort((cols[top_unsorted], -vals[top_unsorted]))]
        else:
            order = np.lexsort((cols, -vals))
        for rank, col in enumerate(cols[order], start=1):
            for set_name in pos_to_sets.get(col, []):
                if rank <= auc_max_rank:
                    auc_values[set_name][cell_idx] += auc_max_rank - rank + 1
                if rank <= max_rank:
                    rank_sums[set_name][cell_idx] += rank
                    hit_counts[set_name][cell_idx] += 1

    out: dict[str, tuple[pd.Series, pd.Series, list[str]]] = {}
    for gene_set in gene_sets:
        present = present_by_set[gene_set.name]
        n_present = len(present)
        if n_present == 0:
            out[gene_set.name] = (
                pd.Series(auc_values[gene_set.name], index=universe.cells),
                pd.Series(np.nan, index=universe.cells),
                present,
            )
            continue
        auc = auc_values[gene_set.name] / (n_present * auc_max_rank) if auc_max_rank > 0 else np.full(len(universe.cells), np.nan)
        max_u = n_present * max_rank - (n_present * (n_present + 1)) / 2
        if max_u > 0:
            rank_sum = rank_sums[gene_set.name] + (n_present - hit_counts[gene_set.name]) * (max_rank + 1)
            u_stat = rank_sum - (n_present * (n_present + 1)) / 2
            ucell = np.clip(1 - (u_stat / max_u), 0, 1)
        else:
            ucell = np.full(len(universe.cells), np.nan)
        out[gene_set.name] = (
            pd.Series(np.clip(auc, 0, 1), index=universe.cells),
            pd.Series(ucell, index=universe.cells),
            present,
        )
    return out


def marker_info(gene_set: GeneSet, matrix_genes: Iterable[str]) -> dict[str, object]:
    genes = set(matrix_genes)
    present = [g for g in gene_set.genes if g in genes]
    missing = [g for g in gene_set.genes if g not in genes]
    return {
        "markers_requested": ";".join(gene_set.genes),
        "markers_present": ";".join(present),
        "markers_missing": ";".join(missing),
        "n_markers_requested": len(gene_set.genes),
        "n_markers_present": len(present),
        "marker_coverage_fraction": len(present) / len(gene_set.genes) if gene_set.genes else np.nan,
        "present_list": present,
    }


def infer_state_scope(state_name: str, tissues: Iterable[str], cell_types: Iterable[str]) -> tuple[str, str, str]:
    norm_state = norm_name(state_name)
    tissue_match = ""
    for tissue in sorted({norm_name(t) for t in tissues}, key=len, reverse=True):
        if norm_state.startswith(tissue + "_"):
            tissue_match = tissue
            remainder = norm_state[len(tissue) + 1 :]
            break
    else:
        parts = norm_state.split("_", 1)
        tissue_match = parts[0]
        remainder = parts[1] if len(parts) > 1 else ""

    cell_match = ""
    state_suffix = remainder
    for cell_type in sorted({norm_name(c) for c in cell_types}, key=len, reverse=True):
        if remainder.startswith(cell_type + "_"):
            cell_match = cell_type
            state_suffix = remainder[len(cell_type) + 1 :]
            break
    return tissue_match, cell_match, state_suffix


def state_kind(state_name: str) -> str:
    name = norm_name(state_name)
    if any(token in name for token in PROCESS_TOKENS):
        return "process"
    return "biological"


def is_composite_state(state_name: str) -> bool:
    name = norm_name(state_name)
    return any(token in name for token in COMPOSITE_TOKENS)


def infer_state_class(state_name: str, overrides: dict[str, str]) -> str:
    if state_name in overrides:
        return overrides[state_name]
    name = norm_name(state_name)
    if any(token in name for token in QC_TOKENS):
        return "qc_or_contamination"
    if any(token in name for token in COMPOSITE_TOKENS):
        return "composite_required"
    if any(token in name for token in RARE_PROCESS_TOKENS):
        return "rare_process"
    if any(token in name for token in PROCESS_TOKENS):
        return "process_gradient"
    if any(token in name for token in FUNCTION_TOKENS):
        return "broad_function_gradient"
    if any(token in name for token in IDENTITY_TOKENS):
        return "broad_identity_gradient"
    return "unknown"


def gene_set_scope(gene_set: GeneSet, tissues: Iterable[str], cell_types: Iterable[str]) -> tuple[str, str, str]:
    if gene_set.meta.get("tissue") or gene_set.meta.get("cell_type"):
        tissue = norm_name(gene_set.meta.get("tissue", ""))
        cell_type = norm_name(gene_set.meta.get("cell_type", ""))
        return tissue, cell_type, norm_name(gene_set.meta.get("state_label", gene_set.name))
    return infer_state_scope(gene_set.name, tissues, cell_types)


def gene_set_class(gene_set: GeneSet, overrides: dict[str, str]) -> str:
    if gene_set.meta.get("state_class"):
        return str(gene_set.meta["state_class"])
    return infer_state_class(gene_set.name, overrides)


def score_biological_states(
    matrix: pd.DataFrame,
    ucell_ranks: pd.DataFrame,
    aucell_ranks: pd.DataFrame,
    rank_universe: SparseRankUniverse | None,
    metadata: pd.DataFrame,
    gene_sets: list[GeneSet],
    args: argparse.Namespace,
    sparse_scores: dict[str, tuple[pd.Series, pd.Series, list[str]]] | None = None,
) -> pd.DataFrame:
    rows = []
    tissues = metadata[args.tissue_col].unique()
    cell_types = metadata[args.cell_type_col].unique()
    group_keys = [args.map_id_col, args.tissue_col, args.cell_type_col]
    sparse_scores = sparse_scores if sparse_scores is not None else (
        sparse_rank_scores_for_gene_sets(rank_universe, gene_sets, args.aucell_max_rank, args.max_rank)
        if rank_universe is not None
        else {}
    )
    for gene_set in gene_sets:
        scope_tissue, scope_cell_type, _ = gene_set_scope(gene_set, tissues, cell_types)
        matrix_genes = rank_universe.genes if rank_universe is not None else matrix.columns
        info = marker_info(gene_set, matrix_genes)
        if rank_universe is not None:
            aucell_score, ucell_score, present = sparse_scores[gene_set.name]
            missing = [g for g in gene_set.genes if g not in set(present)]
            info["present_list"] = present
            info["markers_present"] = ";".join(present)
            info["markers_missing"] = ";".join(missing)
            info["n_markers_present"] = len(present)
            info["marker_coverage_fraction"] = len(present) / len(gene_set.genes) if gene_set.genes else np.nan
        else:
            ucell_score = ucell_score_from_ranks(ucell_ranks, info["present_list"], args.max_rank)
            aucell_score = aucell_score_from_ranks(aucell_ranks, info["present_list"], args.aucell_max_rank)
        score_frame = metadata[["cell_id", args.map_id_col, args.tissue_col, args.cell_type_col]].copy()
        score_frame["state_name"] = gene_set.name
        score_frame["ucell_score"] = ucell_score.reindex(score_frame["cell_id"]).to_numpy()
        score_frame["aucell_score"] = aucell_score.reindex(score_frame["cell_id"]).to_numpy()
        score_frame["scope_tissue"] = scope_tissue
        score_frame["scope_cell_type"] = scope_cell_type
        score_frame["state_kind"] = state_kind(gene_set.name)
        score_frame["is_composite_state"] = gene_set.meta.get("is_composite_required", "false") == "true" or is_composite_state(gene_set.name)
        score_frame["state_class"] = gene_set_class(gene_set, args.state_class_overrides)
        score_frame["manifest_present"] = gene_set.meta.get("manifest_present", "false")
        for key, value in info.items():
            if key != "present_list":
                score_frame[key] = value

        relevant = (metadata[args.tissue_col].map(norm_name) == scope_tissue)
        if scope_cell_type:
            relevant &= metadata[args.cell_type_col].map(norm_name) == scope_cell_type
        score_frame = score_frame.loc[relevant.to_numpy()].copy()
        if score_frame.empty:
            continue
        score_frame["score_percentile_within_calibration_group"] = (
            score_frame.groupby(group_keys + ["state_name"])["aucell_score"].rank(pct=True, method="average")
        )
        rows.append(score_frame)
    if not rows:
        raise SystemExit("No biological state scores were relevant to the provided metadata")
    out = pd.concat(rows, ignore_index=True)
    return out.rename(
        columns={
            args.map_id_col: "map_id",
            args.tissue_col: "tissue",
            args.cell_type_col: "annotated_cell_type",
        }
    )


def score_qc_signatures(
    matrix: pd.DataFrame,
    ucell_ranks: pd.DataFrame,
    aucell_ranks: pd.DataFrame,
    rank_universe: SparseRankUniverse | None,
    metadata: pd.DataFrame,
    qc_sets: list[GeneSet],
    args: argparse.Namespace,
    sparse_scores: dict[str, tuple[pd.Series, pd.Series, list[str]]] | None = None,
) -> pd.DataFrame:
    rows = []
    sparse_scores = sparse_scores if sparse_scores is not None else (
        sparse_rank_scores_for_gene_sets(rank_universe, qc_sets, args.aucell_max_rank, args.max_rank)
        if rank_universe is not None
        else {}
    )
    for gene_set in qc_sets:
        matrix_genes = rank_universe.genes if rank_universe is not None else matrix.columns
        info = marker_info(gene_set, matrix_genes)
        if rank_universe is not None:
            aucell_score, ucell_score, present = sparse_scores[gene_set.name]
            missing = [g for g in gene_set.genes if g not in set(present)]
            info["present_list"] = present
            info["markers_present"] = ";".join(present)
            info["markers_missing"] = ";".join(missing)
            info["n_markers_present"] = len(present)
            info["marker_coverage_fraction"] = len(present) / len(gene_set.genes) if gene_set.genes else np.nan
        else:
            ucell_score = ucell_score_from_ranks(ucell_ranks, info["present_list"], args.max_rank)
            aucell_score = aucell_score_from_ranks(aucell_ranks, info["present_list"], args.aucell_max_rank)
        frame = metadata[["cell_id", args.map_id_col, args.tissue_col, args.cell_type_col, args.sample_col]].copy()
        frame["qc_signature_name"] = gene_set.name
        frame["ucell_score"] = ucell_score.reindex(frame["cell_id"]).to_numpy()
        frame["aucell_score"] = aucell_score.reindex(frame["cell_id"]).to_numpy()
        frame["state_class"] = "qc_or_contamination"
        frame["qc_tier"] = gene_set.meta.get("tier", "")
        frame["qc_category"] = gene_set.meta.get("category", "")
        frame["score_percentile_within_sample"] = (
            frame.groupby([args.map_id_col, args.sample_col, "qc_signature_name"])["ucell_score"].rank(pct=True, method="average")
        )
        for key, value in info.items():
            if key != "present_list":
                frame[key] = value
        rows.append(frame)
    out = pd.concat(rows, ignore_index=True)
    return out.rename(
        columns={
            args.map_id_col: "map_id",
            args.tissue_col: "tissue",
            args.cell_type_col: "annotated_cell_type",
        }
    )


def mad(series: pd.Series) -> float:
    med = series.median()
    return float((series - med).abs().median())


def adaptive_high_flag(values: pd.Series) -> pd.Series:
    med = values.median()
    spread = mad(values)
    return values > med + 3 * spread if spread > 0 else pd.Series(False, index=values.index)


def adaptive_low_flag(values: pd.Series) -> pd.Series:
    med = values.median()
    spread = mad(values)
    return values < med - 3 * spread if spread > 0 else pd.Series(False, index=values.index)


def sparse_gene_sum(universe: SparseRankUniverse, mask: np.ndarray) -> pd.Series:
    if not mask.any():
        return pd.Series(0.0, index=universe.cells)
    return pd.Series(np.asarray(universe.matrix[:, mask].sum(axis=1)).ravel(), index=universe.cells)


def qc_metric_inputs(
    matrix: pd.DataFrame,
    raw_counts: SparseRankUniverse | None,
) -> tuple[pd.Series, pd.Series, pd.Series, pd.Series, pd.Series, pd.Series, str]:
    if raw_counts is not None:
        genes = pd.Series(raw_counts.genes, index=raw_counts.genes)
        total_counts = pd.Series(np.asarray(raw_counts.matrix.sum(axis=1)).ravel(), index=raw_counts.cells)
        n_genes = pd.Series(np.diff(raw_counts.matrix.indptr), index=raw_counts.cells)
        mito = sparse_gene_sum(raw_counts, genes.str.startswith("MT-").to_numpy())
        ribo = sparse_gene_sum(raw_counts, genes.str.match(r"RP[LS]").to_numpy())
        hemo = sparse_gene_sum(raw_counts, genes.str.startswith(HEMOGLOBIN_PREFIXES).to_numpy())
        malat = sparse_gene_sum(raw_counts, (genes == "MALAT1").to_numpy())
        return total_counts, n_genes, mito, ribo, hemo, malat, "full_sparse_counts"
    genes = pd.Series(matrix.columns, index=matrix.columns)
    total_counts = matrix.sum(axis=1)
    n_genes = (matrix > 0).sum(axis=1)
    mito = matrix.loc[:, genes.str.startswith("MT-").to_numpy()].sum(axis=1) if genes.str.startswith("MT-").any() else pd.Series(0.0, index=matrix.index)
    ribo = matrix.loc[:, genes.str.match(r"RP[LS]").to_numpy()].sum(axis=1) if genes.str.match(r"RP[LS]").any() else pd.Series(0.0, index=matrix.index)
    hemo_mask = genes.str.startswith(HEMOGLOBIN_PREFIXES)
    hemo = matrix.loc[:, hemo_mask.to_numpy()].sum(axis=1) if hemo_mask.any() else pd.Series(0.0, index=matrix.index)
    malat = matrix["MALAT1"] if "MALAT1" in matrix.columns else pd.Series(0.0, index=matrix.index)
    return total_counts, n_genes, mito, ribo, hemo, malat, "legacy_expression_matrix"


def qc_metrics(matrix: pd.DataFrame, metadata: pd.DataFrame, qc_scores: pd.DataFrame, args: argparse.Namespace, raw_counts: SparseRankUniverse | None = None) -> pd.DataFrame:
    total_counts, n_genes, mito, ribo, hemo, malat, qc_metric_source = qc_metric_inputs(matrix, raw_counts)
    out = metadata[["cell_id", args.map_id_col, args.tissue_col, args.cell_type_col, args.sample_col]].copy()
    out["total_counts"] = total_counts.reindex(out["cell_id"]).to_numpy()
    out["n_genes_detected"] = n_genes.reindex(out["cell_id"]).to_numpy()
    denom = out["total_counts"].replace(0, np.nan)
    out["percent_mitochondrial"] = 100 * np.asarray(mito.reindex(out["cell_id"]) if hasattr(mito, "reindex") else 0) / denom
    out["percent_ribosomal"] = 100 * np.asarray(ribo.reindex(out["cell_id"]) if hasattr(ribo, "reindex") else 0) / denom
    out["percent_hemoglobin"] = 100 * np.asarray(hemo.reindex(out["cell_id"]) if hasattr(hemo, "reindex") else 0) / denom
    out["percent_malat1"] = 100 * np.asarray(malat.reindex(out["cell_id"]) if hasattr(malat, "reindex") else 0) / denom
    out["qc_metric_source"] = qc_metric_source
    if args.doublet_col and args.doublet_col in metadata.columns:
        out["doublet_score"] = pd.to_numeric(metadata[args.doublet_col], errors="coerce")
    else:
        out["doublet_score"] = np.nan

    out["hard_exclusion_flag"] = False
    out["review_flag"] = False
    reasons = {cell: [] for cell in out["cell_id"]}
    for _, group in out.groupby([args.map_id_col, args.sample_col], sort=False):
        idx = group.index
        high_total = adaptive_high_flag(group["total_counts"])
        low_total = adaptive_low_flag(group["total_counts"])
        high_genes = adaptive_high_flag(group["n_genes_detected"])
        low_genes = adaptive_low_flag(group["n_genes_detected"])
        high_mito = adaptive_high_flag(group["percent_mitochondrial"].fillna(0))
        high_ribo = adaptive_high_flag(group["percent_ribosomal"].fillna(0))
        high_hemo = adaptive_high_flag(group["percent_hemoglobin"].fillna(0))
        out.loc[idx[high_mito | high_hemo], "hard_exclusion_flag"] = True
        out.loc[idx[high_total | low_total | high_genes | low_genes | high_ribo], "review_flag"] = True
        for flag, label in [
            (high_mito, "high_percent_mitochondrial"),
            (high_hemo, "high_percent_hemoglobin"),
            (high_ribo, "high_percent_ribosomal"),
            (high_total, "high_total_counts"),
            (low_total, "low_total_counts"),
            (high_genes, "high_n_genes_detected"),
            (low_genes, "low_n_genes_detected"),
        ]:
            for cell in group.loc[flag, "cell_id"]:
                reasons[cell].append(label)

    high_qc = qc_scores.loc[qc_scores["score_percentile_within_sample"] >= args.qc_extreme_percentile].copy()
    hard_qc = high_qc["qc_tier"].str.contains("hard_exclude", na=False) & bool(getattr(args, "allow_qc_signature_hard_exclusion", False))
    review_qc = high_qc["qc_tier"].str.contains("review", na=False)
    hard_cells = set(high_qc.loc[hard_qc, "cell_id"])
    review_cells = set(high_qc.loc[review_qc, "cell_id"])
    out.loc[out["cell_id"].isin(hard_cells), "hard_exclusion_flag"] = True
    out.loc[out["cell_id"].isin(review_cells), "review_flag"] = True
    for _, row in high_qc.loc[hard_qc | review_qc].iterrows():
        reasons[row["cell_id"]].append(row["qc_signature_name"])

    identity = high_qc.loc[high_qc["qc_signature_name"].str.contains("identity", case=False, na=False)]
    ambient = high_qc.loc[high_qc["qc_signature_name"].str.contains("ambient", case=False, na=False)]
    top_identity = identity.sort_values("ucell_score").groupby("cell_id").tail(1).set_index("cell_id")["qc_signature_name"]
    top_ambient = ambient.sort_values("ucell_score").groupby("cell_id").tail(1).set_index("cell_id")["qc_signature_name"]
    out["bad_cell_reason"] = out["cell_id"].map(lambda c: ";".join(sorted(set(reasons[c]))) if reasons[c] else "none")
    out["top_offtarget_identity_signature"] = out["cell_id"].map(top_identity).fillna("none")
    out["top_ambient_signature"] = out["cell_id"].map(top_ambient).fillna("none")
    out["parent_identity_score"] = np.nan
    return out.rename(
        columns={
            args.map_id_col: "map_id",
            args.tissue_col: "tissue",
            args.cell_type_col: "annotated_cell_type",
        }
    )


def explore_aucell_threshold(scores: pd.Series, args: argparse.Namespace, max_active_fraction: float) -> tuple[float, str, str, float]:
    values = pd.to_numeric(scores, errors="coerce").dropna().clip(0, 1)
    if len(values) < args.min_calibration_cells:
        return np.nan, "insufficient_cells", "calibration_group_below_min_cells", np.nan
    if values.quantile(0.75) - values.quantile(0.25) < args.min_score_iqr:
        return np.nan, "continuous_only", "aucell_iqr_below_minimum", np.nan
    if values.nunique() < 4:
        return np.nan, "continuous_only", "too_few_unique_aucell_scores", np.nan
    edges = np.linspace(0, 1, args.aucell_threshold_bins + 1)
    hist, _ = np.histogram(values, bins=edges, density=False)
    smooth = gaussian_filter1d(hist.astype(float), sigma=args.aucell_threshold_smoothing)
    centers = (edges[:-1] + edges[1:]) / 2
    median = float(values.median())
    left = np.where(centers <= median)[0]
    right = np.where(centers > median)[0]
    if len(left) == 0 or len(right) == 0 or smooth[left].max() <= 0 or smooth[right].max() <= 0:
        return np.nan, "continuous_only", "no_two_density_regions", np.nan
    left_peak = left[np.argmax(smooth[left])]
    right_peak = right[np.argmax(smooth[right])]
    if right_peak <= left_peak + 1:
        return np.nan, "continuous_only", "no_density_valley_between_peaks", np.nan
    valley_region = np.arange(left_peak + 1, right_peak)
    valley = valley_region[np.argmin(smooth[valley_region])]
    threshold = float(centers[valley])
    active_fraction = float((values >= threshold).mean())
    if active_fraction < args.aucell_min_active_fraction or active_fraction > max_active_fraction:
        return np.nan, "continuous_only", "active_fraction_outside_bounds", active_fraction
    if values.quantile(0.99) <= threshold:
        return np.nan, "continuous_only", "q99_not_above_threshold", active_fraction
    return threshold, "hard_callable", "minimum_density_threshold", active_fraction


def class_allows_hard_call(state_class: str) -> bool:
    return state_class in {"process_gradient", "rare_process"}


def class_max_active_fraction(state_class: str, args: argparse.Namespace) -> float:
    if state_class == "process_gradient":
        return args.aucell_process_gradient_max_active_fraction
    if state_class == "rare_process":
        return args.aucell_rare_process_max_active_fraction
    return 0.0


def class_continuous_status(state_class: str) -> tuple[str, str]:
    if state_class == "composite_required":
        return "composite_required", "composite_state_requires_additional_logic"
    if state_class in {"broad_identity_gradient", "broad_function_gradient"}:
        return "continuous_only", f"{state_class}_defaults_to_gradient"
    if state_class == "qc_or_contamination":
        return "continuous_only", "qc_signature_not_biological_hard_call"
    if state_class == "unknown":
        return "continuous_only_unknown_class", "unknown_state_class_requires_manifest_or_yaml_for_hard_calls"
    return "continuous_only", "no_supported_aucell_hard_threshold"


def calibrate_thresholds(
    scores: pd.DataFrame,
    args: argparse.Namespace,
    state_thresholds: dict[str, float],
) -> pd.DataFrame:
    rows = []
    for keys, group in scores.groupby(["map_id", "tissue", "annotated_cell_type", "state_name"], sort=False):
        map_id, tissue, cell_type, state_name = keys
        aucell = pd.to_numeric(group["aucell_score"], errors="coerce")
        score_iqr = aucell.quantile(0.75) - aucell.quantile(0.25)
        q50 = float(aucell.quantile(0.50)) if aucell.notna().any() else np.nan
        q75 = float(aucell.quantile(0.75)) if aucell.notna().any() else np.nan
        q90 = float(aucell.quantile(0.90)) if aucell.notna().any() else np.nan
        q95 = float(aucell.quantile(0.95)) if aucell.notna().any() else np.nan
        q99 = float(aucell.quantile(0.99)) if aucell.notna().any() else np.nan
        base = {
            "map_id": map_id,
            "tissue": tissue,
            "annotated_cell_type": cell_type,
            "state_name": state_name,
            "mixture_threshold": np.nan,
            "n_cells_in_calibration_group": len(group),
            "score_iqr": score_iqr,
            "q50_score": q50,
            "q75_score": q75,
            "q90_score_diagnostic": q90,
            "q95_score_diagnostic": q95,
            "q99_score": q99,
            "n_markers_requested": group["n_markers_requested"].iloc[0],
            "n_markers_present": group["n_markers_present"].iloc[0],
            "marker_coverage_fraction": group["marker_coverage_fraction"].iloc[0],
            "state_class": group["state_class"].iloc[0] if "state_class" in group.columns else "unknown",
        }
        if len(group) < args.min_calibration_cells:
            rows.append({**base, "threshold_method": "insufficient_cells", "threshold_value": np.nan, "threshold_status": "insufficient_cells", "threshold_reason": "calibration_group_below_min_cells", "active_fraction_at_threshold": np.nan})
            continue
        if score_iqr < args.min_score_iqr:
            rows.append({**base, "threshold_method": "continuous_only", "threshold_value": np.nan, "threshold_status": "continuous_only", "threshold_reason": "aucell_iqr_below_minimum", "active_fraction_at_threshold": np.nan})
            continue
        if not class_allows_hard_call(base["state_class"]) and state_name not in state_thresholds:
            status, reason = class_continuous_status(base["state_class"])
            rows.append({**base, "threshold_method": status, "threshold_value": np.nan, "threshold_status": status, "threshold_reason": reason, "active_fraction_at_threshold": np.nan})
            continue
        if state_name in state_thresholds:
            threshold = state_thresholds[state_name]
            rows.append(
                {
                    **base,
                    "threshold_method": "yaml_aucell_threshold",
                    "threshold_value": threshold,
                    "threshold_status": "hard_callable",
                    "threshold_reason": "yaml_override",
                    "active_fraction_at_threshold": float((aucell >= threshold).mean()),
                }
            )
            continue
        threshold, status, reason, active_fraction = explore_aucell_threshold(aucell, args, class_max_active_fraction(base["state_class"], args))
        rows.append(
            {
                **base,
                "threshold_method": "aucell_minimum_density" if status == "hard_callable" else "continuous_only",
                "threshold_value": threshold,
                "threshold_status": status,
                "threshold_reason": reason,
                "active_fraction_at_threshold": active_fraction,
            }
        )
    return pd.DataFrame(rows)


def compute_activity(
    scores: pd.DataFrame,
    state_col: str,
    state_type: str,
    thresholds: pd.DataFrame | None = None,
) -> pd.DataFrame:
    rows = []
    score_state_col = "state_name" if state_col == "state_name" else "qc_signature_name"
    threshold_lookup = None
    if thresholds is not None and not thresholds.empty:
        threshold_lookup = thresholds.set_index(["map_id", "tissue", "annotated_cell_type", "state_name"])
    if {"map_id", "tissue", "annotated_cell_type"}.issubset(scores.columns):
        group_cols = ["map_id", "tissue", "annotated_cell_type", score_state_col]
    else:
        group_cols = [score_state_col]
    for group_key, group in scores.groupby(group_cols, sort=False):
        state_name = group_key[-1] if isinstance(group_key, tuple) else group_key
        out = group[["cell_id"]].copy()
        if {"map_id", "tissue", "annotated_cell_type"}.issubset(group.columns):
            out[["map_id", "tissue", "annotated_cell_type"]] = group[["map_id", "tissue", "annotated_cell_type"]]
        out["state_type"] = state_type
        out["state_name"] = state_name
        out["ucell_score"] = group["ucell_score"].to_numpy()
        out["aucell_score"] = group["aucell_score"].to_numpy()
        out["threshold_value"] = np.nan
        out["threshold_status"] = "not_applicable_qc" if state_type == "qc" else "continuous_only"
        out["threshold_method"] = "not_applicable_qc" if state_type == "qc" else "continuous_only"
        out["threshold_reason"] = "not_applicable_qc" if state_type == "qc" else "no_threshold_metadata"
        out["q50_score"] = np.nan
        out["q75_score"] = np.nan
        out["q90_score_diagnostic"] = np.nan
        out["q95_score_diagnostic"] = np.nan
        out["q99_score"] = np.nan
        out["aucell_percentile_within_group"] = 0.0
        out["state_activity_weight_gradient"] = 0.0
        out["state_activity_weight_hightail"] = 0.0
        out["state_class"] = group["state_class"].iloc[0] if "state_class" in group.columns else ("qc_or_contamination" if state_type == "qc" else "unknown")
        if state_type == "biological" and threshold_lookup is not None:
            key_cols = ["map_id", "tissue", "annotated_cell_type", "state_name"]
            merged = group[key_cols].drop_duplicates().merge(
                thresholds[
                    key_cols
                    + [
                        "threshold_value",
                        "threshold_status",
                        "threshold_method",
                        "threshold_reason",
                        "q75_score",
                        "q50_score",
                        "q90_score_diagnostic",
                        "q95_score_diagnostic",
                        "q99_score",
                        "state_class",
                    ]
                ],
                on=key_cols,
                how="left",
            )
            for _, threshold_row in merged.iterrows():
                row_mask = (
                    group["map_id"].eq(threshold_row["map_id"])
                    & group["tissue"].eq(threshold_row["tissue"])
                    & group["annotated_cell_type"].eq(threshold_row["annotated_cell_type"])
                    & group["state_name"].eq(threshold_row["state_name"])
                )
                idx = group.index[row_mask]
                if len(idx) == 0:
                    continue
                out_idx = out.index.intersection(idx)
                threshold_value = threshold_row["threshold_value"]
                threshold_status = threshold_row["threshold_status"]
                threshold_method = threshold_row["threshold_method"]
                threshold_reason = threshold_row["threshold_reason"]
                q50 = threshold_row["q50_score"]
                q75 = threshold_row["q75_score"]
                q90 = threshold_row["q90_score_diagnostic"]
                q95 = threshold_row["q95_score_diagnostic"]
                q99 = threshold_row["q99_score"]
                out.loc[out_idx, "threshold_value"] = threshold_value
                out.loc[out_idx, "threshold_status"] = threshold_status
                out.loc[out_idx, "threshold_method"] = threshold_method
                out.loc[out_idx, "threshold_reason"] = threshold_reason
                out.loc[out_idx, "q50_score"] = q50
                out.loc[out_idx, "q75_score"] = q75
                out.loc[out_idx, "q90_score_diagnostic"] = q90
                out.loc[out_idx, "q95_score_diagnostic"] = q95
                out.loc[out_idx, "q99_score"] = q99
                scores_for_weight = pd.to_numeric(out.loc[out_idx, "aucell_score"], errors="coerce")
                percentile = group.loc[out_idx, "score_percentile_within_calibration_group"].fillna(0.0).to_numpy(dtype=float)
                out.loc[out_idx, "aucell_percentile_within_group"] = percentile
                out.loc[out_idx, "state_activity_weight_gradient"] = np.square(percentile)
                out.loc[out_idx, "state_activity_weight_hightail"] = np.clip((percentile - 0.90) / 0.10, 0, 1)
                out.loc[out_idx, "state_class"] = threshold_row["state_class"]
                if threshold_status == "hard_callable" and pd.notna(threshold_value):
                    denom = q99 - threshold_value if pd.notna(q99) and pd.notna(threshold_value) else np.nan
                    method = "aucell_threshold_to_q99"
                    baseline = threshold_value
                else:
                    denom = q99 - q75 if pd.notna(q99) and pd.notna(q75) else np.nan
                    method = "aucell_q75_to_q99_continuous_only"
                    baseline = q75
                if pd.notna(denom) and denom > 0 and pd.notna(baseline):
                    out.loc[out_idx, "state_activity_weight"] = ((scores_for_weight - baseline) / denom).clip(0, 1).fillna(0.0).to_numpy()
                else:
                    out.loc[out_idx, "state_activity_weight"] = 0.0
                out.loc[out_idx, "soft_weight_method"] = method
            if merged.empty:
                out["state_activity_weight"] = 0.0
                out["state_activity_weight_gradient"] = 0.0
                out["state_activity_weight_hightail"] = 0.0
                out["aucell_percentile_within_group"] = 0.0
                out["state_class"] = "unknown"
                out["threshold_reason"] = "missing_threshold_metadata"
                out["q50_score"] = np.nan
                out["soft_weight_method"] = "missing_threshold_metadata"
        else:
            percentile_col = "score_percentile_within_sample" if "score_percentile_within_sample" in group.columns else None
            percentile = group[percentile_col].fillna(0.0).to_numpy(dtype=float) if percentile_col else np.zeros(len(group))
            out["aucell_percentile_within_group"] = percentile
            out["state_activity_weight_gradient"] = np.square(percentile)
            out["state_activity_weight_hightail"] = np.clip((percentile - 0.90) / 0.10, 0, 1)
            out["state_activity_weight"] = out["state_activity_weight_gradient"]
            out["state_class"] = group["state_class"].iloc[0] if "state_class" in group.columns else "qc_or_contamination"
            out["threshold_reason"] = "not_applicable_qc" if state_type == "qc" else "no_threshold_metadata"
            out["q50_score"] = np.nan
            out["soft_weight_method"] = "aucell_percentile_gradient_qc_score" if state_type == "qc" else "aucell_percentile_gradient"
        out["marker_coverage_fraction"] = group["marker_coverage_fraction"].iloc[0]
        out["n_markers_present"] = group["n_markers_present"].iloc[0]
        out["n_markers_total"] = group["n_markers_requested"].iloc[0]
        rows.append(out)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def build_ucell_scores_output(bio_scores: pd.DataFrame, qc_scores: pd.DataFrame) -> pd.DataFrame:
    bio = bio_scores.rename(columns={"state_name": "state_name", "n_markers_requested": "n_markers_total"}).copy()
    bio["state_type"] = "biological"
    bio_out = bio[
        [
            "cell_id",
            "state_type",
            "state_name",
            "ucell_score",
            "aucell_score",
            "n_markers_present",
            "n_markers_total",
            "marker_coverage_fraction",
            "markers_present",
            "markers_missing",
        ]
    ]
    qc = qc_scores.rename(columns={"qc_signature_name": "state_name", "n_markers_requested": "n_markers_total"}).copy()
    qc["state_type"] = "qc"
    if "markers_present" not in qc.columns:
        qc["markers_present"] = ""
    if "markers_missing" not in qc.columns:
        qc["markers_missing"] = ""
    qc_out = qc[
        [
            "cell_id",
            "state_type",
            "state_name",
            "ucell_score",
            "aucell_score",
            "n_markers_present",
            "n_markers_total",
            "marker_coverage_fraction",
            "markers_present",
            "markers_missing",
        ]
    ]
    return pd.concat([bio_out, qc_out], ignore_index=True)


def confidence_from_margin(margin: float, exploratory: bool = False) -> str:
    if exploratory:
        return "exploratory"
    if pd.isna(margin):
        return "low"
    if margin >= 0.05:
        return "high"
    if margin >= 0:
        return "medium"
    return "low"


def hard_assignments_from_activity(
    activity: pd.DataFrame,
    qc_thresholds: dict[str, float],
    args: argparse.Namespace,
) -> pd.DataFrame:
    frame = activity.copy()
    qc_threshold = frame["state_name"].map(qc_thresholds)
    frame["threshold"] = np.where(
        frame["state_type"].eq("qc"),
        qc_threshold.fillna(args.default_qc_threshold),
        frame["threshold_value"],
    )
    frame["threshold_source"] = np.where(
        frame["state_type"].eq("qc"),
        np.where(qc_threshold.notna(), "yaml", "default_qc"),
        frame["threshold_method"],
    )
    frame["marker_coverage_pass"] = (frame["n_markers_present"] >= args.min_markers_present) & (
        frame["marker_coverage_fraction"] >= args.min_marker_coverage
    )
    callable_or_qc = frame["state_type"].eq("qc") | frame["threshold_status"].eq("hard_callable")
    score_for_call = np.where(frame["state_type"].eq("qc"), frame["state_activity_weight"], frame["aucell_score"])
    frame["hard_call"] = frame["marker_coverage_pass"] & callable_or_qc & (score_for_call >= frame["threshold"])
    frame["reason"] = np.select(
        [
            ~frame["marker_coverage_pass"],
            ~callable_or_qc,
            frame["hard_call"],
        ],
        [
            "insufficient_marker_coverage",
            "continuous_only_no_hard_threshold",
            "activity_above_threshold",
        ],
        default="activity_below_threshold",
    )
    for col in ["map_id", "tissue", "annotated_cell_type"]:
        if col not in frame.columns:
            frame[col] = ""
    return frame[
        [
            "cell_id",
            "map_id",
            "tissue",
            "annotated_cell_type",
            "state_type",
            "state_name",
            "aucell_score",
            "ucell_score",
            "state_activity_weight",
            "threshold",
            "hard_call",
            "threshold_source",
            "marker_coverage_pass",
            "reason",
        ]
    ].copy()


def qc_exclusions(activity: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    qc = activity.loc[activity["state_type"] == "qc"].copy()
    cells = pd.DataFrame({"cell_id": sorted(activity["cell_id"].unique())})
    if args.exclude_qc_above is None:
        cells["excluded"] = False
        cells["exclusion_reason"] = "qc_exclusion_not_requested"
        cells["triggering_qc_states"] = "none"
        cells["max_qc_activity"] = cells["cell_id"].map(qc.groupby("cell_id")["state_activity_weight"].max()).fillna(0.0)
        return cells
    selected = set(args.exclude_qc_states.split(",")) if args.exclude_qc_states else set(qc["state_name"].unique())
    qc = qc.loc[qc["state_name"].isin(selected)].copy()
    triggers = qc.loc[qc["state_activity_weight"] >= args.exclude_qc_above]
    grouped = triggers.groupby("cell_id")["state_name"].agg(lambda x: ";".join(sorted(set(x))))
    cells["triggering_qc_states"] = cells["cell_id"].map(grouped).fillna("none")
    cells["excluded"] = cells["triggering_qc_states"].ne("none")
    cells["exclusion_reason"] = np.where(cells["excluded"], "qc_activity_above_threshold", "not_excluded")
    cells["max_qc_activity"] = cells["cell_id"].map(qc.groupby("cell_id")["state_activity_weight"].max()).fillna(0.0)
    return cells


def qc_signature_review_flags(activity: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    qc = activity.loc[activity["state_type"] == "qc"].copy()
    if qc.empty:
        return empty_table(["cell_id", "state_name", "review_flag", "qc_activity_weight", "reason"])
    qc["review_flag"] = qc["state_activity_weight_gradient"] >= args.default_qc_threshold
    qc["reason"] = np.where(qc["review_flag"], "qc_signature_activity_above_review_threshold", "below_review_threshold")
    return qc[["cell_id", "state_name", "review_flag", "state_activity_weight_gradient", "reason"]].rename(
        columns={"state_activity_weight_gradient": "qc_activity_weight"}
    )


def qc_legacy_fixed_tail_flags(activity: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    qc = activity.loc[activity["state_type"] == "qc"].copy()
    if qc.empty:
        return empty_table(["cell_id", "state_name", "legacy_tail_flag", "aucell_percentile_within_group", "reason"])
    qc["legacy_tail_flag"] = qc["aucell_percentile_within_group"] >= args.qc_extreme_percentile
    qc["reason"] = np.where(qc["legacy_tail_flag"], "legacy_fixed_tail_percentile_flag", "below_legacy_tail")
    return qc[["cell_id", "state_name", "legacy_tail_flag", "aucell_percentile_within_group", "reason"]]


def state_activity_lookup(activity: pd.DataFrame) -> dict[tuple[str, str], float]:
    bio = activity.loc[activity["state_type"] == "biological"]
    return {(row.cell_id, row.state_name): float(row.state_activity_weight) for row in bio.itertuples(index=False)}


def loo_activity_for_gene(
    gene: str,
    state_name: str,
    gene_set: GeneSet,
    matrix: pd.DataFrame,
    aucell_ranks: pd.DataFrame,
    rank_universe: SparseRankUniverse | None,
    args: argparse.Namespace,
    threshold_status: str,
    threshold_value: float,
    q75: float,
    q99: float,
) -> pd.Series | None:
    if gene not in gene_set.genes:
        return None
    loo_genes = [g for g in gene_set.genes if g != gene]
    present = [g for g in loo_genes if g in matrix.columns]
    if rank_universe is not None:
        auc_score, _, present_sparse = sparse_rank_scores_for_gene_set(rank_universe, loo_genes, args.aucell_max_rank, args.max_rank)
        score = auc_score
        if not present_sparse:
            return pd.Series(0.0, index=matrix.index)
    else:
        if not present:
            return pd.Series(0.0, index=matrix.index)
        score = aucell_score_from_ranks(aucell_ranks, present, args.aucell_max_rank)
    if score.empty:
        return pd.Series(0.0, index=matrix.index)
    if threshold_status == "hard_callable" and pd.notna(threshold_value):
        baseline = threshold_value
    else:
        baseline = q75
    denom = q99 - baseline if pd.notna(q99) and pd.notna(baseline) else np.nan
    if pd.isna(denom) or denom <= 0:
        return pd.Series(0.0, index=matrix.index)
    return ((score - baseline) / denom).clip(0, 1).fillna(0.0)


EXPECTED_EXPRESSION_COLUMNS = [
    "gene",
    "state_name",
    "state_type",
    "weighted_mean_expression",
    "weighted_detection_fraction",
    "whole_cell_type_mean_expression",
    "log2_weighted_vs_all",
    "n_cells",
    "sum_state_weight",
    "leave_one_gene_out_used",
]


def expected_expression_summary(
    matrix: pd.DataFrame,
    activity: pd.DataFrame,
    biological_sets: dict[str, GeneSet],
    aucell_ranks: pd.DataFrame,
    rank_universe: SparseRankUniverse | None,
    args: argparse.Namespace,
    query_genes: list[str],
) -> pd.DataFrame:
    rows = []
    if args.mode == "hard":
        return empty_table(EXPECTED_EXPRESSION_COLUMNS)
    all_mean = matrix[query_genes].mean(axis=0)
    bio_activity = activity.loc[activity["state_type"] == "biological"].copy()
    for state_name, group in bio_activity.groupby("state_name", sort=False):
        weights = group.set_index("cell_id")["state_activity_weight"].reindex(matrix.index).fillna(0.0)
        threshold_status = group["threshold_status"].iloc[0]
        threshold_value = group["threshold_value"].iloc[0]
        q75 = group["q75_score"].iloc[0]
        q99 = group["q99_score"].iloc[0]
        for gene in query_genes:
            used_loo = False
            gene_weights = weights
            if args.leave_one_gene_out and state_name in biological_sets:
                loo = loo_activity_for_gene(gene, state_name, biological_sets[state_name], matrix, aucell_ranks, rank_universe, args, threshold_status, threshold_value, q75, q99)
                if loo is not None:
                    gene_weights = loo.reindex(matrix.index).fillna(0.0)
                    used_loo = True
            denom = gene_weights.sum()
            if denom > 0:
                weighted_mean = float((gene_weights * matrix[gene]).sum() / denom)
                weighted_detection = float((gene_weights * (matrix[gene] > 0)).sum() / denom)
            else:
                weighted_mean = np.nan
                weighted_detection = np.nan
            rows.append(
                {
                    "gene": gene,
                    "state_name": state_name,
                    "state_type": "biological",
                    "weighted_mean_expression": weighted_mean,
                    "weighted_detection_fraction": weighted_detection,
                    "whole_cell_type_mean_expression": float(all_mean[gene]),
                    "log2_weighted_vs_all": float(np.log2((weighted_mean + 1e-6) / (all_mean[gene] + 1e-6))) if pd.notna(weighted_mean) else np.nan,
                    "n_cells": int(len(matrix)),
                    "sum_state_weight": float(denom),
                    "leave_one_gene_out_used": used_loo,
                }
            )
    return pd.DataFrame(rows, columns=EXPECTED_EXPRESSION_COLUMNS)


HARD_EXPRESSION_COLUMNS = [
    "gene",
    "state_name",
    "state_type",
    "mean_expression_state_positive",
    "mean_expression_state_negative",
    "detection_fraction_positive",
    "detection_fraction_negative",
    "log2_fc_positive_vs_negative",
    "p_value",
    "q_value",
    "n_positive_cells",
    "n_negative_cells",
]


def hard_expression_summary(matrix: pd.DataFrame, hard: pd.DataFrame, args: argparse.Namespace, query_genes: list[str]) -> pd.DataFrame:
    rows = []
    if args.mode == "expected":
        return empty_table(HARD_EXPRESSION_COLUMNS)
    bio_hard = hard.loc[hard["state_type"] == "biological"].copy()
    for state_name, group in bio_hard.groupby("state_name", sort=False):
        calls = group.set_index("cell_id")["hard_call"].reindex(matrix.index).fillna(False).astype(bool)
        pos = matrix.loc[calls]
        neg = matrix.loc[~calls]
        for gene in query_genes:
            if len(pos) > 0 and len(neg) > 0:
                stat = stats.mannwhitneyu(pos[gene], neg[gene], alternative="two-sided").pvalue
            else:
                stat = np.nan
            rows.append(
                {
                    "gene": gene,
                    "state_name": state_name,
                    "state_type": "biological",
                    "mean_expression_state_positive": float(pos[gene].mean()) if len(pos) else np.nan,
                    "mean_expression_state_negative": float(neg[gene].mean()) if len(neg) else np.nan,
                    "detection_fraction_positive": float((pos[gene] > 0).mean()) if len(pos) else np.nan,
                    "detection_fraction_negative": float((neg[gene] > 0).mean()) if len(neg) else np.nan,
                    "log2_fc_positive_vs_negative": float(np.log2((pos[gene].mean() + 1e-6) / (neg[gene].mean() + 1e-6))) if len(pos) and len(neg) else np.nan,
                    "p_value": stat,
                    "n_positive_cells": int(len(pos)),
                    "n_negative_cells": int(len(neg)),
                }
            )
    out = pd.DataFrame(rows, columns=[c for c in HARD_EXPRESSION_COLUMNS if c != "q_value"])
    out["q_value"] = bh_fdr(out["p_value"]) if not out.empty else []
    return out[HARD_EXPRESSION_COLUMNS]


def donor_weighted_expression(matrix: pd.DataFrame, metadata: pd.DataFrame, weights: pd.Series, donor_col: str, query_genes: list[str]) -> pd.DataFrame:
    rows = []
    meta = metadata.set_index("cell_id").reindex(matrix.index)
    for donor, idx in meta.groupby(donor_col).groups.items():
        donor_cells = list(idx)
        w = weights.reindex(donor_cells).fillna(0.0)
        denom = w.sum()
        values = {"donor_id": donor, "sum_state_weight": float(denom)}
        if denom > 0:
            weighted = matrix.loc[donor_cells, query_genes].multiply(w, axis=0).sum(axis=0) / denom
        else:
            weighted = pd.Series(np.nan, index=query_genes)
        values.update(weighted.to_dict())
        rows.append(values)
    return pd.DataFrame(rows)


def fit_simple_lm(y: pd.Series, design: pd.DataFrame, phenotype: str) -> tuple[float, float, int]:
    frame = pd.concat([y.rename("y"), design], axis=1).dropna()
    if len(frame) < 3 or phenotype not in frame.columns or frame[phenotype].nunique() < 2:
        return np.nan, np.nan, len(frame)
    x = frame.drop(columns=["y"])
    x = pd.get_dummies(x, drop_first=True, dtype=float)
    x.insert(0, "intercept", 1.0)
    if phenotype not in x.columns:
        phen_cols = [c for c in x.columns if c.startswith(phenotype + "_")]
        if not phen_cols:
            return np.nan, np.nan, len(frame)
        coef_col = phen_cols[0]
    else:
        coef_col = phenotype
    X = x.to_numpy(float)
    Y = frame["y"].to_numpy(float)
    coef, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
    resid = Y - X @ coef
    dof = len(Y) - X.shape[1]
    if dof <= 0:
        return coef[list(x.columns).index(coef_col)], np.nan, len(frame)
    sigma2 = float((resid @ resid) / dof)
    cov = sigma2 * np.linalg.pinv(X.T @ X)
    idx = list(x.columns).index(coef_col)
    se = float(np.sqrt(cov[idx, idx])) if cov[idx, idx] >= 0 else np.nan
    p = 2 * stats.t.sf(abs(coef[idx] / se), dof) if se and np.isfinite(se) and se > 0 else np.nan
    return float(coef[idx]), float(p), len(frame)


def de_summaries(
    matrix: pd.DataFrame,
    metadata: pd.DataFrame,
    activity: pd.DataFrame,
    hard: pd.DataFrame,
    phenotypes: pd.DataFrame | None,
    args: argparse.Namespace,
    query_genes: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    expected_cols = ["gene", "state_name", "phenotype", "coefficient", "coefficient_units", "p_value", "q_global", "q_by_trait", "q_by_gene", "n_donors", "sum_state_weight", "model_formula"]
    hard_cols = ["gene", "state_name", "phenotype", "coefficient", "coefficient_units", "p_value", "q_global", "q_by_trait", "q_by_gene", "n_donors", "n_state_positive_cells", "model_formula"]
    if args.mode == "expected":
        hard_requested = False
        expected_requested = True
    elif args.mode == "hard":
        hard_requested = True
        expected_requested = False
    else:
        hard_requested = True
        expected_requested = True
    if phenotypes is None or phenotypes.empty:
        return empty_table(expected_cols), empty_table(hard_cols)
    if args.donor_col not in phenotypes.columns:
        raise SystemExit(f"Phenotype table is missing donor column: {args.donor_col}")
    pheno_cols = [c for c in phenotypes.columns if c != args.donor_col]
    if not pheno_cols:
        return empty_table(expected_cols), empty_table(hard_cols)
    phenotype = pheno_cols[0]
    donor_design = phenotypes.set_index(args.donor_col)[pheno_cols]
    exp_rows = []
    hard_rows = []
    bio_activity = activity.loc[activity["state_type"] == "biological"]
    for state_name, group in bio_activity.groupby("state_name", sort=False):
        weights = group.set_index("cell_id")["state_activity_weight"].reindex(matrix.index).fillna(0.0)
        if expected_requested:
            donor_expr = donor_weighted_expression(matrix, metadata, weights, args.donor_col, query_genes).set_index("donor_id")
            design = donor_design.reindex(donor_expr.index)
            for gene in query_genes:
                coef, p, n = fit_simple_lm(donor_expr[gene], design, phenotype)
                exp_rows.append(
                    {
                        "gene": gene,
                        "state_name": state_name,
                        "phenotype": phenotype,
                        "coefficient": coef,
                        "coefficient_units": "expression_per_phenotype_unit",
                        "p_value": p,
                        "n_donors": n,
                        "sum_state_weight": float(weights.sum()),
                        "model_formula": f"E_d[{gene}|{state_name}] ~ {phenotype}",
                    }
                )
        if hard_requested:
            calls = hard.loc[(hard["state_type"] == "biological") & (hard["state_name"] == state_name)].set_index("cell_id")["hard_call"].reindex(matrix.index).fillna(False).astype(bool)
            meta = metadata.set_index("cell_id").reindex(matrix.index)
            for gene in query_genes:
                donor_values = []
                for donor, idx in meta.groupby(args.donor_col).groups.items():
                    donor_cells = list(idx)
                    active_cells = [c for c in donor_cells if calls.loc[c]]
                    value = float(matrix.loc[active_cells, gene].mean()) if active_cells else np.nan
                    donor_values.append({"donor_id": donor, "value": value})
                donor_frame = pd.DataFrame(donor_values).set_index("donor_id")
                design = donor_design.reindex(donor_frame.index)
                coef, p, n = fit_simple_lm(donor_frame["value"], design, phenotype)
                hard_rows.append(
                    {
                        "gene": gene,
                        "state_name": state_name,
                        "phenotype": phenotype,
                        "coefficient": coef,
                        "coefficient_units": "expression_per_phenotype_unit",
                        "p_value": p,
                        "n_donors": n,
                        "n_state_positive_cells": int(calls.sum()),
                        "model_formula": f"mean_hard_positive[{gene}|{state_name}] ~ {phenotype}",
                    }
                )
    exp = pd.DataFrame(exp_rows, columns=[c for c in expected_cols if not c.startswith("q_")])
    hard_de = pd.DataFrame(hard_rows, columns=[c for c in hard_cols if not c.startswith("q_")])
    for frame in [exp, hard_de]:
        if frame.empty:
            continue
        frame["q_global"] = bh_fdr(frame["p_value"])
        frame["q_by_trait"] = frame.groupby("phenotype", group_keys=False)["p_value"].apply(bh_fdr)
        frame["q_by_gene"] = frame.groupby("gene", group_keys=False)["p_value"].apply(bh_fdr)
    return exp.reindex(columns=expected_cols), hard_de.reindex(columns=hard_cols)


def call_states(scores: pd.DataFrame, thresholds: pd.DataFrame, bad_flags: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    threshold_cols = [
        "map_id",
        "tissue",
        "annotated_cell_type",
        "state_name",
        "threshold_method",
        "threshold_value",
        "threshold_status",
        "score_iqr",
        "n_cells_in_calibration_group",
        "threshold_reason",
    ]
    frame = scores.merge(thresholds[threshold_cols], on=["map_id", "tissue", "annotated_cell_type", "state_name"], how="left")
    frame = frame.merge(bad_flags[["cell_id", "hard_exclusion_flag", "review_flag"]], on="cell_id", how="left")
    calls = []
    for _, row in frame.iterrows():
        requires_composite = bool(row["is_composite_state"])
        composite_rule = "none"
        if row["n_markers_present"] < args.min_markers_present or row["marker_coverage_fraction"] < args.min_marker_coverage:
            call = "insufficient_marker_coverage"
            reason = "insufficient_marker_coverage"
        elif row["threshold_status"] == "insufficient_cells":
            call = "not_called_insufficient_cells"
            reason = "calibration_group_below_min_cells"
        elif row["threshold_status"] == "composite_required":
            call = "composite_required"
            reason = "composite_state_requires_additional_logic"
        elif row["threshold_status"] == "continuous_only":
            call = "continuous_only"
            reason = row.get("threshold_reason", "no_supported_aucell_hard_threshold")
        elif pd.isna(row["threshold_value"]):
            call = "not_called_missing_threshold"
            reason = "missing_threshold"
        elif row["aucell_score"] >= row["threshold_value"]:
            if bool(row["hard_exclusion_flag"]):
                call = "ambiguous_qc_flagged"
                reason = "score_high_but_hard_exclusion_qc_flag"
            elif requires_composite:
                call = "exploratory_marker_high"
                reason = "composite_state_requires_validation"
            else:
                call = "active"
                reason = "score_above_threshold"
        else:
            call = "inactive"
            reason = "score_below_threshold"
        margin = row["aucell_score"] - row["threshold_value"] if not pd.isna(row["threshold_value"]) else np.nan
        calls.append(
            {
                "cell_id": row["cell_id"],
                "map_id": row["map_id"],
                "tissue": row["tissue"],
                "annotated_cell_type": row["annotated_cell_type"],
                "state_name": row["state_name"],
                "ucell_score": row["ucell_score"],
                "aucell_score": row["aucell_score"],
                "threshold_value": row["threshold_value"],
                "call": call,
                "confidence": confidence_from_margin(margin, call == "exploratory_marker_high"),
                "reason": reason,
                "hard_exclusion_flag": bool(row["hard_exclusion_flag"]),
                "review_flag": bool(row["review_flag"]),
                "composite_rule_used": composite_rule,
                "requires_composite_validation": requires_composite and composite_rule == "none",
                "state_kind": row["state_kind"],
                "state_class": row["state_class"],
            }
        )
    return pd.DataFrame(calls)


def multilabel_summary(calls: pd.DataFrame, bad_flags: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for cell_id, group in calls.groupby("cell_id", sort=False):
        active_bio = sorted(group.loc[(group["call"] == "active") & (group["state_kind"] == "biological"), "state_name"].unique())
        process = sorted(group.loc[(group["call"] == "active") & (group["state_kind"] == "process"), "state_name"].unique())
        exploratory = sorted(group.loc[group["call"] == "exploratory_marker_high", "state_name"].unique())
        flags = bad_flags.loc[bad_flags["cell_id"] == cell_id].iloc[0]
        use_program = "exclude_hard_flagged" if flags["hard_exclusion_flag"] else "include"
        use_de = "exclude_hard_flagged" if flags["hard_exclusion_flag"] else "include"
        rows.append(
            {
                "cell_id": cell_id,
                "map_id": flags["map_id"],
                "tissue": flags["tissue"],
                "annotated_cell_type": flags["annotated_cell_type"],
                "active_biological_states": ";".join(active_bio) if active_bio else "none",
                "active_process_flags": ";".join(process) if process else "none",
                "exploratory_states": ";".join(exploratory) if exploratory else "none",
                "qc_flags": flags["bad_cell_reason"],
                "recommended_use_for_program_inference": use_program,
                "recommended_use_for_de": use_de,
                "summary_label": "; ".join(active_bio + process + exploratory) if (active_bio or process or exploratory) else "no_active_state_calls",
            }
        )
    return pd.DataFrame(rows)


def state_call_summary(calls: pd.DataFrame, metadata: pd.DataFrame, thresholds: pd.DataFrame, args: argparse.Namespace) -> pd.DataFrame:
    active = calls["call"].isin(["active", "exploratory_marker_high"])
    frame = calls.copy()
    frame["is_active"] = active
    frame = frame.merge(metadata[["cell_id", args.donor_col]], on="cell_id", how="left")
    rows = []
    for keys, group in frame.groupby(["map_id", "tissue", "annotated_cell_type", "state_name"], sort=False):
        threshold = thresholds.set_index(["map_id", "tissue", "annotated_cell_type", "state_name"]).loc[keys]
        donor_active = group.groupby(args.donor_col)["is_active"].sum()
        rows.append(
            {
                "map_id": keys[0],
                "tissue": keys[1],
                "annotated_cell_type": keys[2],
                "state_name": keys[3],
                "n_cells_scored": len(group),
                "n_cells_active": int(group["is_active"].sum()),
                "active_fraction": float(group["is_active"].mean()) if len(group) else np.nan,
                "n_donors_with_active_cells": int((donor_active > 0).sum()),
                "median_active_cells_per_donor": float(donor_active.median()) if len(donor_active) else np.nan,
                "threshold_method": threshold["threshold_method"],
                "threshold_value": threshold["threshold_value"],
                "marker_coverage_fraction": threshold["marker_coverage_fraction"],
                "confidence_summary": ";".join(sorted(group["confidence"].astype(str).unique())),
            }
        )
    return pd.DataFrame(rows)


def acceptance_checks(summary: pd.DataFrame, scores: pd.DataFrame, qc_scores: pd.DataFrame, calls: pd.DataFrame) -> list[str]:
    failures = []
    for _, group in summary.groupby(["map_id", "tissue", "annotated_cell_type"], sort=False):
        if len(group) > 1 and group["active_fraction"].nunique(dropna=True) == 1:
            failures.append(
                "All biological states have the same active fraction for "
                f"{group.iloc[0]['map_id']}/{group.iloc[0]['tissue']}/{group.iloc[0]['annotated_cell_type']}"
            )
    if scores["marker_coverage_fraction"].isna().any():
        failures.append("Marker coverage is missing for at least one biological score row")
    if qc_scores.empty:
        failures.append("Auxiliary bad-cell QC signatures were not scored")
    if not {"hard_exclusion_flag", "review_flag"}.issubset(set(calls.columns)):
        failures.append("Hard-exclusion and review flags are not both present in calls")
    composite = calls.loc[calls["requires_composite_validation"]]
    composite_called = composite.loc[composite["call"].isin(["active", "ambiguous_qc_flagged"])]
    if not composite_called.empty:
        failures.append("Composite states without rules were not labeled exploratory marker-high")
    return failures


def write_methods(out_dir: Path, args: argparse.Namespace, mapping_info: dict[str, str], failures: list[str]) -> None:
    lines = [
        "# CMDKP cell-state scoring method",
        "",
        f"- Input matrix type: `{args.expression_kind}`",
        f"- Rank universe: `{args.rank_10x_dir or args.rank_matrix_mtx or 'expression_matrix_fallback'}`",
        "- UCell implementation: `local_ucell_style_rank_statistic`",
        "- UCell package version: `not_used`",
        f"- UCell parameters: `maxRank={args.max_rank}`, `ties.method=average`, `chunk.size=not_applicable_python_runner`, `missing_genes=skip`, `knn_smoothing=none`",
        "- AUCell implementation: `local_aucell_style_rank_statistic`",
        f"- AUCell parameters: `aucMaxRank={args.aucell_max_rank}`, `ties.method=average`, threshold exploration bins `{args.aucell_threshold_bins}`",
        f"- Biological GMT: `{args.biological_gmt}`",
        f"- Auxiliary bad-cell QC GMT: `{args.qc_gmt}`",
        f"- State manifest: `{args.state_manifest or 'not_supplied_state_name_fallback'}`",
        f"- Gene ID handling: `{mapping_info['gene_id_type']}`",
        f"- Duplicate collapse method: `{mapping_info['duplicate_collapse_method']}`",
        f"- Marker coverage rules: score with at least 1 marker present; confident calls require `n_markers_present >= {args.min_markers_present}` and `marker_coverage_fraction >= {args.min_marker_coverage}`",
        "- Activity method: AUCell threshold-to-q99 state-excess weights for hard-callable biological states; AUCell q75-to-q99 weights for continuous-only biological states",
        "- Q90 and Q95 are distribution diagnostics only; they are not used as default hard-call thresholds.",
        f"- Thresholding method: AUCell minimum-density threshold exploration with YAML overrides; default QC threshold `{args.default_qc_threshold}`",
        f"- Threshold eligibility group: `map_id + tissue + annotated_cell_type + state_name`; minimum cells `{args.min_calibration_cells}`; minimum score IQR `{args.min_score_iqr}`",
        f"- Summary query genes: `{args.query_genes or args.query_gene or 'all_expression_matrix_genes'}`",
        f"- Summary mode: `{args.mode}`",
        "- Bad-cell hard exclusions: technical metric outliers and QC signatures marked with hard-exclude tiers in the auxiliary GMT.",
        "- Bad-cell review flags: ribosomal/translation, heat-shock, immediate-early, and other review-tier QC signatures when extreme.",
        "- Scores are primary and calls are secondary. State calls are multi-label and not mutually exclusive.",
        "",
        "## Acceptance checks",
    ]
    if failures:
        lines.extend(f"- FAILED: {failure}" for failure in failures)
    else:
        lines.append("- All implemented acceptance checks passed.")
    (out_dir / "state_scoring_method.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--expression", default="", help="Long or wide expression TSV/TSV.GZ")
    parser.add_argument("--expression-matrix", default="", help="Alias for --expression")
    parser.add_argument("--rank-10x-dir", default="", help="10x directory for full sparse rank-universe AUCell scoring")
    parser.add_argument("--rank-matrix-mtx", default="", help="Matrix Market sparse matrix for full rank-universe AUCell scoring")
    parser.add_argument("--rank-genes", default="", help="Gene/features TSV for --rank-matrix-mtx")
    parser.add_argument("--rank-cells", default="", help="Cell/barcode TSV for --rank-matrix-mtx")
    parser.add_argument("--metadata", default="", help="Cell metadata TSV/TSV.GZ")
    parser.add_argument("--cell-metadata", default="", help="Alias for --metadata")
    parser.add_argument("--biological-gmt", default="")
    parser.add_argument("--states-gmt", default="", help="Alias for --biological-gmt")
    parser.add_argument("--qc-gmt", default="")
    parser.add_argument("--qc-states-gmt", default="", help="Alias for --qc-gmt")
    parser.add_argument("--state-manifest", default="", help="Optional TSV with state_name, tissue, cell_type, state_class, and related metadata")
    parser.add_argument("--require-state-manifest", action="store_true", help="Fail unless every GMT state has production manifest metadata")
    parser.add_argument("--qc-raw-10x-dir", default="", help="Optional full raw 10x counts for direct QC metrics; defaults to --rank-10x-dir when available")
    parser.add_argument("--out-dir", default="")
    parser.add_argument("--expression-kind", choices=["raw_counts", "log1p_normalized"], default="log1p_normalized")
    parser.add_argument("--phenotypes", default="")
    parser.add_argument("--state-thresholds-yaml", default="")
    parser.add_argument("--state-threshold-config", default="", help="Alias for --state-thresholds-yaml")
    parser.add_argument("--state-class-config", default="", help="YAML mapping state names to state_class")
    parser.add_argument("--qc-thresholds-yaml", default="")
    parser.add_argument("--mode", choices=["expected", "hard", "both"], default="both")
    parser.add_argument("--exclude-qc-above", type=float, default=None)
    parser.add_argument("--exclude-qc-states", default="")
    parser.add_argument("--gene-map", default="")
    parser.add_argument("--duplicate-collapse", choices=["auto", "sum", "mean", "max"], default="auto")
    parser.add_argument("--map-id-col", default="map_id")
    parser.add_argument("--tissue-col", default="tissue")
    parser.add_argument("--cell-type-col", default="cell_type")
    parser.add_argument("--donor-col", default="donor_id")
    parser.add_argument("--sample-col", default="sample_id")
    parser.add_argument("--doublet-col", default="")
    parser.add_argument("--max-rank", type=int, default=1500)
    parser.add_argument("--aucell-max-rank", type=int, default=0, help="Top-ranked genes used for AUCell-style scoring; default 5%% of genes")
    parser.add_argument("--aucell-threshold-bins", type=int, default=64)
    parser.add_argument("--aucell-threshold-smoothing", type=float, default=1.0)
    parser.add_argument("--aucell-min-active-fraction", type=float, default=0.01)
    parser.add_argument("--aucell-max-active-fraction", type=float, default=0.30, help="Compatibility alias for process-gradient max active fraction")
    parser.add_argument("--aucell-process-gradient-max-active-fraction", type=float, default=0.30)
    parser.add_argument("--aucell-rare-process-max-active-fraction", type=float, default=0.10)
    parser.add_argument("--min-markers-present", type=int, default=5)
    parser.add_argument("--min-marker-coverage", type=float, default=0.5)
    parser.add_argument("--default-qc-threshold", type=float, default=0.95)
    parser.add_argument("--leave-one-gene-out", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--query-genes", default="", help="Optional newline-delimited gene list restricting expression and DE summaries")
    parser.add_argument("--query-gene", default="", help="Optional comma-separated genes restricting expression and DE summaries")
    parser.add_argument("--min-calibration-cells", type=int, default=100)
    parser.add_argument("--min-score-iqr", type=float, default=0.01)
    parser.add_argument("--qc-extreme-percentile", type=float, default=0.99)
    parser.add_argument("--allow-qc-signature-hard-exclusion", action="store_true", help="Allow QC signature percentile tails to set hard_exclusion_flag; otherwise signature tails are review diagnostics only")
    parser.add_argument("--legacy-selected-gene-summaries", choices=["skip", "write"], default="skip", help="Write compatibility selected-gene expression/de summaries from this runner")
    parser.add_argument("--min-rank-genes", type=int, default=5000)
    parser.add_argument("--allow-small-rank-universe", action="store_true")
    parser.add_argument("--parent-cell-filter", default="", help="Metadata filter expression column=value1,value2 applied before scoring")
    parser.add_argument("--ucell-scores-out", default="")
    parser.add_argument("--aucell-state-activity-out", default="")
    parser.add_argument("--cell-state-activity-out", default="")
    parser.add_argument("--cell-state-hard-assignments-out", default="")
    parser.add_argument("--qc-exclusions-out", default="")
    parser.add_argument("--expression-expected-assignments-out", default="")
    parser.add_argument("--expression-hard-assignments-out", default="")
    parser.add_argument("--de-expected-assignments-out", default="")
    parser.add_argument("--de-hard-assignments-out", default="")
    parser.add_argument("--state-summary-out", default="")
    parser.add_argument("--run-summary-out", default="")
    parser.add_argument("--timing-log-out", default="")
    parser.add_argument("--allow-acceptance-failures", action="store_true")
    args = parser.parse_args()

    timer = StepTimer()

    args.expression = args.expression or args.expression_matrix
    args.metadata = args.metadata or args.cell_metadata
    args.biological_gmt = args.biological_gmt or args.states_gmt
    args.state_thresholds_yaml = args.state_thresholds_yaml or args.state_threshold_config
    args.qc_gmt = args.qc_gmt or args.qc_states_gmt or str(Path(__file__).resolve().parents[1] / "dat" / "qc" / "cmdkp_all_tissues_minimal_bad_cell_qc_signatures.gmt")
    if not args.expression or not args.metadata or not args.biological_gmt:
        raise SystemExit("--expression-matrix, --cell-metadata, and --states-gmt are required")
    if not args.out_dir:
        if not all(
            [
                args.ucell_scores_out,
                args.aucell_state_activity_out,
                args.cell_state_activity_out,
                args.cell_state_hard_assignments_out,
                args.qc_exclusions_out,
                args.expression_expected_assignments_out,
                args.expression_hard_assignments_out,
                args.de_expected_assignments_out,
                args.de_hard_assignments_out,
                args.state_summary_out,
                args.run_summary_out,
            ]
        ):
            raise SystemExit("Provide --out-dir or all explicit output paths")
        args.out_dir = str(Path(args.ucell_scores_out).parent)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    metadata = read_table(args.metadata)
    timer.mark("load_metadata", n_rows=len(metadata))
    if args.cell_type_col not in metadata.columns and "annotated_cell_type" in metadata.columns:
        args.cell_type_col = "annotated_cell_type"
    required_meta = ["cell_id", args.map_id_col, args.tissue_col, args.cell_type_col, args.donor_col]
    if args.sample_col not in metadata.columns:
        metadata[args.sample_col] = metadata[args.donor_col]
    missing = [c for c in required_meta if c not in metadata.columns]
    if missing:
        raise SystemExit(f"Metadata is missing required column(s): {', '.join(missing)}")
    metadata = metadata.drop_duplicates("cell_id").copy()
    metadata = apply_metadata_filter(metadata, args.parent_cell_filter)
    args.state_class_overrides = load_state_class_config(args.state_class_config)
    if args.aucell_max_active_fraction != 0.30:
        args.aucell_process_gradient_max_active_fraction = args.aucell_max_active_fraction
    state_manifest = load_state_manifest(args.state_manifest)
    expression, mapping_info = harmonize_expression(read_expression_input(args.expression), args.expression_kind, args.gene_map, args.duplicate_collapse)
    matrix = expression_matrix(expression, metadata["cell_id"])
    timer.mark("load_expression_matrix", n_cells=matrix.shape[0], n_genes=matrix.shape[1])
    if args.aucell_max_rank <= 0:
        rank_gene_count = len(read_10x_features(resolve_10x_path(args.rank_10x_dir, ["features.tsv.gz", "features.tsv", "genes.tsv.gz", "genes.tsv"]))) if args.rank_10x_dir and not args.rank_genes else None
        if rank_gene_count is None and args.rank_genes:
            rank_gene_count = len(read_10x_features(args.rank_genes))
        args.aucell_max_rank = max(1, int(np.ceil((rank_gene_count or matrix.shape[1]) * 0.05)))
    query_genes = resolve_query_genes(args, matrix)
    rank_universe = load_sparse_rank_universe(args, metadata)
    timer.mark(
        "load_rank_universe",
        n_cells=rank_universe.matrix.shape[0] if rank_universe is not None else None,
        n_genes=rank_universe.matrix.shape[1] if rank_universe is not None else None,
    )
    n_genes_ranked = int(rank_universe.matrix.shape[1]) if rank_universe is not None else int(matrix.shape[1])
    if n_genes_ranked < args.min_rank_genes:
        message = (
            f"Rank universe has {n_genes_ranked} genes, below --min-rank-genes {args.min_rank_genes}. "
            "AUCell should be run on a broad/full rank universe."
        )
        if args.allow_small_rank_universe:
            print("Warning: " + message, file=sys.stderr)
        else:
            raise SystemExit(message + " Use --allow-small-rank-universe for tests or exploratory runs.")
    if rank_universe is None:
        print(
            "Warning: no sparse rank universe was supplied; AUCell/UCell scoring will use "
            "the expression matrix. For production runs, pass --rank-10x-dir or "
            "--rank-matrix-mtx --rank-genes --rank-cells.",
            file=sys.stderr,
        )
    ucell_ranks = rank_matrix_for_ucell(matrix, args.max_rank) if rank_universe is None else pd.DataFrame(index=matrix.index)
    aucell_ranks = rank_matrix_for_aucell(matrix) if rank_universe is None else pd.DataFrame(index=matrix.index)

    raw_biological_sets = read_gmt(args.biological_gmt)
    validate_state_manifest(state_manifest, raw_biological_sets, args.require_state_manifest)
    biological_sets = apply_state_manifest(raw_biological_sets, state_manifest)
    qc_sets = apply_state_manifest(read_gmt(args.qc_gmt), state_manifest) if args.qc_gmt and Path(args.qc_gmt).exists() else []
    combined_sparse_scores = None
    if rank_universe is not None:
        combined_sparse_scores = sparse_rank_scores_for_gene_sets(
            rank_universe,
            biological_sets + qc_sets,
            args.aucell_max_rank,
            args.max_rank,
        )
        timer.mark("score_sparse_rank_universe_once", n_gene_sets=len(biological_sets) + len(qc_sets))
    biological_scores = score_biological_states(
        matrix,
        ucell_ranks,
        aucell_ranks,
        rank_universe,
        metadata,
        biological_sets,
        args,
        sparse_scores=combined_sparse_scores,
    )
    qc_scores = (
        score_qc_signatures(matrix, ucell_ranks, aucell_ranks, rank_universe, metadata, qc_sets, args, sparse_scores=combined_sparse_scores)
        if qc_sets
        else empty_table(["cell_id", "map_id", "tissue", "annotated_cell_type", "qc_signature_name", "ucell_score", "aucell_score", "state_class", "qc_tier", "qc_category", "score_percentile_within_sample", "n_markers_requested", "n_markers_present", "marker_coverage_fraction", "markers_present", "markers_missing"])
    )
    timer.mark("assemble_score_tables", n_biological_score_rows=len(biological_scores), n_qc_score_rows=len(qc_scores))
    state_thresholds = load_threshold_yaml(args.state_thresholds_yaml)
    qc_thresholds = load_threshold_yaml(args.qc_thresholds_yaml)
    thresholds = calibrate_thresholds(biological_scores, args, state_thresholds)
    timer.mark("calibrate_thresholds", n_threshold_rows=len(thresholds))
    qc_raw_counts = load_sparse_10x_dir(args.qc_raw_10x_dir, metadata) if args.qc_raw_10x_dir else rank_universe
    bad_flags = qc_metrics(matrix, metadata, qc_scores, args, raw_counts=qc_raw_counts)
    bio_activity = compute_activity(biological_scores, "state_name", "biological", thresholds)
    qc_activity = compute_activity(qc_scores, "qc_signature_name", "qc")
    activity = pd.concat([bio_activity, qc_activity], ignore_index=True)
    hard = hard_assignments_from_activity(activity, qc_thresholds, args)
    exclusions = qc_exclusions(activity, args)
    qc_review = qc_signature_review_flags(activity, args)
    qc_legacy = qc_legacy_fixed_tail_flags(activity, args)
    timer.mark("compute_activity_and_qc", n_activity_rows=len(activity))
    if args.legacy_selected_gene_summaries == "write":
        expression_expected = expected_expression_summary(matrix, activity, {s.name: s for s in biological_sets}, aucell_ranks, rank_universe, args, query_genes)
        expression_hard = hard_expression_summary(matrix, hard, args, query_genes)
        phenotypes = read_table(args.phenotypes) if args.phenotypes else None
        de_expected, de_hard = de_summaries(matrix, metadata, activity, hard, phenotypes, args, query_genes)
    else:
        expression_expected = empty_table(EXPECTED_EXPRESSION_COLUMNS)
        expression_hard = empty_table(HARD_EXPRESSION_COLUMNS)
        de_expected = empty_table(["gene", "state_name", "phenotype", "coefficient", "coefficient_units", "p_value", "q_global", "q_by_trait", "q_by_gene", "n_donors", "sum_state_weight", "model_formula"])
        de_hard = empty_table(["gene", "state_name", "phenotype", "coefficient", "coefficient_units", "p_value", "q_global", "q_by_trait", "q_by_gene", "n_donors", "n_state_positive_cells", "model_formula"])
    timer.mark("legacy_selected_gene_summaries", mode=args.legacy_selected_gene_summaries, n_expected_rows=len(expression_expected), n_hard_rows=len(expression_hard))

    calls = call_states(biological_scores, thresholds, bad_flags, args)
    multilabel = multilabel_summary(calls, bad_flags)
    summary = state_call_summary(calls, metadata, thresholds, args)
    state_summary_keys = ["map_id", "tissue", "annotated_cell_type", "state_type", "state_name"]
    state_summary = activity.groupby(state_summary_keys, sort=False).agg(
        n_cells_scored=("cell_id", "nunique"),
        mean_ucell_score=("ucell_score", "mean"),
        mean_aucell_score=("aucell_score", "mean"),
        mean_activity_weight=("state_activity_weight", "mean"),
        marker_coverage_fraction=("marker_coverage_fraction", "first"),
    ).reset_index()
    hard_counts = hard.groupby(state_summary_keys, sort=False).agg(
        n_hard_assigned=("hard_call", "sum"),
        threshold=("threshold", "first"),
    ).reset_index()
    state_summary = state_summary.merge(hard_counts, on=state_summary_keys, how="left")
    state_summary["hard_assigned_fraction"] = state_summary["n_hard_assigned"] / state_summary["n_cells_scored"]
    failures = acceptance_checks(summary, biological_scores, qc_scores, calls)
    timer.mark("build_calls_and_summaries", n_failures=len(failures))

    ucell_scores = build_ucell_scores_output(biological_scores, qc_scores)
    aucell_state_activity = activity.loc[activity["state_type"] == "biological"].merge(
        hard.loc[hard["state_type"] == "biological", ["cell_id", "state_name", "hard_call"]],
        on=["cell_id", "state_name"],
        how="left",
    )
    write_table(ucell_scores, output_path(args, "ucell_scores_out", "ucell_scores.tsv.gz"))
    write_table(
        aucell_state_activity[["cell_id", "state_name", "aucell_score", "threshold_status", "hard_call", "state_activity_weight"]],
        output_path(args, "aucell_state_activity_out", "aucell_state_activity.tsv.gz"),
    )
    activity_cols = ["cell_id", "map_id", "tissue", "annotated_cell_type", "state_type", "state_name", "aucell_score", "ucell_score", "aucell_percentile_within_group", "state_activity_weight", "state_activity_weight_gradient", "state_activity_weight_hightail", "state_class", "soft_weight_method", "threshold_value", "threshold_status", "threshold_reason", "threshold_method", "q50_score", "q75_score", "q90_score_diagnostic", "q95_score_diagnostic", "q99_score", "n_markers_present", "n_markers_total", "marker_coverage_fraction"]
    for col in activity_cols:
        if col not in activity.columns:
            activity[col] = ""
    write_table(activity[activity_cols], output_path(args, "cell_state_activity_out", "cell_state_activity.tsv.gz"))
    write_table(hard, output_path(args, "cell_state_hard_assignments_out", "cell_state_hard_assignments.tsv.gz"))
    write_table(exclusions, output_path(args, "qc_exclusions_out", "qc_exclusions.tsv.gz"))
    write_table(exclusions, out_dir / "qc_applied_exclusions.tsv.gz")
    write_table(bad_flags, out_dir / "qc_direct_metric_flags.tsv.gz")
    write_table(qc_review, out_dir / "qc_signature_review_flags.tsv.gz")
    write_table(qc_legacy, out_dir / "qc_legacy_fixed_tail_flags.tsv.gz")
    write_table(expression_expected, output_path(args, "expression_expected_assignments_out", "expression_expected_assignments.tsv.gz"))
    write_table(expression_hard, output_path(args, "expression_hard_assignments_out", "expression_hard_assignments.tsv.gz"))
    write_table(de_expected, output_path(args, "de_expected_assignments_out", "de_expected_assignments.tsv.gz"))
    write_table(de_hard, output_path(args, "de_hard_assignments_out", "de_hard_assignments.tsv.gz"))
    write_table(state_summary[["map_id", "tissue", "annotated_cell_type", "state_type", "state_name", "n_cells_scored", "mean_ucell_score", "mean_aucell_score", "mean_activity_weight", "n_hard_assigned", "hard_assigned_fraction", "marker_coverage_fraction", "threshold"]], output_path(args, "state_summary_out", "state_summary.tsv.gz"))
    write_table(biological_scores.drop(columns=["scope_tissue", "scope_cell_type", "state_kind", "is_composite_state"]), out_dir / "cell_state_scores.tsv.gz")
    write_table(thresholds, out_dir / "cell_state_thresholds.tsv.gz")
    write_table(calls.drop(columns=["state_kind"]), out_dir / "cell_state_calls.tsv.gz")
    write_table(qc_scores.drop(columns=["qc_tier", "qc_category"]), out_dir / "qc_signature_scores.tsv.gz")
    write_table(bad_flags, out_dir / "bad_cell_qc_flags.tsv.gz")
    write_table(multilabel, out_dir / "cell_multilabel_state_summary.tsv.gz")
    write_table(summary, out_dir / "state_call_summary.tsv.gz")
    write_methods(out_dir, args, mapping_info, failures)
    timer.mark("write_outputs")

    run_summary = {
        "input_files": {
            "expression_matrix": args.expression,
            "rank_10x_dir": args.rank_10x_dir or None,
            "rank_matrix_mtx": args.rank_matrix_mtx or None,
            "cell_metadata": args.metadata,
            "states_gmt": args.biological_gmt,
            "qc_states_gmt": args.qc_gmt,
            "state_manifest": args.state_manifest or None,
            "qc_raw_10x_dir": args.qc_raw_10x_dir or args.rank_10x_dir or None,
            "phenotypes": args.phenotypes or None,
        },
        "parameters": {
            "mode": args.mode,
            "max_rank": args.max_rank,
            "aucell_max_rank": args.aucell_max_rank,
            "aucMaxRank": args.aucell_max_rank,
            "n_genes_ranked": n_genes_ranked,
            "min_rank_genes": args.min_rank_genes,
            "query_genes": args.query_genes or None,
            "query_gene": args.query_gene or None,
            "n_query_genes": len(query_genes),
            "default_qc_threshold": args.default_qc_threshold,
            "leave_one_gene_out": args.leave_one_gene_out,
            "exclude_qc_above": args.exclude_qc_above,
            "allow_qc_signature_hard_exclusion": args.allow_qc_signature_hard_exclusion,
            "legacy_selected_gene_summaries": args.legacy_selected_gene_summaries,
            "require_state_manifest": args.require_state_manifest,
        },
        "n_cells": int(matrix.shape[0]),
        "n_genes": int(matrix.shape[1]),
        "n_rank_universe_cells": int(rank_universe.matrix.shape[0]) if rank_universe is not None else None,
        "n_rank_universe_genes": int(rank_universe.matrix.shape[1]) if rank_universe is not None else None,
        "n_states": int(len(biological_sets)),
        "n_qc_states": int(len(qc_sets)),
        "n_excluded_cells": int(exclusions["excluded"].sum()),
        "qc_metric_source": str(bad_flags["qc_metric_source"].iloc[0]) if "qc_metric_source" in bad_flags.columns and not bad_flags.empty else None,
        "timing": timer.rows,
        "software_versions": {"python": sys.version.split()[0], "pandas": pd.__version__, "numpy": np.__version__},
        "timestamp": datetime.now().isoformat(timespec="seconds"),
    }
    output_path(args, "run_summary_out", "run_summary.json").write_text(json.dumps(run_summary, indent=2), encoding="utf-8")
    write_table(pd.DataFrame(timer.rows), output_path(args, "timing_log_out", "timing_log.tsv"))

    if failures and not args.allow_acceptance_failures:
        raise SystemExit("Acceptance checks failed:\n" + "\n".join(failures))


if __name__ == "__main__":
    main()
