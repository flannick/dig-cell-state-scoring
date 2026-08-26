#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import sys
from dataclasses import dataclass
from datetime import datetime

import numpy as np
import pandas as pd
from scipy.io import mmread

VALUE_TYPES = {'raw_counts', 'linear_cp10k', 'log1p_cp10k', 'linear_normalized', 'log1p_normalized'}
min_rank_genes = 5000
max_rank = 1500


@dataclass
class GeneSet:
    name: str
    genes: list


@dataclass
class SparseRankUniverse:
    matrix: object
    cells: list
    genes: list


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'r')


def read_10x_features(path):
    with open_text(path) as handle:
        return [line.rstrip('\n').split('\t')[1] for line in handle]


def read_one_column(path):
    with open_text(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def read_gmt(path):
    sets = []
    with open_text(path) as handle:
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            genes = list(dict.fromkeys(g for g in parts[2:] if g))
            sets.append(GeneSet(name=parts[0], genes=genes))
    return sets


def load_rank_universe(rank_10x_dir, metadata):
    matrix = mmread(f'{rank_10x_dir}/matrix.mtx.gz').tocsr()
    genes = read_10x_features(f'{rank_10x_dir}/features.tsv.gz')
    cells = read_one_column(f'{rank_10x_dir}/barcodes.tsv.gz')
    if matrix.shape == (len(genes), len(cells)):
        matrix = matrix.T.tocsr()
    elif matrix.shape != (len(cells), len(genes)):
        raise SystemExit(f'Rank matrix shape {matrix.shape} does not match cells x genes ({len(cells)}, {len(genes)})')
    cell_to_pos = {cell: i for i, cell in enumerate(cells)}
    missing = [cell for cell in metadata['cell_id'] if cell not in cell_to_pos]
    if missing:
        raise SystemExit(f'Rank universe is missing {len(missing)} metadata cell IDs')
    order = [cell_to_pos[cell] for cell in metadata['cell_id']]
    matrix = matrix[order, :].tocsr()
    return SparseRankUniverse(matrix=matrix, cells=metadata['cell_id'].astype(str).tolist(), genes=genes)


def sparse_rank_scores_for_gene_sets(universe, gene_sets, aucell_max_rank, progress_every_cells):
    gene_to_pos = {gene: i for i, gene in enumerate(universe.genes)}
    set_positions = {}
    pos_to_sets = {}
    for gene_set in gene_sets:
        positions = {gene_to_pos[g] for g in gene_set.genes if g in gene_to_pos}
        set_positions[gene_set.name] = positions
        for position in positions:
            pos_to_sets.setdefault(position, []).append(gene_set.name)

    n_cells = len(universe.cells)
    auc_values = {gene_set.name: np.zeros(n_cells) for gene_set in gene_sets}
    rank_sums = {gene_set.name: np.zeros(n_cells) for gene_set in gene_sets}
    hit_counts = {gene_set.name: np.zeros(n_cells, dtype=np.int32) for gene_set in gene_sets}
    top_rank = max(max_rank, aucell_max_rank)

    if progress_every_cells:
        print(f'[sparse-score] starting {n_cells} cells x {len(gene_sets)} gene sets; top_rank={top_rank}', file=sys.stderr, flush=True)
    for cell_idx in range(n_cells):
        if progress_every_cells and (cell_idx == 0 or cell_idx + 1 == n_cells or (cell_idx + 1) % progress_every_cells == 0):
            print(f'[sparse-score] cells {cell_idx + 1}/{n_cells} ({100.0 * (cell_idx + 1) / n_cells:.1f}%)', file=sys.stderr, flush=True)
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
                if rank <= aucell_max_rank:
                    auc_values[set_name][cell_idx] += aucell_max_rank - rank + 1
                if rank <= max_rank:
                    rank_sums[set_name][cell_idx] += rank
                    hit_counts[set_name][cell_idx] += 1

    scores = {}
    for gene_set in gene_sets:
        n_present = len(set_positions[gene_set.name])
        if n_present == 0:
            scores[gene_set.name] = (pd.Series(np.nan, index=universe.cells), pd.Series(np.nan, index=universe.cells))
            continue
        auc = np.clip(auc_values[gene_set.name] / (n_present * aucell_max_rank), 0, 1)
        max_u = n_present * max_rank - (n_present * (n_present + 1)) / 2
        if max_u > 0:
            rank_sum = rank_sums[gene_set.name] + (n_present - hit_counts[gene_set.name]) * (max_rank + 1)
            u_stat = rank_sum - (n_present * (n_present + 1)) / 2
            ucell = np.clip(1 - (u_stat / max_u), 0, 1)
        else:
            ucell = np.full(n_cells, np.nan)
        scores[gene_set.name] = (pd.Series(auc, index=universe.cells), pd.Series(ucell, index=universe.cells))
    return scores


def score_activity(metadata, gene_sets, scores, metadata_cols, group_cols, score_col, state_type):
    if not gene_sets:
        return pd.DataFrame(columns=metadata_cols + ['cell_id', 'state_name', 'aucell_score', 'ucell_score', 'state_activity_weight_gradient', 'state_activity_weight_hightail', 'state_type'])
    rows = []
    for gene_set in gene_sets:
        aucell, ucell = scores[gene_set.name]
        frame = metadata[['cell_id'] + metadata_cols].copy()
        frame['state_name'] = gene_set.name
        frame['aucell_score'] = aucell.reindex(frame['cell_id']).to_numpy()
        frame['ucell_score'] = ucell.reindex(frame['cell_id']).to_numpy()
        rows.append(frame)
    frame = pd.concat(rows, ignore_index=True)
    percentile = frame.groupby(group_cols + ['state_name'])[score_col].rank(pct=True, method='average').fillna(0)
    frame['state_activity_weight_gradient'] = np.square(percentile)
    frame['state_activity_weight_hightail'] = np.clip((percentile - 0.90) / 0.10, 0, 1)
    frame['state_type'] = state_type
    return frame


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--rank-10x-dir', required=True)
    parser.add_argument('--rank-value-type', required=True, choices=sorted(VALUE_TYPES))
    parser.add_argument('--cell-metadata', required=True)
    parser.add_argument('--states-gmt', required=True)
    parser.add_argument('--state-manifest', required=True)
    parser.add_argument('--require-state-manifest', action='store_true')
    parser.add_argument('--qc-gmt', default='')
    parser.add_argument('--allow-small-rank-universe', action='store_true')
    parser.add_argument('--tissue-col', required=True)
    parser.add_argument('--cell-type-col', required=True)
    parser.add_argument('--donor-col', required=True)
    parser.add_argument('--sample-col', required=True)
    parser.add_argument('--progress-every-cells', type=int, required=True)
    parser.add_argument('--legacy-selected-gene-summaries', choices=['skip'], required=True)
    parser.add_argument('--api-minimal-output', action='store_true')
    parser.add_argument('--out-dir', required=True)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    print(f'[state-scoring] loading metadata: {args.cell_metadata}', file=sys.stderr, flush=True)
    metadata = pd.read_csv(args.cell_metadata, sep='\t', compression='infer')
    print(f'[state-scoring] loaded metadata: {len(metadata)} rows', file=sys.stderr, flush=True)
    required = ['cell_id', args.tissue_col, args.cell_type_col, args.donor_col, args.sample_col]
    missing = [c for c in required if c not in metadata.columns]
    if missing:
        raise SystemExit(f'Metadata is missing required column(s): {", ".join(missing)}')
    metadata = metadata.drop_duplicates('cell_id').copy()

    print('[state-scoring] loading sparse rank universe', file=sys.stderr, flush=True)
    rank_universe = load_rank_universe(args.rank_10x_dir, metadata)
    print(f'[state-scoring] loaded sparse rank universe: {rank_universe.matrix.shape[0]} cells x {rank_universe.matrix.shape[1]} genes; value_type={args.rank_value_type}', file=sys.stderr, flush=True)
    if len(rank_universe.genes) < min_rank_genes:
        message = f'Rank universe has {len(rank_universe.genes)} genes, below {min_rank_genes}.'
        if args.allow_small_rank_universe:
            print('Warning: ' + message, file=sys.stderr)
        else:
            raise SystemExit(message)
    aucell_max_rank = max(1, math.ceil(len(rank_universe.genes) * 0.05))

    print(f'[state-scoring] loading biological GMT: {args.states_gmt}', file=sys.stderr, flush=True)
    biological_sets = read_gmt(args.states_gmt)
    if args.require_state_manifest:
        manifest_states = set(pd.read_csv(args.state_manifest, sep='\t', compression='infer')['state_name'])
        missing_states = [gene_set.name for gene_set in biological_sets if gene_set.name not in manifest_states]
        if missing_states:
            raise SystemExit(f'State manifest is missing {len(missing_states)} GMT state(s): {", ".join(missing_states[:10])}')
    qc_sets = read_gmt(args.qc_gmt) if args.qc_gmt and os.path.exists(args.qc_gmt) else []

    scores = sparse_rank_scores_for_gene_sets(rank_universe, biological_sets + qc_sets, aucell_max_rank, args.progress_every_cells)

    bio = score_activity(metadata, biological_sets, scores, [args.tissue_col, args.cell_type_col], [args.tissue_col, args.cell_type_col], 'aucell_score', 'biological')
    qc = score_activity(metadata, qc_sets, scores, [args.tissue_col, args.cell_type_col, args.sample_col], [args.tissue_col, args.sample_col], 'ucell_score', 'qc')
    activity = pd.concat([bio, qc], ignore_index=True).rename(columns={args.tissue_col: 'tissue', args.cell_type_col: 'cell_type'})
    activity = activity[['cell_id', 'tissue', 'cell_type', 'state_type', 'state_name', 'aucell_score', 'ucell_score', 'state_activity_weight_gradient', 'state_activity_weight_hightail']]
    activity.to_csv(f'{args.out_dir}/cell_state_activity.tsv.gz', sep='\t', index=False, compression='gzip')

    with open(f'{args.out_dir}/run_summary.json', 'w') as f:
        json.dump({
            'rank_10x_dir': args.rank_10x_dir,
            'rank_value_type': args.rank_value_type,
            'cell_metadata': args.cell_metadata,
            'states_gmt': args.states_gmt,
            'qc_gmt': args.qc_gmt or None,
            'n_cells': len(metadata),
            'n_rank_universe_genes': len(rank_universe.genes),
            'n_states': len(biological_sets),
            'n_qc_states': len(qc_sets),
            'aucell_max_rank': aucell_max_rank,
            'timestamp': datetime.now().isoformat(timespec='seconds'),
        }, f, indent=2)


if __name__ == '__main__':
    main()
