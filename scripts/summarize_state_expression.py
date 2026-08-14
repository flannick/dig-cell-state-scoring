#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
import json
import os
from datetime import datetime

import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.io import mmread

from matrix_value_types import VALUE_TYPES, linearize_expression_matrix

pseudocount = 0.05
min_mean_for_log2fc = 0.01
weight_specs = [
    ('gradient_percentile_squared', 'state_activity_weight_gradient'),
    ('high_tail_percentile_90_100', 'state_activity_weight_hightail'),
]


def open_text(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path, 'r')


def read_one_column(path):
    with open_text(path) as handle:
        return [line.strip() for line in handle if line.strip()]


def read_features(path):
    with open_text(path) as handle:
        return [line.rstrip('\n').split('\t')[1] for line in handle]


def load_10x(directory):
    genes = read_features(f'{directory}/features.tsv.gz')
    cells = read_one_column(f'{directory}/barcodes.tsv.gz')
    matrix = mmread(f'{directory}/matrix.mtx.gz').tocsr()
    if matrix.shape == (len(genes), len(cells)):
        matrix = matrix.T.tocsr()
    elif matrix.shape != (len(cells), len(genes)):
        raise SystemExit(f'10x matrix shape {matrix.shape} does not match feature/cell files')
    return matrix.astype(float).tocsr(), cells, genes


def sparse_square(matrix):
    out = matrix.copy()
    out.data = np.square(out.data)
    return out


def weighted_vs_parent_p_values(state_mean, state_second, parent_mean, parent_second, n_parent, eff_n):
    state_var = np.maximum(state_second - np.square(state_mean), 0)
    parent_var = np.maximum(parent_second - np.square(parent_mean), 0)
    if not np.isfinite(eff_n) or eff_n <= 1 or n_parent <= 1:
        return np.full(len(state_mean), np.nan)
    se = np.sqrt((state_var / eff_n) + (parent_var / n_parent))
    z = np.full(len(state_mean), np.nan)
    valid = np.isfinite(se) & (se > 0)
    z[valid] = (state_mean[valid] - parent_mean[valid]) / se[valid]
    p_value = np.full(len(state_mean), np.nan)
    p_value[valid] = 2 * stats.norm.sf(np.abs(z[valid]))
    return p_value


def state_weight_matrix(activity, cells):
    cell_pos = {cell: i for i, cell in enumerate(cells)}
    columns = []
    col_meta = []
    for state, group in activity.groupby('state_name', sort=True):
        positions = group['cell_id'].map(cell_pos).dropna().astype(int).to_numpy()
        for weight_label, weight_col in weight_specs:
            weights = pd.to_numeric(group[weight_col], errors='coerce').fillna(0.0).to_numpy()
            keep = weights > 0
            if keep.any():
                columns.append(sparse.csr_matrix((weights[keep], (positions[keep], np.zeros(int(keep.sum()), dtype=int))), shape=(len(cells), 1)))
            else:
                columns.append(sparse.csr_matrix((len(cells), 1)))
            col_meta.append((state, weight_label))
    return sparse.hstack(columns, format='csr'), col_meta


def state_expression_table(cp10k, cp10k_sq, parent_mean, parent_second, n_parent, weight_matrix, col_meta, genes, tissue, cell_type):
    denom = np.asarray(weight_matrix.sum(axis=0)).ravel()
    denom_safe = np.where(denom > 0, denom, np.nan)
    state_means = np.asarray((weight_matrix.T @ cp10k).multiply(1 / denom_safe[:, None]).todense())
    state_seconds = np.asarray((weight_matrix.T @ cp10k_sq).multiply(1 / denom_safe[:, None]).todense())
    weight_square_sums = np.asarray(weight_matrix.power(2).sum(axis=0)).ravel()
    eff_n = np.full(len(denom), np.nan)
    valid = weight_square_sums > 0
    eff_n[valid] = denom[valid] ** 2 / weight_square_sums[valid]

    frames = []
    for col_idx, (state, weight_label) in enumerate(col_meta):
        state_mean = state_means[col_idx]
        state_second = state_seconds[col_idx]
        with np.errstate(invalid='ignore', divide='ignore'):
            log2fc = np.log2((state_mean + pseudocount) / (parent_mean + pseudocount))
        low_mean = (state_mean < min_mean_for_log2fc) & (parent_mean < min_mean_for_log2fc)
        log2fc[low_mean] = np.nan
        p_value = weighted_vs_parent_p_values(state_mean, state_second, parent_mean, parent_second, n_parent, eff_n[col_idx])
        frames.append(pd.DataFrame({
            'gene': genes,
            'tissue': tissue,
            'cell_type': cell_type,
            'state_name': state,
            'state_weight_type': weight_label,
            'weighted_mean_expression': state_mean,
            'log10_cpk': np.log10(np.clip(state_mean, 0, None) + 1.0),
            'log2fc_weighted_vs_all_parent': log2fc,
            'p_value': p_value,
        }))
    return pd.concat(frames, ignore_index=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--raw-10x-dir', required=True)
    parser.add_argument('--expression-value-type', required=True, choices=sorted(VALUE_TYPES))
    parser.add_argument('--metadata', required=True)
    parser.add_argument('--cell-state-activity', required=True)
    parser.add_argument('--cell-type-col', required=True)
    parser.add_argument('--api-minimal-output', action='store_true')
    parser.add_argument('--out-dir', required=True)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    counts, cells, genes = load_10x(args.raw_10x_dir)
    cp10k = linearize_expression_matrix(counts, args.expression_value_type)
    cp10k_sq = sparse_square(cp10k)

    metadata = pd.read_csv(args.metadata, sep='\t', compression='infer').drop_duplicates('cell_id').set_index('cell_id').reindex(cells)
    tissue = str(metadata['tissue'].iloc[0])
    cell_type = str(metadata[args.cell_type_col].iloc[0])

    parent_mean = np.asarray(cp10k.mean(axis=0)).ravel()
    parent_second = np.asarray(cp10k_sq.mean(axis=0)).ravel()
    n_parent = len(cells)

    cell_type_expression = pd.DataFrame({
        'gene': genes,
        'tissue': tissue,
        'cell_type': cell_type,
        'weighted_mean_expression': parent_mean,
        'log10_cpk': np.log10(np.clip(parent_mean, 0, None) + 1.0),
        'log2fc_weighted_vs_all_parent': np.nan,
        'p_value': np.nan,
    })

    activity = pd.read_csv(args.cell_state_activity, sep='\t', compression='infer')
    activity = activity.loc[activity['state_type'].eq('biological') & activity['cell_id'].isin(cells)]
    weight_matrix, col_meta = state_weight_matrix(activity, cells)
    state_expression = state_expression_table(cp10k, cp10k_sq, parent_mean, parent_second, n_parent, weight_matrix, col_meta, genes, tissue, cell_type)

    cell_type_expression.to_csv(f'{args.out_dir}/all_gene_cell_type_expression_cp10k.tsv.gz', sep='\t', index=False, compression='gzip')
    state_expression.to_csv(f'{args.out_dir}/all_gene_state_expression_specificity_cp10k.tsv.gz', sep='\t', index=False, compression='gzip')

    with open(f'{args.out_dir}/state_expression_summary.json', 'w') as f:
        json.dump({
            'raw_10x_dir': args.raw_10x_dir,
            'cell_state_activity': args.cell_state_activity,
            'n_cells': len(cells),
            'n_genes': len(genes),
            'n_states': int(activity['state_name'].nunique()),
            'expression_value_type': args.expression_value_type,
            'weight_types': [label for label, _ in weight_specs],
            'timestamp': datetime.now().isoformat(timespec='seconds'),
        }, f, indent=2)


if __name__ == '__main__':
    main()
