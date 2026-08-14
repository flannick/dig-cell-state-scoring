#!/usr/bin/env python3
import argparse
import gzip
import os
import re
import subprocess

import pandas as pd

program_output_columns = ['dataset', 'cell_type', 'model', 'factor', 'trait', 'beta', 'beta_uncorrected']
curated_output_columns = ['tissue', 'cell_type', 'state_name', 'trait', 'beta', 'beta_uncorrected']


def read_gmt_names(path):
    with open(path) as handle:
        return [line.split('\t', 1)[0] for line in handle if line.strip()]


def write_pigean_x_from_gmt(gmt_path, out_path):
    with open(gmt_path) as src, open(out_path, 'w') as dst:
        for line in src:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            genes = [g for g in parts[2:] if g]
            if genes:
                dst.write('\t'.join([parts[0]] + genes) + '\n')


def display_factor(value):
    m = re.search(r'(?:program[_-]?)?factor[_-]?([0-9]+)$', str(value), re.I)
    return f'Factor{m.group(1)}' if m else str(value)


def read_table(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame()
    return pd.read_csv(path, sep='\t', compression='infer', low_memory=False)


def pick_col(frame, names):
    by_lower = {str(c).lower(): c for c in frame.columns}
    for name in names:
        col = by_lower.get(name.lower())
        if col is not None:
            return col
    return None


def trait_is_auto_blacklisted(trait):
    return trait.startswith('HP_') or trait.startswith('exomes_') or 'gcat_' in trait or 'Orphanet' in trait


def generate_auto_trait_blacklist(multi_y_in, out_path, pheno_col):
    opener = gzip.open if multi_y_in.endswith('.gz') else open
    traits = set()
    with opener(multi_y_in, 'rt') as handle:
        header = handle.readline().rstrip('\n').split('\t')
        idx = header.index(pheno_col)
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if idx < len(parts) and trait_is_auto_blacklisted(parts[idx]):
                traits.add(parts[idx])
    with open(out_path, 'w') as f:
        f.writelines(f'{trait}\n' for trait in sorted(traits))
    return out_path


def keep_positive_beta_uncorrected(frame):
    if frame.empty or 'beta_uncorrected' not in frame.columns:
        return frame
    return frame[pd.to_numeric(frame['beta_uncorrected'], errors='coerce') > 0].copy()


def normalize_one(frame, cell_type, tissue, dataset, model, kind):
    if frame.empty:
        return pd.DataFrame()
    gene_set_col = pick_col(frame, ['Gene_Set', 'gene_set', 'gene_set_name', 'state_name', 'factor', 'set', 'name'])
    if gene_set_col is None:
        return pd.DataFrame()
    trait_col = pick_col(frame, ['trait', 'Trait', 'phenotype', 'Trait_Internal', 'y'])
    beta_col = pick_col(frame, ['beta', 'Beta'])
    beta_unc_col = pick_col(frame, ['beta_uncorrected', 'Beta_uncorrected', 'beta_uncorrected_orig'])
    out = pd.DataFrame(index=frame.index)
    out['gene_set_id'] = frame[gene_set_col].astype(str)
    out['trait'] = frame[trait_col].astype(str) if trait_col else ''
    out['beta'] = pd.to_numeric(frame[beta_col], errors='coerce') if beta_col else pd.NA
    out['beta_uncorrected'] = pd.to_numeric(frame[beta_unc_col], errors='coerce') if beta_unc_col else pd.NA
    if kind == 'program':
        out['dataset'] = dataset
        out['cell_type'] = cell_type
        out['model'] = model
        out['factor'] = out['gene_set_id'].map(display_factor)
        return out[program_output_columns]
    out['tissue'] = tissue
    out['cell_type'] = cell_type
    out['state_name'] = out['gene_set_id']
    return out[curated_output_columns]


def run_pigean_for_gmt(gmt_path, cell_type, out_dir, args, blacklist):
    run_dir = f'{out_dir}/{cell_type}'
    os.makedirs(run_dir, exist_ok=True)
    pigean_x = f'{run_dir}/input.x.tsv'
    write_pigean_x_from_gmt(gmt_path, pigean_x)
    stats_out = f'{run_dir}/gene_set_stats.debug.out.gz'
    cmd = [
        args.python, '-m', 'pigean', 'betas',
        '--X-in', pigean_x,
        '--multi-y-in', args.multi_y_in,
        '--multi-y-id-col', args.multi_y_id_col,
        '--multi-y-pheno-col', args.multi_y_pheno_col,
        '--multi-y-log-bf-col', args.multi_y_log_bf_col,
        '--multi-y-combined-col', args.multi_y_combined_col,
        '--multi-y-prior-col', args.multi_y_prior_col,
        '--multi-y-trait-blacklist-in', blacklist,
        '--gene-universe-in', args.gene_universe_in,
        '--gene-universe-id-col', '6',
        '--gene-universe-no-header',
        '--gene-set-stats-out', stats_out,
        '--params-out', f'{run_dir}/params.out.gz',
        '--log-file', f'{run_dir}/run.log',
        '--warnings-file', f'{run_dir}/warnings.log',
        '--output-detail', 'debug',
        '--deterministic',
        '--hide-progress',
        '--min-gene-set-size', '1',
        '--filter-gene-set-p', '1',
        '--max-gene-set-read-p', '1',
        '--no-filter-negative',
        '--max-no-write-gene-set-beta-uncorrected', '0',
        '--prune-gene-sets', '1.1',
        '--weighted-prune-gene-sets', '1.1',
    ]
    env = os.environ.copy()
    env['PYTHONPATH'] = args.pythonpath
    with open(f'{run_dir}/stdout.log', 'w') as stdout, open(f'{run_dir}/stderr.log', 'w') as stderr:
        subprocess.run(cmd, env=env, stdout=stdout, stderr=stderr, check=True)
    return read_table(stats_out)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gmt-dir', required=True)
    parser.add_argument('--out-dir', required=True)
    parser.add_argument('--combined-out', required=True)
    parser.add_argument('--kind', required=True, choices=['curated', 'program'])
    parser.add_argument('--tissue', required=True)
    parser.add_argument('--dataset', required=True)
    parser.add_argument('--model', required=True)
    parser.add_argument('--python', required=True)
    parser.add_argument('--pythonpath', required=True)
    parser.add_argument('--multi-y-in', required=True)
    parser.add_argument('--multi-y-id-col', required=True)
    parser.add_argument('--multi-y-pheno-col', required=True)
    parser.add_argument('--multi-y-log-bf-col', required=True)
    parser.add_argument('--multi-y-combined-col', required=True)
    parser.add_argument('--multi-y-prior-col', required=True)
    parser.add_argument('--trait-blacklist-in', required=True, choices=['auto'])
    parser.add_argument('--gene-universe-in', required=True)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(os.path.dirname(args.combined_out), exist_ok=True)

    blacklist = generate_auto_trait_blacklist(args.multi_y_in, f'{args.out_dir}/trait_blacklist_hp_exomes_gcat_orphanet.txt', args.multi_y_pheno_col)
    print(f'Using PIGEAN trait blacklist {blacklist}', flush=True)

    frames = []
    for gmt_name in sorted(f for f in os.listdir(args.gmt_dir) if f.endswith('.gmt')):
        gmt_path = f'{args.gmt_dir}/{gmt_name}'
        if os.path.getsize(gmt_path) == 0 or not read_gmt_names(gmt_path):
            continue
        cell_type = gmt_name[:-len('.gmt')]
        print(f'PIGEAN {args.kind} cell_type={cell_type}', flush=True)
        frame = run_pigean_for_gmt(gmt_path, cell_type, args.out_dir, args, blacklist)
        norm = keep_positive_beta_uncorrected(normalize_one(frame, cell_type, args.tissue, args.dataset, args.model, args.kind))
        print(f'{cell_type}: positive_beta_uncorrected_rows={len(norm)}', flush=True)
        if not norm.empty:
            frames.append(norm)

    columns = program_output_columns if args.kind == 'program' else curated_output_columns
    combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=columns)
    combined.to_csv(args.combined_out, sep='\t', index=False, compression='infer')


if __name__ == '__main__':
    main()
