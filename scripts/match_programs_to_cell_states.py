#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import re
from datetime import datetime

import numpy as np
import pandas as pd
from scipy import stats

random_seed = 1
gsea_weight = 1.0
min_marker_overlap = 3
min_marker_coverage = 0.2


def norm_col(value):
    return re.sub(r'_+', '_', re.sub(r'[^a-z0-9]+', '_', str(value).strip().lower())).strip('_')


def display_label(value):
    text = re.sub(r'\s+', ' ', str(value).replace('_', ' ')).strip()
    return text[:1].upper() + text[1:] if text else ''


def bh_fdr(values):
    p = pd.to_numeric(values, errors='coerce')
    out = pd.Series(np.nan, index=p.index, dtype=float)
    valid = p.notna()
    if valid.sum() == 0:
        return out
    ranked = p.loc[valid].sort_values()
    n = len(ranked)
    q = ranked.to_numpy(float) * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out.loc[ranked.index] = np.clip(q, 0, 1)
    return out


def read_gmt(path, state_type):
    rows = []
    with open(path) as handle:
        for line in handle:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 3:
                continue
            markers = list(dict.fromkeys(g for g in parts[2:] if g))
            rows.append({'state_id': parts[0], 'state_label': display_label(parts[0]), 'state_type': state_type, 'markers': markers, 'n_state_markers': len(markers)})
    return pd.DataFrame(rows)


def load_program_loadings(path):
    frame = pd.read_csv(path, sep='\t', compression='infer')
    gene_col = frame.columns[0]
    out = frame.melt(id_vars=[gene_col], var_name='program_id', value_name='loading').rename(columns={gene_col: 'gene'})
    out['loading'] = pd.to_numeric(out['loading'], errors='coerce')
    return out.dropna(subset=['loading'])


def load_program_cell_activity(path):
    frame = pd.read_csv(path, sep='\t', compression='infer')
    return frame[['cell_id', 'program_id', 'program_activity']].dropna(subset=['program_activity'])


def gsea_es_for_hit_indices(loadings, hit_indices):
    n = len(loadings)
    nh = len(hit_indices)
    hit_weights = np.abs(loadings[hit_indices]) ** gsea_weight
    hit_total = hit_weights.sum()
    hit_weights = hit_weights / hit_total if hit_total > 0 else np.ones(nh) / nh
    increments = np.full(n, -1.0 / (n - nh))
    increments[hit_indices] = hit_weights
    running = np.cumsum(increments)
    return float(max(0.0, running.max()))


def marker_enrichment_for_pair(program, loadings, state, rng, tissue, cell_type, gsea_permutations):
    loadings = loadings.sort_values(['loading', 'gene'], ascending=[False, True])
    genes = loadings['gene'].tolist()
    loading_values = loadings['loading'].to_numpy(float)
    universe = set(genes)
    markers = state['markers']
    hits = set(g for g in markers if g in universe)
    hit_positions = np.array([i for i, gene in enumerate(genes) if gene in hits], dtype=int)
    n_program = len(genes)
    n_markers = len(markers)
    nh = len(hit_positions)
    coverage = nh / n_markers if n_markers else np.nan
    row = {
        'tissue': tissue,
        'cell_type': cell_type,
        'program_id': program,
        'state_id': state['state_id'],
        'state_label': state['state_label'],
        'state_type': state['state_type'],
        'n_program_genes': n_program,
        'n_state_markers': n_markers,
        'n_state_markers_in_program_universe': nh,
        'marker_coverage_fraction': coverage,
        'gsea_nes': np.nan,
        'gsea_p': np.nan,
        'match_status': 'ok',
    }
    if nh < min_marker_overlap or (not np.isnan(coverage) and coverage < min_marker_coverage) or nh == n_program:
        row['match_status'] = 'insufficient_marker_coverage'
        return row
    es = gsea_es_for_hit_indices(loading_values, hit_positions)
    null_es = np.array([gsea_es_for_hit_indices(loading_values, np.sort(rng.choice(n_program, size=nh, replace=False))) for _ in range(gsea_permutations)])
    mean_null = float(np.nanmean(null_es))
    row['gsea_nes'] = es / mean_null if mean_null > 0 else np.nan
    row['gsea_p'] = (1 + int(np.sum(null_es >= es))) / (1 + gsea_permutations)
    return row


def compute_marker_enrichment(program_loadings, states, tissue, cell_type, gsea_permutations):
    rng = np.random.default_rng(random_seed)
    rows = [
        marker_enrichment_for_pair(program, loadings, state, rng, tissue, cell_type, gsea_permutations)
        for program, loadings in program_loadings.groupby('program_id', sort=False)
        for _, state in states.iterrows()
    ]
    out = pd.DataFrame(rows)
    out['gsea_q'] = bh_fdr(out['gsea_p'])
    return out


def compute_cell_correlations(program_activity, cell_state_activity):
    activity = cell_state_activity[['cell_id', 'state_name', 'state_activity_weight_gradient']].rename(columns={'state_name': 'state_id', 'state_activity_weight_gradient': 'state_activity'})
    rows = []
    for program, prog in program_activity.groupby('program_id', sort=False):
        for state_id, st in activity.groupby('state_id', sort=False):
            merged = prog.merge(st[['cell_id', 'state_activity']], on='cell_id', how='inner').dropna()
            row = {'program_id': program, 'state_id': state_id, 'cell_spearman_r': np.nan, 'cell_spearman_p': np.nan}
            if len(merged) >= 3 and merged['program_activity'].nunique() > 1 and merged['state_activity'].nunique() > 1:
                corr = stats.spearmanr(merged['program_activity'], merged['state_activity'])
                row['cell_spearman_r'] = float(corr.statistic)
                row['cell_spearman_p'] = float(corr.pvalue)
            rows.append(row)
    out = pd.DataFrame(rows)
    out['cell_spearman_q'] = bh_fdr(out['cell_spearman_p'])
    return out


def neglog10_q(q):
    if pd.isna(q):
        return 0.0
    return min(50.0, -np.log10(max(float(q), 1e-300)))


def qc_caveat_for_state(state_id):
    text = norm_col(state_id)
    if any(x in text for x in ['ribosomal', 'translation']):
        return 'ribosomal_or_translation'
    if any(x in text for x in ['mitochondrial', 'apoptosis', 'cell_death', 'dying']):
        return 'mitochondrial_or_dying_cell'
    if any(x in text for x in ['offtarget', 'off_target', 'lineage']):
        return 'off_target_identity'
    if any(x in text for x in ['ambient', 'contamination', 'doublet']):
        return 'ambient_or_contamination'
    if any(x in text for x in ['heat_shock', 'dissociation']):
        return 'heat_shock_or_dissociation'
    if any(x in text for x in ['immediate_early', 'fos', 'jun']):
        return 'immediate_early'
    return 'none'


def interpretation_for(match_class):
    return {
        'strong_state_match': 'Program loading genes are enriched for state markers and program activity tracks state activity.',
        'gene_only_state_match': 'Program matches state markers but cell-level coactivity is weak or unavailable.',
        'cell_only_coactivity': 'Program activity tracks state activity but marker overlap is weak.',
        'qc_dominated': 'Program is QC-dominated; do not label as biological without review.',
        'mixed_state_qc': 'Program has both biological state and QC/artifact evidence.',
        'insufficient_marker_coverage': 'Program has insufficient marker coverage for this state.',
        'unmatched': 'Program has no strong curated-state or QC match.',
    }.get(match_class, 'Program has no strong curated-state or QC match.')


def build_summary(marker, corr):
    summary = marker.merge(corr.rename(columns={'cell_spearman_r': 'cell_spearman_r_gradient', 'cell_spearman_q': 'cell_spearman_q_gradient'}), on=['program_id', 'state_id'], how='left')
    summary['best_gene_level_score'] = pd.to_numeric(summary['gsea_nes'], errors='coerce').clip(lower=0).fillna(0) * summary['gsea_q'].map(neglog10_q)
    summary['best_cell_level_score'] = pd.to_numeric(summary['cell_spearman_r_gradient'], errors='coerce').clip(lower=0).fillna(0) * summary['cell_spearman_q_gradient'].map(neglog10_q)
    summary['combined_match_score'] = summary['best_gene_level_score'] + summary['best_cell_level_score']
    summary['qc_caveat'] = np.where(summary['state_type'].eq('qc_state'), summary['state_id'].map(qc_caveat_for_state), 'none')

    best_curated = summary.loc[summary['state_type'].eq('curated_state')].sort_values('combined_match_score', ascending=False).drop_duplicates('program_id')
    best_qc = summary.loc[summary['state_type'].eq('qc_state')].sort_values('combined_match_score', ascending=False).drop_duplicates('program_id')
    best_curated_score = best_curated.set_index('program_id')['combined_match_score'].to_dict()
    best_qc_score = best_qc.set_index('program_id')['combined_match_score'].to_dict()
    best_qc_sig = best_qc.set_index('program_id')['gsea_q'].to_dict()
    best_qc_cell_q = best_qc.set_index('program_id')['cell_spearman_q_gradient'].to_dict()

    classes = []
    for _, row in summary.iterrows():
        program = row['program_id']
        gsea_q = row.get('gsea_q')
        cell_q = row.get('cell_spearman_q_gradient')
        cell_r = row.get('cell_spearman_r_gradient')
        marker_ok = row.get('match_status') != 'insufficient_marker_coverage'
        is_qc = row['state_type'] == 'qc_state'
        qc_dominated = best_qc_score.get(program, -np.inf) > best_curated_score.get(program, -np.inf) and (
            (not pd.isna(best_qc_sig.get(program, np.nan)) and best_qc_sig[program] <= 0.05)
            or (not pd.isna(best_qc_cell_q.get(program, np.nan)) and best_qc_cell_q[program] <= 0.05)
        )
        strong_curated = (
            row['state_type'] == 'curated_state'
            and marker_ok
            and not pd.isna(gsea_q)
            and gsea_q <= 0.05
            and row.get('gsea_nes', 0) > 0
            and not pd.isna(cell_r)
            and cell_r >= 0.20
        )
        strong_qc_for_program = qc_dominated and best_curated_score.get(program, 0) > 0 and best_qc_score.get(program, 0) > 0
        if row.get('match_status') == 'insufficient_marker_coverage':
            cls = 'insufficient_marker_coverage'
        elif is_qc and qc_dominated:
            cls = 'qc_dominated'
        elif strong_curated and strong_qc_for_program:
            cls = 'mixed_state_qc'
        elif strong_curated:
            cls = 'strong_state_match'
        elif row['state_type'] == 'curated_state' and marker_ok and not pd.isna(gsea_q) and gsea_q <= 0.05:
            cls = 'gene_only_state_match'
        elif row['state_type'] == 'curated_state' and not pd.isna(cell_q) and cell_q <= 0.05 and not pd.isna(cell_r) and cell_r >= 0.30 and (pd.isna(gsea_q) or gsea_q > 0.05 or not marker_ok):
            cls = 'cell_only_coactivity'
        else:
            cls = 'unmatched'
        classes.append(cls)
    summary['match_class'] = classes
    summary['interpretation'] = summary['match_class'].map(interpretation_for)
    return summary


def build_qc_summary(summary):
    qc = summary.loc[summary['state_type'].eq('qc_state')]
    if qc.empty:
        return pd.DataFrame(columns=['program_id', 'best_qc_state_id', 'best_qc_label', 'best_qc_gsea_q', 'best_qc_gsea_nes', 'best_qc_cell_spearman_r', 'best_qc_cell_spearman_q', 'qc_combined_match_score', 'qc_caveat', 'qc_recommendation'])
    rows = []
    for program, group in qc.groupby('program_id', sort=False):
        best = group.sort_values('combined_match_score', ascending=False).iloc[0]
        caveats = sorted(set(x for x in group.loc[group['combined_match_score'] > 0, 'qc_caveat'].dropna() if x != 'none'))
        caveat = 'mixed_qc' if len(caveats) > 1 else caveats[0] if caveats else best['qc_caveat']
        significant = (not pd.isna(best['gsea_q']) and best['gsea_q'] <= 0.05) or (not pd.isna(best['cell_spearman_q_gradient']) and best['cell_spearman_q_gradient'] <= 0.05)
        recommendation = 'suppress_or_hide_by_default' if best['match_class'] == 'qc_dominated' and significant else 'review' if significant else 'pass'
        rows.append({
            'program_id': program,
            'best_qc_state_id': best['state_id'],
            'best_qc_label': best['state_label'],
            'best_qc_gsea_q': best['gsea_q'],
            'best_qc_gsea_nes': best['gsea_nes'],
            'best_qc_cell_spearman_r': best['cell_spearman_r_gradient'],
            'best_qc_cell_spearman_q': best['cell_spearman_q_gradient'],
            'qc_combined_match_score': best['combined_match_score'],
            'qc_caveat': caveat,
            'qc_recommendation': recommendation,
        })
    return pd.DataFrame(rows)


def build_label_suggestions(summary):
    rows = []
    curated = summary.loc[summary['state_type'].eq('curated_state')]
    qc = summary.loc[summary['state_type'].eq('qc_state')]
    for program in sorted(summary['program_id'].unique()):
        c = curated.loc[curated['program_id'].eq(program)].sort_values('combined_match_score', ascending=False)
        q = qc.loc[qc['program_id'].eq(program)].sort_values('combined_match_score', ascending=False)
        best_c = c.iloc[0] if not c.empty else None
        best_q = q.iloc[0] if not q.empty else None
        qc_caveat = best_q['qc_caveat'] if best_q is not None else 'none'
        if best_q is not None and best_q['match_class'] == 'qc_dominated':
            label = f"{qc_caveat.replace('_', ' ')}/QC program" if qc_caveat != 'none' else f"{best_q['state_label']} QC program"
            quality = 'qc_or_artifact'
        elif best_c is not None and best_c['match_class'] == 'strong_state_match':
            label = f"{best_c['state_label']}-like program"
            quality = 'high_confidence_biological'
        elif best_c is not None and best_c['match_class'] == 'gene_only_state_match':
            label = f"{best_c['state_label']}-enriched program"
            quality = 'exploratory_biological'
        elif best_c is not None and best_c['match_class'] == 'cell_only_coactivity':
            label = f"{best_c['state_label']}-coactive program"
            quality = 'exploratory_biological'
        elif best_c is not None and best_c['match_class'] == 'mixed_state_qc':
            label = f"{best_c['state_label']}-like mixed QC program"
            quality = 'mixed_state_qc'
        else:
            label = 'unmatched data-driven program'
            quality = 'unmatched'
        rows.append({
            'program_id': program,
            'best_curated_state_id': best_c['state_id'] if best_c is not None else '',
            'best_curated_state_label': best_c['state_label'] if best_c is not None else '',
            'best_curated_match_class': best_c['match_class'] if best_c is not None else '',
            'best_qc_state_id': best_q['state_id'] if best_q is not None else '',
            'best_qc_label': best_q['state_label'] if best_q is not None else '',
            'qc_caveat': qc_caveat,
            'suggested_program_label': label,
            'suggested_program_quality_class': quality,
        })
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--program-loadings', required=True)
    parser.add_argument('--state-gmt', required=True)
    parser.add_argument('--qc-gmt', required=True)
    parser.add_argument('--program-cell-activity', required=True)
    parser.add_argument('--cell-state-activity', required=True)
    parser.add_argument('--tissue', required=True)
    parser.add_argument('--cell-type', required=True)
    parser.add_argument('--gsea-permutations', type=int, required=True)
    parser.add_argument('--out-dir', required=True)
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    program_loadings = load_program_loadings(args.program_loadings)
    states = pd.concat([read_gmt(args.state_gmt, 'curated_state'), read_gmt(args.qc_gmt, 'qc_state')], ignore_index=True)
    program_activity = load_program_cell_activity(args.program_cell_activity)
    cell_state_activity = pd.read_csv(args.cell_state_activity, sep='\t', compression='infer')

    marker = compute_marker_enrichment(program_loadings, states, args.tissue, args.cell_type, args.gsea_permutations)
    corr = compute_cell_correlations(program_activity, cell_state_activity)
    summary = build_summary(marker, corr)
    qc_summary = build_qc_summary(summary)
    labels = build_label_suggestions(summary)

    marker.to_csv(f'{args.out_dir}/program_state_marker_enrichment.tsv.gz', sep='\t', index=False, compression='gzip')
    summary.to_csv(f'{args.out_dir}/program_state_match_summary.tsv.gz', sep='\t', index=False, compression='gzip')
    qc_summary.to_csv(f'{args.out_dir}/program_qc_match_summary.tsv.gz', sep='\t', index=False, compression='gzip')
    labels.to_csv(f'{args.out_dir}/program_label_suggestions.tsv.gz', sep='\t', index=False, compression='gzip')

    with open(f'{args.out_dir}/run_summary.json', 'w') as f:
        json.dump({
            'program_loadings': args.program_loadings,
            'state_gmt': args.state_gmt,
            'qc_gmt': args.qc_gmt,
            'tissue': args.tissue,
            'cell_type': args.cell_type,
            'n_programs': int(program_loadings['program_id'].nunique()),
            'n_curated_states': int(states['state_type'].eq('curated_state').sum()),
            'n_qc_states': int(states['state_type'].eq('qc_state').sum()),
            'timestamp': datetime.now().isoformat(timespec='seconds'),
        }, f, indent=2)


if __name__ == '__main__':
    main()
