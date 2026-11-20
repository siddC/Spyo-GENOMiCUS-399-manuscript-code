#!/usr/bin/env python3
"""
Spyo-GENOMiCUS-399 Results Section 3 reproducibility script.
"Phylons are concordant with established classification metrics"

For each phylon, computes ARI, NMI, purity (PPV), sensitivity (recall),
and F1 against the relevant ground-truth classification:

    ST28-emm1, ST52-emm28, ST39-emm4, ST63-emm77, ST101-emm89:
        vs emm-type AND vs MLST sequence type
    ST36:
        vs MLST only (lineage spans emm12 and emm82 via emm-switching)
    emm3/6/5:
        vs emm-type membership in {emm3, emm5, emm6}
    emm-patternD-other:
        vs emm-pattern D (clusters D1-D5 plus clade Y)
    emm-patternE-other:
        vs emm-pattern E (clusters E1-E6)

Treats each phylon as a binary classifier ("in phylon" vs "not"); ground
truth label is binary in the same way ("has truth label" vs "not"). ARI
and NMI are computed as 2x2 partition-level metrics over all 399 genomes.

Inputs (in INPUT_DIR; default "data/", override with argv[1] or $SPYO_DATA):
    A_bin_panta_cdhit_Spyo.csv     binarized A matrix (11 phylons x 399)
    L_bin_panta_cdhit_Spyo.csv     binarized L matrix (1167 genes x 11)
    car_metadata_Sypo.csv          per-genome metadata (399 rows)
    gene_annot_panta_Spyo.csv      PanTA gene annotations
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.metrics import (adjusted_rand_score, normalized_mutual_info_score,
                             confusion_matrix)

import os, sys

# Input/output directories. Defaults assume the script is run from the repository
# root with the Zenodo matrices unpacked into data/.
INPUT_DIR  = Path(sys.argv[1] if len(sys.argv) > 1
                  else os.environ.get('SPYO_DATA', 'data'))
OUTPUT_DIR = Path(sys.argv[2] if len(sys.argv) > 2
                  else os.environ.get('SPYO_RESULTS', 'results'))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

for _f in ('A_bin_panta_cdhit_Spyo.csv', 'L_bin_panta_cdhit_Spyo.csv',
           'car_metadata_Sypo.csv', 'gene_annot_panta_Spyo.csv'):
    if not (INPUT_DIR / _f).exists():
        sys.exit(f'Missing input: {INPUT_DIR / _f}\n'
                 f'Download the Zenodo deposit first:\n'
                 f'    python scripts/zenodo_fetch.py 20403112 -o {INPUT_DIR}')

A   = pd.read_csv(INPUT_DIR / 'A_bin_panta_cdhit_Spyo.csv', index_col=0)
L   = pd.read_csv(INPUT_DIR / 'L_bin_panta_cdhit_Spyo.csv', index_col=0)
md  = pd.read_csv(INPUT_DIR / 'car_metadata_Sypo.csv', low_memory=False)
ann = pd.read_csv(INPUT_DIR / 'gene_annot_panta_Spyo.csv', low_memory=False)

md['genome_id'] = md['genome_id'].astype(str)
md = md.set_index('genome_id').loc[A.columns]   # align with A's 399 columns

# ============================================================================
# Sanity checks on input dimensions
# ============================================================================
print('=' * 78)
print('SANITY CHECKS')
print('=' * 78)
print(f'A matrix shape:    {A.shape}     (expect (11, 399))')
print(f'L matrix shape:    {L.shape}     (expect (1167, 11))')
print(f'Metadata rows:     {len(md)}     (expect 399)')
print(f'Annotation rows:   {len(ann)}    (expect 8864 within partition)')

# Pangenome breakdown by frequency
n_genomes = A.shape[1]
core_thr  = int(np.ceil(0.979 * n_genomes))   # >=391
rare_thr  = int(np.ceil(0.082 * n_genomes))   # <33

ann_freq = ann['No. isolates'].astype(int)
n_core  = (ann_freq >= core_thr).sum()
n_acc   = ((ann_freq >= rare_thr) & (ann_freq < core_thr)).sum()
n_rare  = ((ann_freq >= 1) & (ann_freq < rare_thr)).sum()
n_zero  = (ann_freq == 0).sum()
total_in_partition = n_core + n_acc + n_rare

print(f'\nPangenome breakdown (from gene_annot No. isolates):')
print(f'  Core genes (>={core_thr} strains):       {n_core:>5}   (expect 1292)')
print(f'  Accessory genes ({rare_thr}-{core_thr-1} strains):  {n_acc:>5}   '
      f'(expect 1167)')
print(f'  Rare genes (1-{rare_thr-1} strains):        {n_rare:>5}   '
      f'(expect 6405)')
print(f'  Genes with No. isolates == 0:        {n_zero:>5}')
print(f'  Sum (core+accessory+rare):           {total_in_partition:>5}   '
      f'(expect 8864)')
print(f'  Annotation rows TOTAL:               {len(ann):>5}')
print(f'  Annotation rows OUTSIDE partition:   '
      f'{len(ann) - total_in_partition:>5}')

# Sanity: do L matrix gene IDs all live in the accessory bucket of ann?
L_in_ann = L.index.isin(ann['Gene'])
print(f'\nL matrix gene-ID coverage in annotation file: '
      f'{L_in_ann.sum()} / {len(L)}')

# ============================================================================
# Build ground-truth labels
# ============================================================================
def parse_emm(s):
    """Extract 'emmNN' canonical token from a Predicted emm-type string."""
    if pd.isna(s):
        return None
    m = pd.Series([str(s)]).str.extract(r'[Ee][Mm][Mm]_?(\d+)',
                                         expand=False).iloc[0]
    return f'emm{m}' if pd.notna(m) else None

def cluster_to_pattern(c):
    """Map an emmtyper emm-cluster string to an emm-pattern code."""
    if pd.isna(c) or not isinstance(c, str):
        return None
    c = c.strip()
    if c.startswith('A-C'):
        return 'A-C'
    cl = c.lower()
    if 'clade y' in cl or c == 'Y':
        return 'Y'
    if c.startswith('D'):
        return 'D'
    if c.startswith('E'):
        return 'E'
    return None

md['emm_typed']    = md['Predicted emm-type'].apply(parse_emm)
md['emm_pattern']  = md['emm cluster'].apply(cluster_to_pattern)
md['ST']           = md['mlst'].astype(str).replace('nan', None)

print(f'\nGround-truth coverage across 399 genomes:')
print(f'  emm-typed (Predicted emm-type parsed): {md["emm_typed"].notna().sum()}')
print(f'  emm-pattern (from emm cluster):        '
      f'{md["emm_pattern"].notna().sum()}')
print(f'  MLST ST assignment:                    {md["ST"].notna().sum()}')

# ============================================================================
# Per-phylon binary concordance
# ============================================================================
print()
print('=' * 78)
print('PER-PHYLON CONCORDANCE')
print('=' * 78)

def binary_concordance(phylon_member, truth_member, label):
    """
    Compute ARI, NMI, and confusion-matrix metrics between two binary
    vectors. Strains with missing truth label are dropped.
    """
    mask = phylon_member.notna() & truth_member.notna()
    a = phylon_member[mask].astype(int)
    b = truth_member[mask].astype(int)
    n = len(a)
    if n == 0 or a.sum() == 0 or b.sum() == 0:
        return None
    ari = adjusted_rand_score(a, b)
    nmi = normalized_mutual_info_score(a, b)
    tn, fp, fn, tp = confusion_matrix(b, a, labels=[0, 1]).ravel()
    sens = tp / (tp + fn) if (tp + fn) else 0.0   # recall on truth label
    ppv  = tp / (tp + fp) if (tp + fp) else 0.0   # purity of phylon
    f1   = 2 * sens * ppv / (sens + ppv) if (sens + ppv) else 0.0
    return {
        'comparison': label,
        'n_evaluable': n,
        'phylon_n': int(a.sum()),
        'truth_n': int(b.sum()),
        'TP': int(tp), 'FP': int(fp), 'FN': int(fn), 'TN': int(tn),
        'PPV_purity_pct': round(ppv * 100, 1),
        'sensitivity_pct': round(sens * 100, 1),
        'F1_pct': round(f1 * 100, 1),
        'ARI': round(ari, 3),
        'NMI': round(nmi, 3),
    }

# Display labels (manuscript uses "patternD/E", A/L matrices use "clusterD/E")
DISPLAY = {
    'emm-clusterD-other': 'emm-patternD-other',
    'emm-clusterE-other': 'emm-patternE-other',
}

COMPARISONS = [
    ('ST28-emm1',          'emm_typed',   {'emm1'},                 'vs emm1 (emm-type)'),
    ('ST28-emm1',          'ST',          {'28'},                   'vs ST28 (MLST)'),
    ('ST52-emm28',         'emm_typed',   {'emm28'},                'vs emm28 (emm-type)'),
    ('ST52-emm28',         'ST',          {'52'},                   'vs ST52 (MLST)'),
    ('ST39-emm4',          'emm_typed',   {'emm4'},                 'vs emm4 (emm-type)'),
    ('ST39-emm4',          'ST',          {'39'},                   'vs ST39 (MLST)'),
    ('ST63-emm77',         'emm_typed',   {'emm77'},                'vs emm77 (emm-type)'),
    ('ST63-emm77',         'ST',          {'63'},                   'vs ST63 (MLST)'),
    ('ST101-emm89',        'emm_typed',   {'emm89'},                'vs emm89 (emm-type)'),
    ('ST101-emm89',        'ST',          {'101'},                  'vs ST101 (MLST)'),
    ('ST36',               'ST',          {'36'},                   'vs ST36 (MLST)'),
    ('emm3/6/5',           'emm_typed',   {'emm3','emm5','emm6'},   'vs emm3/5/6 (emm-type)'),
    ('emm-clusterD-other', 'emm_pattern', {'D','Y'},                'vs pattern D (incl. clade Y)'),
    ('emm-clusterE-other', 'emm_pattern', {'E'},                    'vs pattern E'),
]

results = []
for phylon, truth_col, truth_val, label in COMPARISONS:
    if phylon not in A.index:
        print(f'SKIP: {phylon} not in A.index')
        continue
    phylon_member = pd.Series(A.loc[phylon].astype(int).values,
                              index=A.columns, dtype='Int64')
    truth_member  = md[truth_col].apply(
        lambda x: 1 if x in truth_val else (0 if pd.notna(x) else None))
    res = binary_concordance(phylon_member, truth_member, label)
    if res is not None:
        res['phylon'] = DISPLAY.get(phylon, phylon)
        results.append(res)

results_df = pd.DataFrame(results)
cols_order = ['phylon', 'comparison', 'n_evaluable', 'phylon_n', 'truth_n',
              'TP', 'FP', 'FN', 'TN',
              'PPV_purity_pct', 'sensitivity_pct', 'F1_pct',
              'ARI', 'NMI']
results_df = results_df[cols_order]

# Pretty print
with pd.option_context('display.max_rows', None, 'display.width', 200,
                       'display.max_colwidth', 40):
    print(results_df.to_string(index=False))

results_df.to_csv(OUTPUT_DIR / 'TableS5_per_phylon_concordance.csv', index=False)

# ============================================================================
# Global multi-class concordance (matches manuscript Results text)
# Manuscript Methods: "only strains assigned to exactly one taxonomic phylon
# (excluding mobilons) were included."
# ============================================================================
print()
print('=' * 78)
print('GLOBAL MULTI-CLASS CONCORDANCE (matches manuscript ARI/NMI text)')
print('=' * 78)

LINEAGE_PHYLONS = ['ST28-emm1', 'ST52-emm28', 'ST39-emm4', 'emm3/6/5',
                   'ST36', 'ST63-emm77', 'ST101-emm89']
COMPOSITE_PHYLONS = ['emm-clusterD-other', 'emm-clusterE-other']
ALL_TAXON_PHYLONS = LINEAGE_PHYLONS + COMPOSITE_PHYLONS

def assign_partition(A_sub, strains):
    """Return a Series mapping strain -> phylon name (single assignment)."""
    sub = A_sub.loc[:, strains]
    out = {}
    for s in strains:
        col = sub[s]
        assigned = col.index[col == 1].tolist()
        out[s] = assigned[0] if len(assigned) == 1 else None
    return pd.Series(out)

def global_concordance(scope_label, scope_phylons, truth_col, truth_label):
    A_sub = A.loc[scope_phylons]
    # strains assigned to exactly one phylon in this scope
    counts = A_sub.sum(axis=0)
    eligible = counts[counts == 1].index.tolist()
    phylon_assignment = assign_partition(A_sub, eligible)
    truth = md.loc[eligible, truth_col]
    mask = truth.notna()
    if mask.sum() == 0:
        return None
    p = phylon_assignment[mask].astype(str).values
    t = truth[mask].astype(str).values
    return {
        'scope': scope_label,
        'truth': truth_label,
        'n_strains': int(mask.sum()),
        'n_phylon_classes': len(set(p)),
        'n_truth_classes':  len(set(t)),
        'ARI': round(adjusted_rand_score(t, p), 3),
        'NMI': round(normalized_mutual_info_score(t, p), 3),
    }

global_results = []
for scope_label, scope in [
    ('7 lineage-specific phylons', LINEAGE_PHYLONS),
    ('all 9 taxonomic phylons',     ALL_TAXON_PHYLONS),
]:
    for truth_col, truth_label in [('emm_typed',   'emm-type'),
                                    ('ST',          'MLST'),
                                    ('emm_pattern', 'emm-pattern')]:
        res = global_concordance(scope_label, scope, truth_col, truth_label)
        if res:
            global_results.append(res)

global_df = pd.DataFrame(global_results)
print('\nManuscript-reported values:')
print('  Lineage-specific only:  ARI=0.90/NMI=0.88 (MLST), ARI=0.92/NMI=0.89 (emm)')
print('  All 9 taxonomic phylons: ARI=0.74/NMI=0.75 (MLST), ARI=0.76/NMI=0.76 (emm)')
print()
print(global_df.to_string(index=False))
global_df.to_csv(OUTPUT_DIR / 'TableS6_global_concordance.csv', index=False)

print()
print('Files written:')
print(f'  {OUTPUT_DIR / "TableS5_per_phylon_concordance.csv"}')
print(f'  {OUTPUT_DIR / "TableS6_global_concordance.csv"}')
