#!/usr/bin/env python3
"""
Verification of all claims in the Results section:
"The Sda1 prophage defines the primary taxonomic split and is identified
 in ST63-emm77"

Checks every numerical claim against the binarized A matrix (phylon
assignments), P matrix (gene presence/absence in 399 genomes), and
tBLASTn results (sda1_hits.tsv).

Input files (expected in working directory or adjust paths):
    - A_bin_panta_cdhit_Spyo.csv
    - L_bin_panta_cdhit_Spyo.csv
    - df_genes_panta_cdhit_final.pickle
    - car_metadata_Sypo.csv
    - gene_annot_panta_Spyo.csv
    - sda1_hits.tsv  (tBLASTn outfmt 6, query = Q675N6, db = 399 genomes)

Author: Siddharth M. Chauhan
"""

import os
from pathlib import Path
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
# ---------------------------------------------------------------------------
# Input/output directories.  Defaults assume this script is run from the
# repository root with the Zenodo matrices unpacked into data/:
#     python scripts/sda1_results_verification.py
# Override either with an environment variable:
#     SPYO_DATA=/path/to/data SPYO_RESULTS=/path/to/out python scripts/sda1_results_verification.py
# ---------------------------------------------------------------------------
INPUT_DIR  = Path(os.environ.get("SPYO_DATA", "data"))
OUTPUT_DIR = Path(os.environ.get("SPYO_RESULTS", "results"))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

A_PATH     = INPUT_DIR / "A_bin_panta_cdhit_Spyo.csv"
L_PATH     = INPUT_DIR / "L_bin_panta_cdhit_Spyo.csv"
P_PATH     = INPUT_DIR / "df_genes_panta_cdhit_final.pickle"
META_PATH  = INPUT_DIR / "car_metadata_Sypo.csv"
ANNOT_PATH = INPUT_DIR / "gene_annot_panta_Spyo.csv"
BLAST_PATH = INPUT_DIR / "sda1_hits.tsv"

SDA1_PHYLONS = ["ST28-emm1", "ST36", "ST63-emm77"]
QUERY_LEN = 390   # Sda1/SdaD2 protein length (UniProt Q675N6)

# ══════════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════
print("Loading data...")
A = pd.read_csv(A_PATH, index_col=0)
L = pd.read_csv(L_PATH, index_col=0)
P = pd.read_pickle(P_PATH)
meta = pd.read_csv(META_PATH, low_memory=False)

blast_cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
              "qstart", "qend", "sstart", "send", "evalue", "bitscore"]
blast = pd.read_csv(BLAST_PATH, sep="\t", header=None, names=blast_cols)
blast["qcov"] = (blast["qend"] - blast["qstart"] + 1) / QUERY_LEN * 100

n_genomes = P.shape[1]
print(f"  {n_genomes} genomes, {P.shape[0]} genes in P matrix")

# ══════════════════════════════════════════════════════════════════════════════
# 1. PHYLON SIZES
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("[1] PHYLON SIZES")
print("=" * 72)
for p in SDA1_PHYLONS:
    n = int(A.loc[p].sum())
    print(f"    {p}: {n} strains")
total_sda1 = sum(int(A.loc[p].sum()) for p in SDA1_PHYLONS)
print(f"    Total Sda1+ phylon strains: {total_sda1}")

# ══════════════════════════════════════════════════════════════════════════════
# 2. SDA1 CARRIAGE WITHIN SDA1+ PHYLONS (P matrix)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("[2] sda1 PanTA CLUSTER CARRIAGE IN SDA1+ PHYLONS")
print("=" * 72)
sda1_row = P.loc["sda1"]
for p in SDA1_PHYLONS:
    strains = list(A.columns[A.loc[p] == 1])
    n = len(strains)
    c = int(sda1_row[strains].sum())
    print(f"    {p}: {c}/{n} ({c/n*100:.1f}%)")

# ══════════════════════════════════════════════════════════════════════════════
# 3. tBLASTn RESULTS
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("[3] tBLASTn ORTHOLOG SUMMARY")
print("=" * 72)
ortho = blast[(blast["pident"] >= 80) & (blast["qcov"] >= 80)]
n_ortho_contigs = ortho["sseqid"].nunique()
n_100pct = len(ortho[ortho["pident"] == 100])
n_sub100 = len(ortho[ortho["pident"] < 100])

print(f"    Orthologous hits (≥80% id, ≥80% qcov): {len(ortho)}")
print(f"    Unique contigs with orthologs: {n_ortho_contigs}")
print(f"    At 100% identity: {n_100pct}")
print(f"    At <100% identity: {n_sub100}")
if n_sub100 > 0:
    print(f"    Non-100% hits:")
    for _, r in ortho[ortho["pident"] < 100].iterrows():
        print(f"      {r['sseqid']}: {r['pident']:.1f}% id, "
              f"{r['length']}aa aligned")

# Contigs with multiple ortholog hits (double-prophage)
multi = ortho.groupby("sseqid").size()
doubles = multi[multi > 1]
if len(doubles) > 0:
    print(f"    Contigs with >1 orthologous hit: {len(doubles)}")
    for c in doubles.index:
        print(f"      {c} (n={doubles[c]} hits)")

# ══════════════════════════════════════════════════════════════════════════════
# 4. NON-SDA1 PHYLON CARRIAGE
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("[4] sda1 CARRIAGE IN NON-SDA1 PHYLONS")
print("=" * 72)
non_sda1 = [p for p in A.index if p not in SDA1_PHYLONS]
max_pct = 0
for p in non_sda1:
    strains = list(A.columns[A.loc[p] == 1])
    c = int(sda1_row[strains].sum())
    pct = c / len(strains) * 100
    max_pct = max(max_pct, pct)
    print(f"    {p:25s}: {c}/{len(strains)} ({pct:.1f}%)")
print(f"    Maximum: {max_pct:.1f}%")

# ══════════════════════════════════════════════════════════════════════════════
# 5. TOTAL CARRIER COUNTS AND RESTRICTION TO SDA1+ PHYLONS
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("[5] TOTAL sda1 CARRIERS AND PHYLON RESTRICTION")
print("=" * 72)
total_carriers = int(sda1_row.sum())
sda1_phylon_strains = set()
for p in SDA1_PHYLONS:
    sda1_phylon_strains |= set(A.columns[A.loc[p] == 1])
carriers = set(P.columns[P.loc["sda1"] == 1])
in_phylons = carriers & sda1_phylon_strains
outside = carriers - sda1_phylon_strains

print(f"    Total sda1 carriers (P matrix): {total_carriers}/{n_genomes}")
print(f"    In Sda1+ phylons: {len(in_phylons)}/{total_carriers} "
      f"({len(in_phylons)/total_carriers*100:.1f}%)")
print(f"    Sporadic (outside Sda1+ phylons): {len(outside)}")

# ══════════════════════════════════════════════════════════════════════════════
# 6. IDENTIFY 27-GENE SDA1 PROPHAGE CORE (from L matrix)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "=" * 72)
print("[6] 27-GENE SDA1 PROPHAGE CORE")
print("=" * 72)
other_phylons = [p for p in L.columns if p not in SDA1_PHYLONS]
in_all = (L[SDA1_PHYLONS] == 1).all(axis=1)
out_all = (L[other_phylons] == 0).all(axis=1)
sda1_genes = L.index[in_all & out_all].tolist()
print(f"    Genes in all 3 Sda1+ phylons, absent from all others: {len(sda1_genes)}")

# Carriage stats for the 27 genes in Sda1+ vs non-Sda1 strains
sda1_strain_list = list(sda1_phylon_strains)
non_strain_list = list(set(A.columns) - sda1_phylon_strains)
pcts_sda1 = []
pcts_non = []
for g in sda1_genes:
    c_s = int(P.loc[g, sda1_strain_list].sum())
    c_n = int(P.loc[g, non_strain_list].sum())
    pcts_sda1.append(c_s / len(sda1_strain_list) * 100)
    pcts_non.append(c_n / len(non_strain_list) * 100)

print(f"    Median carriage in Sda1+ (n={len(sda1_strain_list)}): "
      f"{np.median(pcts_sda1):.1f}% "
      f"(range {min(pcts_sda1):.1f}% – {max(pcts_sda1):.1f}%)")
print(f"    Median carriage in non-Sda1 (n={len(non_strain_list)}): "
      f"{np.median(pcts_non):.1f}% "
      f"(range {min(pcts_non):.1f}% – {max(pcts_non):.1f}%)")

print("\n" + "=" * 72)
print("VERIFICATION COMPLETE")
print("=" * 72)
