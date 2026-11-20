#!/usr/bin/env python3
"""
Sda1 prophage enrichment analysis for Spyo-GENOMiCUS-399.

Identifies the 27-gene Sda1 prophage core shared by all three Sda1+
phylons (ST28-emm1, ST36, ST63-emm77) and absent from all other phylons
in the binarized L matrix, then tests each gene for enrichment in
Sda1+ strains (n=121) versus non-Sda1 strains (n=278) using Fisher's
exact test with Benjamini-Hochberg correction.

Output: Supplementary Table S5 (27 rows × 12 columns).

Usage:
    python sda1_enrichment_analysis.py

Input files (expected in the same directory or adjust INPUT_DIR):
    - L_bin_panta_cdhit_Spyo.csv        : binarized L matrix (genes × phylons)
    - A_bin_panta_cdhit_Spyo.csv        : binarized A matrix (phylons × strains)
    - df_genes_panta_cdhit_final.pickle : full P matrix (8,864 genes × 399 strains)
    - gene_annot_panta_Spyo.csv         : PanTA gene annotations (9,700 genes)

Author: Siddharth M. Chauhan
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import os
from pathlib import Path

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
# ---------------------------------------------------------------------------
# Input/output directories.  Defaults assume this script is run from the
# repository root with the Zenodo matrices unpacked into data/:
#     python scripts/sda1_enrichment_analysis.py
# Override either with an environment variable:
#     SPYO_DATA=/path/to/data SPYO_RESULTS=/path/to/out python scripts/sda1_enrichment_analysis.py
# ---------------------------------------------------------------------------
INPUT_DIR  = Path(os.environ.get("SPYO_DATA", "data"))
OUTPUT_DIR = Path(os.environ.get("SPYO_RESULTS", "results"))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

L_PATH     = INPUT_DIR / "L_bin_panta_cdhit_Spyo.csv"
A_PATH     = INPUT_DIR / "A_bin_panta_cdhit_Spyo.csv"
P_PATH     = INPUT_DIR / "df_genes_panta_cdhit_final.pickle"
ANNOT_PATH = INPUT_DIR / "gene_annot_panta_Spyo.csv"

SDA1_PHYLONS = ["ST28-emm1", "ST36", "ST63-emm77"]

# Functional module assignments for each of the 27 genes, based on
# canonical temperate-prophage modular architecture (see manuscript Fig. 5).
FUNCTIONAL_MODULE_MAP = {
    # ── Integration module ──
    "groups_4947": "Integration",           # Integrase
    "irrE_00560":  "Integration",           # IrrE-family lysogeny regulator

    # ── DNA modification / replication module ──
    "groups_1725": "DNA modification/replication",   # DNA methyltransferase
    "groups_3351": "DNA modification/replication",   # DNA-binding phage protein
    "groups_2028": "DNA modification/replication",   # SNF2-family helicase

    # ── Phage structural backbone (10 genes) ──
    "groups_107":  "Phage structural",      # DUF1366 domain-containing protein
    "groups_365":  "Phage structural",      # Phage protein
    "groups_715":  "Phage structural",      # Phage protein
    "groups_1281": "Phage structural",      # Phage protein
    "groups_1462": "Phage structural",      # Phage protein
    "groups_2024": "Phage structural",      # Phage protein
    "groups_2526": "Phage structural",      # Phage protein
    "groups_3362": "Phage structural",      # Phage protein
    "groups_4935": "Phage structural",      # Phage protein
    "groups_4950": "Phage structural",      # Phage protein

    # ── Lysis module (3 holins + 1 amidase) ──
    "groups_1105": "Lysis",                 # Holin
    "groups_2032": "Lysis",                 # Holin
    "groups_4948": "Lysis",                 # Holin
    "groups_89":   "Lysis",                 # N-acetylmuramoyl-L-alanine amidase

    # ── Paratox cassette (2 genes) ──
    "groups_141":  "Paratox cassette",      # Paratox
    "groups_2030": "Paratox cassette",      # Paratox / Phage protein

    # ── Accessory virulence cargo (6 genes) ──
    "sda1":        "Accessory virulence cargo",  # Streptodornase Sda1
    "groups_1465": "Accessory virulence cargo",  # Hyaluronidase (gp58-like)
    "groups_4930": "Accessory virulence cargo",  # Hyaluronoglucosaminidase
    "groups_3323": "Accessory virulence cargo",  # Secreted protein (unknown)
    "groups_3349": "Accessory virulence cargo",  # Glycerate kinase
    "groups_343":  "Accessory virulence cargo",  # Nucleoside diphosphate kinase
}


def load_data():
    """Load and return all required matrices and annotations."""
    L = pd.read_csv(L_PATH, index_col=0)
    A = pd.read_csv(A_PATH, index_col=0)
    P = pd.read_pickle(P_PATH)
    gene_annot = pd.read_csv(ANNOT_PATH, low_memory=False)
    return L, A, P, gene_annot


def identify_sda1_exclusive_genes(L):
    """
    Identify genes present in ALL Sda1+ phylons and ABSENT from ALL
    other phylons in the binarized L matrix.

    Returns a list of gene names.
    """
    other_phylons = [p for p in L.columns if p not in SDA1_PHYLONS]

    in_all_sda1 = (L[SDA1_PHYLONS] == 1).all(axis=1)
    out_all_other = (L[other_phylons] == 0).all(axis=1)

    genes = L.index[in_all_sda1 & out_all_other].tolist()
    return genes


def get_strain_sets(A):
    """
    Define Sda1+ and non-Sda1 strain sets from the binarized A matrix.

    Returns (sda1_strains: set, non_sda1_strains: set).
    """
    sda1_strains = set()
    for phylon in SDA1_PHYLONS:
        sda1_strains |= set(A.columns[A.loc[phylon] == 1])

    non_sda1_strains = set(A.columns) - sda1_strains
    return sda1_strains, non_sda1_strains


def build_annotation_lookup(gene_annot):
    """
    Build gene -> primary Bakta product annotation lookup.

    PanTA annotations may contain comma-separated synonyms from multiple
    genomes; we retain only the first (primary) annotation.
    """
    lookup = {}
    for _, row in gene_annot.iterrows():
        gene = row["Gene"]
        annot = str(row["Annotation"])
        # Take only the first annotation if comma-separated
        if "," in annot:
            annot = annot.split(",")[0].strip()
        lookup[gene] = annot
    return lookup


def compute_enrichment(genes, P, sda1_strains, non_sda1_strains, annot_lookup):
    """
    For each gene, compute:
      - total carriers across 399 genomes
      - carriage in Sda1+ strains (n=121)
      - carriage in non-Sda1 strains (n=278)
      - log2 odds ratio (Haldane-Anscombe corrected: +0.5 to all cells)
      - Fisher's exact p-value (one-tailed, alternative='greater')
      - functional module assignment

    Returns a DataFrame with one row per gene.
    """
    n_sda1 = len(sda1_strains)
    n_non = len(non_sda1_strains)
    sda1_list = list(sda1_strains)
    non_list = list(non_sda1_strains)

    records = []
    for gene in genes:
        row = P.loc[gene]

        # Carriage counts
        a = int(row[sda1_list].sum())   # present in Sda1+
        b = int(row[non_list].sum())    # present in non-Sda1
        c = n_sda1 - a                  # absent in Sda1+
        d = n_non - b                   # absent in non-Sda1
        total = a + b

        # Fisher's exact test (one-tailed: gene enriched in Sda1+ group)
        _, pval = fisher_exact([[a, b], [c, d]], alternative="greater")

        # Log2 odds ratio with Haldane-Anscombe correction (+0.5 to all cells)
        a_h, b_h, c_h, d_h = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        log2_or = np.log2((a_h * d_h) / (b_h * c_h))

        records.append({
            "Gene cluster": gene,
            "Bakta product annotation": annot_lookup.get(gene, "Unknown"),
            "Functional module": FUNCTIONAL_MODULE_MAP.get(gene, "Unknown"),
            "Total carriers (n=399)": total,
            f"Carriage in Sda1+ (n={n_sda1})": a,
            f"Carriage in non-Sda1 (n={n_non})": b,
            "Log2 odds ratio (HA-corrected)": round(log2_or, 2),
            "Fisher p-value (one-tailed)": pval,
        })

    df = pd.DataFrame(records)
    return df


def apply_bh_correction(df):
    """Apply Benjamini-Hochberg FDR correction and add enrichment call."""
    _, qvals, _, _ = multipletests(
        df["Fisher p-value (one-tailed)"].values,
        method="fdr_bh",
    )
    df["BH q-value"] = qvals
    df["Enriched (FDR<0.05)"] = ["Yes" if q < 0.05 else "No" for q in qvals]
    return df


def sort_by_module(df):
    """Sort table by functional module in canonical prophage order."""
    module_order = [
        "Integration",
        "DNA modification/replication",
        "Phage structural",
        "Lysis",
        "Paratox cassette",
        "Accessory virulence cargo",
    ]
    order_map = {m: i for i, m in enumerate(module_order)}
    df["_sort_key"] = df["Functional module"].map(order_map)
    df = df.sort_values("_sort_key").drop(columns="_sort_key").reset_index(drop=True)
    return df


def print_summary(df, n_sda1, n_non):
    """Print summary statistics matching manuscript claims."""
    sda1_col = f"Carriage in Sda1+ (n={n_sda1})"
    non_col = f"Carriage in non-Sda1 (n={n_non})"

    pct_sda1 = df[sda1_col] / n_sda1 * 100
    pct_non = df[non_col] / n_non * 100

    print("\n" + "=" * 80)
    print("SUMMARY STATISTICS")
    print("=" * 80)
    print(f"Total Sda1-exclusive genes identified: {len(df)}")
    print(f"All enriched at FDR < 0.05? {(df['Enriched (FDR<0.05)'] == 'Yes').all()}")
    print(f"")
    print(f"Carriage in Sda1+ strains (n={n_sda1}):")
    print(f"  Median: {pct_sda1.median():.1f}%")
    print(f"  Range:  {pct_sda1.min():.1f}% – {pct_sda1.max():.1f}%")
    print(f"")
    print(f"Carriage in non-Sda1 strains (n={n_non}):")
    print(f"  Median: {pct_non.median():.1f}%")
    print(f"  Range:  {pct_non.min():.1f}% – {pct_non.max():.1f}%")
    print(f"")
    print(f"BH-corrected q-values:")
    print(f"  Min:    {df['BH q-value'].min():.2e}")
    print(f"  Median: {df['BH q-value'].median():.2e}")
    print(f"  Max:    {df['BH q-value'].max():.2e}")
    print("=" * 80)


def main():
    # ── Load ──
    print("Loading data...")
    L, A, P, gene_annot = load_data()
    print(f"  L matrix: {L.shape[0]} genes × {L.shape[1]} phylons")
    print(f"  A matrix: {A.shape[0]} phylons × {A.shape[1]} strains")
    print(f"  P matrix: {P.shape[0]} genes × {P.shape[1]} strains")

    # ── Identify Sda1-exclusive genes ──
    genes = identify_sda1_exclusive_genes(L)
    print(f"\nIdentified {len(genes)} genes in all Sda1+ phylons and absent from all others.")

    # ── Define strain sets ──
    sda1_strains, non_sda1_strains = get_strain_sets(A)
    n_sda1, n_non = len(sda1_strains), len(non_sda1_strains)
    print(f"Sda1+ strains: {n_sda1}")
    print(f"Non-Sda1 strains: {n_non}")

    # ── Build annotation lookup ──
    annot_lookup = build_annotation_lookup(gene_annot)

    # ── Compute enrichment ──
    print("\nComputing Fisher's exact test for each gene...")
    df = compute_enrichment(genes, P, sda1_strains, non_sda1_strains, annot_lookup)

    # ── BH correction ──
    df = apply_bh_correction(df)

    # ── Sort by module ──
    df = sort_by_module(df)

    # ── Print summary ──
    print_summary(df, n_sda1, n_non)

    # ── Save ──
    out_path = OUTPUT_DIR / "TableS7_Sda1_prophage_core_enrichment.csv"
    df.to_csv(out_path, index=False)
    print(f"\nTable saved to: {out_path}")
    print(f"Shape: {df.shape}")

    return df


if __name__ == "__main__":
    main()
