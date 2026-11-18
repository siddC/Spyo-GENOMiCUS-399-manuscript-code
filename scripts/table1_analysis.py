#!/usr/bin/env python3
"""
Table 1: Per-phylon virulence and ecological characterization across the
seven lineage-specific S. pyogenes phylons of Spyo-GENOMiCUS-400.

METHODOLOGY:
  Isolation sites are classified from the standardized source_cat metadata
  column with two additional rules:
    (1) ST63-emm77 'Respiratory' entries are reclassified as
        'Throat/Pharyngeal' (Polish upper-respiratory pharyngitis
        surveillance study, refs 51–52).
    (2) Raw metadata fields (isolation_source, disease,
        isolation_comments) are scanned for explicit invasive-disease
        keywords: invasive, STSS, necrotizing fasciitis, bacteremia,
        sepsis, meningitis, cellulitis. Strains matching these keywords
        but lacking a specific anatomical-site classification are
        assigned to 'Invasive (unspecified)' and included in both the
        numerator and denominator of the invasive isolation rate.

  Invasive isolation rate = percentage of strains from invasive sites
  (blood, sterile site, CSF, bone/joint, soft tissue) OR with explicit
  invasive-disease labels among strains with available source metadata.
  Strains with source_cat = 'Not Specified' and no invasive keywords,
  as well as non-anatomical categories (Genital, Urine, Ear), are
  excluded from the denominator.

  Virulence content is derived from the binarized L matrix: for each
  phylon, genes present in that phylon and absent from all others
  (exclusive genes) are categorized by functional annotation.

Input files (expected in working directory):
    - L_bin_panta_cdhit_Spyo.csv
    - A_bin_panta_cdhit_Spyo.csv
    - df_genes_panta_cdhit_final.pickle
    - car_metadata_Sypo.csv
    - gene_annot_panta_Spyo.csv

Output:
    Table1_per_phylon_virulence_ecology.csv

Author: Siddharth M. Chauhan
"""

import pandas as pd
import numpy as np
import os
from pathlib import Path

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
# ---------------------------------------------------------------------------
# Input/output directories.  Defaults assume this script is run from the
# repository root with the Zenodo matrices unpacked into data/:
#     python scripts/table1_analysis.py
# Override either with an environment variable:
#     SPYO_DATA=/path/to/data SPYO_RESULTS=/path/to/out python scripts/table1_analysis.py
# ---------------------------------------------------------------------------
INPUT_DIR  = Path(os.environ.get("SPYO_DATA", "data"))
OUTPUT_DIR = Path(os.environ.get("SPYO_RESULTS", "results"))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

L_PATH     = INPUT_DIR / "L_bin_panta_cdhit_Spyo.csv"
A_PATH     = INPUT_DIR / "A_bin_panta_cdhit_Spyo.csv"
P_PATH     = INPUT_DIR / "df_genes_panta_cdhit_final.pickle"
META_PATH  = INPUT_DIR / "car_metadata_Sypo.csv"
ANNOT_PATH = INPUT_DIR / "gene_annot_panta_Spyo.csv"

LINEAGE_PHYLONS = [
    "ST28-emm1", "ST36", "ST63-emm77", "ST52-emm28",
    "ST39-emm4", "emm3/6/5", "ST101-emm89",
]

SDA1_PHYLONS = ["ST28-emm1", "ST36", "ST63-emm77"]

INVASIVE_SITES = {
    "Blood", "Sterile site", "CSF", "Bone/Joint",
    "Soft Tissue", "Invasive (unspecified)",
}

# Keywords for raw-metadata invasive-disease scanning
INVASIVE_KEYWORDS = (
    "invasive|stss|necrotizing fasciitis|bacteremia|sepsis|meningitis|cellulitis"
)

# Direct source_cat → standardized site mapping
DIRECT_SITE_MAP = {
    "Blood":              "Blood",
    "Throat/Pharyngeal":  "Throat/Pharyngeal",
    "Skin":               "Skin",
    "Soft Tissue":        "Soft Tissue",
    "CSF":                "CSF",
    "Bone/Joint":         "Bone/Joint",
    "Respiratory":        "Respiratory",
    "Wound":              "Soft Tissue",
    "Abscess/Pus":        "Soft Tissue",
}

# ── Per-phylon virulence content (curated from exclusive-gene analysis) ──
VIRULENCE = {
    "ST28-emm1": {
        "Sda1": "+",
        "Other superantigens": "SpeA",
        "DNases (other)": "Sda3",
        "Adhesins/pili": "FCT-2 pilus, Scl1, Scl2, FbpA",
        "Immune evasion": "Sic",
        "Other MGEs": "\u2014",
        "Metabolic/competition": "\u2014",
    },
    "ST36": {
        "Sda1": "+",
        "Other superantigens": "SpeH, SpeI",
        "DNases (other)": "\u2014",
        "Adhesins/pili": "PrtF1/SfbI, distinct pilus locus",
        "Immune evasion": "\u2014",
        "Other MGEs": "\u2014",
        "Metabolic/competition": "\u2014",
    },
    "ST63-emm77": {
        "Sda1": "+",
        "Other superantigens": "\u2014",
        "DNases (other)": "\u2014",
        "Adhesins/pili": "\u2014",
        "Immune evasion": "\u2014",
        "Other MGEs": "2nd prophage (HylA, IsdE); ICE",
        "Metabolic/competition": "\u2014",
    },
    "ST52-emm28": {
        "Sda1": "\u2014",
        "Other superantigens": "SpeC paralog (lineage-specific)",
        "DNases (other)": "\u2014",
        "Adhesins/pili": "FbaA",
        "Immune evasion": "\u2014",
        "Other MGEs": "\u2014",
        "Metabolic/competition": "3 TA systems (MazF, RelE/RelB, HicB)",
    },
    "ST39-emm4": {
        "Sda1": "\u2014",
        "Other superantigens": "SpeC (mobilon)",
        "DNases (other)": "\u2014",
        "Adhesins/pili": "FCT-3 pili",
        "Immune evasion": "\u2014",
        "Other MGEs": "SpeC prophage (mobilon)",
        "Metabolic/competition": "DhaKLM (glycerol/DHA operon)",
    },
    "emm3/6/5": {
        "Sda1": "\u2014",
        "Other superantigens": "\u2014",
        "DNases (other)": "\u2014",
        "Adhesins/pili": "\u2014",
        "Immune evasion": "\u2014",
        "Other MGEs": "HylA/IsdE prophage",
        "Metabolic/competition": "\u2014",
    },
    "ST101-emm89": {
        "Sda1": "\u2014",
        "Other superantigens": "\u2014",
        "DNases (other)": "\u2014",
        "Adhesins/pili": "\u2014",
        "Immune evasion": "\u2014",
        "Other MGEs": "\u2014",
        "Metabolic/competition": "LanM lantibiotic system",
    },
}


# ══════════════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════════════
def load_data():
    L = pd.read_csv(L_PATH, index_col=0)
    A = pd.read_csv(A_PATH, index_col=0)
    P = pd.read_pickle(P_PATH)
    meta = pd.read_csv(META_PATH, low_memory=False)
    gene_annot = pd.read_csv(ANNOT_PATH, low_memory=False)
    return L, A, P, meta, gene_annot


# ══════════════════════════════════════════════════════════════════════════════
# ISOLATION SITE CLASSIFICATION (Option 2: with raw-text keyword scanning)
# ══════════════════════════════════════════════════════════════════════════════
def classify_isolation_sites(meta, A):
    """
    Classify each strain's isolation site using source_cat plus
    raw-text keyword scanning of isolation_source, disease, and
    isolation_comments fields.
    """
    st63_strains = set(A.columns[A.loc["ST63-emm77"] == 1])

    def classify_row(row):
        cat = row["source_cat"]
        iso = str(row["isolation_source"]).lower()
        gid = str(row["genome_id"])

        # ST63-emm77 respiratory → throat reclassification
        if gid in st63_strains and cat == "Respiratory":
            return "Throat/Pharyngeal"

        # Direct mapping
        if cat in DIRECT_SITE_MAP:
            return DIRECT_SITE_MAP[cat]

        # Clinical (Unspecified): check for sterile site or invasive keywords
        if cat == "Clinical (Unspecified)":
            if "sterile" in iso:
                return "Sterile site"
            # Fall through to keyword scan below

        return None  # will be handled by keyword scan

    meta["site"] = meta.apply(classify_row, axis=1)

    # Raw-text keyword scanning for remaining unclassified strains
    for col in ["isolation_source", "disease", "isolation_comments"]:
        if col in meta.columns:
            mask = (
                meta[col].astype(str).str.lower().str.contains(
                    INVASIVE_KEYWORDS, na=False
                )
                & meta["site"].isna()
            )
            meta.loc[mask, "site"] = "Invasive (unspecified)"

    meta["is_invasive"] = meta["site"].isin(INVASIVE_SITES)
    return meta


# ══════════════════════════════════════════════════════════════════════════════
# COMPUTE PER-PHYLON ECOLOGY
# ══════════════════════════════════════════════════════════════════════════════
def compute_ecology(phylon, A, meta):
    """
    For one phylon, compute:
      - n_strains: total members
      - n_with_metadata: members with a classified isolation site
      - site breakdown: count per site category
      - dominant_site: most common site and its percentage
      - invasive_rate: % invasive among strains with metadata
    """
    strains = list(A.columns[A.loc[phylon] == 1])
    sub = meta[meta["genome_id"].astype(str).isin(strains)]
    with_meta = sub[sub["site"].notna()]
    n_total = len(strains)
    n_meta = len(with_meta)

    result = {
        "n_strains": n_total,
        "n_with_metadata": n_meta,
    }

    if n_meta == 0:
        result["dominant_site"] = "N/A"
        result["invasive_rate_pct"] = "N/A"
        result["invasive_n"] = "N/A"
        result["site_breakdown"] = {}
        return result

    # Site breakdown
    site_counts = with_meta["site"].value_counts().to_dict()
    result["site_breakdown"] = site_counts

    # Dominant site
    dom = max(site_counts, key=site_counts.get)
    dom_pct = round(site_counts[dom] / n_meta * 100)
    result["dominant_site"] = f"{dom} ({dom_pct}%)"

    # Invasive rate
    n_inv = int(with_meta["is_invasive"].sum())
    inv_pct = round(n_inv / n_meta * 100)
    result["invasive_rate_pct"] = f"{inv_pct}%"
    result["invasive_n"] = f"{n_inv}/{n_meta}"

    return result


# ══════════════════════════════════════════════════════════════════════════════
# COUNT EXCLUSIVE GENES
# ══════════════════════════════════════════════════════════════════════════════
def count_exclusive(phylon, L):
    others = [q for q in L.columns if q != phylon]
    return int(((L[phylon] == 1) & (L[others] == 0).all(axis=1)).sum())


# ══════════════════════════════════════════════════════════════════════════════
# BUILD TABLE 1
# ══════════════════════════════════════════════════════════════════════════════
def build_table1(L, A, meta):
    rows = []
    for p in LINEAGE_PHYLONS:
        eco = compute_ecology(p, A, meta)
        n_excl = count_exclusive(p, L)
        vir = VIRULENCE[p]

        rows.append({
            "Phylon": p,
            "n_strains": eco["n_strains"],
            "n_exclusive_genes": n_excl,
            "Sda1": vir["Sda1"],
            "Other superantigens": vir["Other superantigens"],
            "DNases (other)": vir["DNases (other)"],
            "Adhesins/pili": vir["Adhesins/pili"],
            "Immune evasion": vir["Immune evasion"],
            "Other MGEs": vir["Other MGEs"],
            "Metabolic/competition": vir["Metabolic/competition"],
            "Dominant isolation site (%)": eco["dominant_site"],
            "Invasive isolation rate (%)": eco["invasive_rate_pct"],
            "Invasive (n/N)": eco["invasive_n"],
        })

    return pd.DataFrame(rows)


# ══════════════════════════════════════════════════════════════════════════════
# PRINT DETAILED REPORT
# ══════════════════════════════════════════════════════════════════════════════
def print_report(table1, A, meta):
    print(f"\n{'=' * 90}")
    print("TABLE 1: Per-phylon virulence and ecological characterization")
    print(f"{'=' * 90}")

    for _, r in table1.iterrows():
        p = r["Phylon"]
        print(f"\n{'─' * 70}")
        print(f"{p} (n={r['n_strains']}, {r['n_exclusive_genes']} exclusive genes)")
        print(f"{'─' * 70}")
        print(f"  Sda1:                {r['Sda1']}")
        print(f"  Other superantigens: {r['Other superantigens']}")
        print(f"  DNases (other):      {r['DNases (other)']}")
        print(f"  Adhesins/pili:       {r['Adhesins/pili']}")
        print(f"  Immune evasion:      {r['Immune evasion']}")
        print(f"  Other MGEs:          {r['Other MGEs']}")
        print(f"  Metabolic/compet.:   {r['Metabolic/competition']}")
        print(f"  Dominant site:       {r['Dominant isolation site (%)']}")
        print(f"  Invasive rate:       {r['Invasive isolation rate (%)']} "
              f"({r['Invasive (n/N)']})")

        # Detailed site breakdown
        strains = list(A.columns[A.loc[p] == 1])
        sub = meta[meta["genome_id"].astype(str).isin(strains)]
        with_meta = sub[sub["site"].notna()]
        if len(with_meta) > 0:
            print(f"  Site breakdown (n={len(with_meta)}):")
            for site, count in with_meta["site"].value_counts().items():
                pct = count / len(with_meta) * 100
                inv_mark = " [invasive]" if site in INVASIVE_SITES else ""
                print(f"    {site}: {count} ({pct:.0f}%){inv_mark}")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    print("Loading data...")
    L, A, P, meta, gene_annot = load_data()
    print(f"  {L.shape[0]} accessory genes, {A.shape[1]} strains, "
          f"{len(meta)} metadata rows")

    print("Classifying isolation sites (Option 2: with keyword scanning)...")
    meta = classify_isolation_sites(meta, A)

    n_classified = meta["site"].notna().sum()
    n_invasive_kw = (meta["site"] == "Invasive (unspecified)").sum()
    print(f"  {n_classified}/{len(meta)} strains with classified site")
    print(f"  {n_invasive_kw} classified via invasive-keyword scanning")

    print("Building Table 1...")
    table1 = build_table1(L, A, meta)

    print_report(table1, A, meta)

    # Save
    out_path = OUTPUT_DIR / "Table1_per_phylon_virulence_ecology.csv"
    table1.to_csv(out_path, index=False)
    print(f"\n\nSaved: {out_path}")
    print(f"Shape: {table1.shape}")


if __name__ == "__main__":
    main()
