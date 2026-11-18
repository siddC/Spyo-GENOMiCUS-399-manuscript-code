#!/usr/bin/env python3
"""
Reproduce Supplementary Table S9: speC/mf2 synteny annotations.

For each of three representative genomes (one per speC paralog pair),
extracts a ±10-feature window around the speC locus, assigns functional
categories, and computes display coordinates anchored at the speC start
codon. The GCF_041430965.1 prophage is integrated in the opposite
orientation relative to the chromosome; its coordinates are
reverse-complemented for display consistency.

Input files (expected in working directory or adjust paths):
    1314_985.gff3           (ST52-emm28, speC + mf2)
    1314_3072.gff3          (ST39-emm4, speC_01049 + mf2_01051)
    GCF_041430965_1.gff3    (emm-patternE-other, speC_01050 + mf2_01052)

Output:
    TableS9_speC_mf2_synteny_annotations.csv  (63 rows × 14 columns)

Author: Siddharth M. Chauhan
"""

import re
import pandas as pd
from pathlib import Path

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
OUTPUT_DIR = Path(".")

REPS = [
    {
        "genome_id":     "1314.985",
        "phylon":        "ST52-emm28",
        "paralog_pair":  "speC + mf2",
        "gff":           "1314_985.gff3",
        "speC_lt":       "GGAFNI_00883",
        "mf2_lt":        "GGAFNI_00884",
        "reverse_comp":  False,
    },
    {
        "genome_id":     "1314.3072",
        "phylon":        "ST39-emm4",
        "paralog_pair":  "speC_01049 + mf2_01051",
        "gff":           "1314_3072.gff3",
        "speC_lt":       "FEFGFE_00607",
        "mf2_lt":        "FEFGFE_00608",
        "reverse_comp":  False,
    },
    {
        "genome_id":     "GCF_041430965.1",
        "phylon":        "emm-patternE-other",
        "paralog_pair":  "speC_01050 + mf2_01052",
        "gff":           "GCF_041430965_1.gff3",
        "speC_lt":       "LABOOI_01147",
        "mf2_lt":        "LABOOI_01146",
        "reverse_comp":  True,
    },
]

# Number of features to include on each side of speC
FLANK_FEATURES = 10


# ══════════════════════════════════════════════════════════════════════════════
# GFF3 PARSER
# ══════════════════════════════════════════════════════════════════════════════
def parse_gff_cds(path):
    """Parse CDS features from a GFF3 file (stops at ##FASTA)."""
    feats = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "CDS":
                continue
            attrs = dict(re.findall(r"([^=;]+)=([^;]+)", parts[8]))
            feats.append({
                "contig":    parts[0],
                "start":     int(parts[3]),
                "end":       int(parts[4]),
                "strand":    parts[6],
                "locus_tag": attrs.get("locus_tag", ""),
                "gene":      attrs.get("gene", "\u2014"),
                "product":   attrs.get("product", ""),
            })
    return pd.DataFrame(feats)


# ══════════════════════════════════════════════════════════════════════════════
# FUNCTIONAL CATEGORY ASSIGNMENT
# ══════════════════════════════════════════════════════════════════════════════
def assign_functional_category(row, speC_lt, mf2_lt):
    """Assign a functional category based on locus tag and product."""
    lt = row["locus_tag"]
    prod = row["product"].lower()

    if lt == speC_lt:
        return "speC"
    if lt == mf2_lt:
        return "mf2"
    if "paratox" in prod:
        return "Paratox"
    if "hyaluron" in prod:
        return "Hyaluronidase / HylA"
    if any(k in prod for k in ("holin", "lysin", "amidase", "endopeptidase")):
        return "Phage lysis"
    if "phage" in prod:
        return "Phage structural"
    if any(k in prod for k in (
        "tail", "capsid", "portal", "terminase", "head",
        "morphogenesis", "duf1366", "vrr-nuc", "tape measure",
    )):
        return "Phage structural"
    if "primase" in prod and "subunit" not in prod:
        return "Phage structural"
    if "dna-directed dna polymerase" in prod:
        return "Phage structural"
    return "Chromosomal"


# ══════════════════════════════════════════════════════════════════════════════
# BUILD TABLE FOR ONE REPRESENTATIVE GENOME
# ══════════════════════════════════════════════════════════════════════════════
def build_genome_table(rep):
    """
    Parse the GFF3, extract a ±FLANK_FEATURES window around speC,
    assign functional categories, and compute display coordinates.
    """
    df = parse_gff_cds(rep["gff"])
    df = df.sort_values("start").reset_index(drop=True)

    # Find speC index
    speC_mask = df["locus_tag"] == rep["speC_lt"]
    speC_idx = df[speC_mask].index[0]
    speC_row = df.loc[speC_idx]

    # Extract window
    lo = max(0, speC_idx - FLANK_FEATURES)
    hi = min(len(df) - 1, speC_idx + FLANK_FEATURES)
    region = df.iloc[lo : hi + 1].copy()

    # Assign functional categories
    region["functional_category"] = region.apply(
        lambda r: assign_functional_category(r, rep["speC_lt"], rep["mf2_lt"]),
        axis=1,
    )

    # Compute display coordinates
    rows = []
    if rep["reverse_comp"]:
        # Anchor at speC end coordinate (for reverse-complement display)
        speC_end = speC_row["end"]
        for _, row in region.iterrows():
            d_start = (speC_end - row["end"]) / 1000
            d_end = (speC_end - row["start"]) / 1000
            disp_strand = "+" if row["strand"] == "-" else "-"
            rows.append({
                "genome_id":                        rep["genome_id"],
                "phylon":                           rep["phylon"],
                "paralog_pair":                     rep["paralog_pair"],
                "locus_tag":                        row["locus_tag"],
                "genomic_start":                    row["start"],
                "genomic_end":                      row["end"],
                "genomic_strand":                   row["strand"],
                "length_bp":                        row["end"] - row["start"] + 1,
                "display_position_kb":              f"{d_start:.3g} to {d_end:.3g}",
                "display_strand":                   disp_strand,
                "gene":                             row["gene"],
                "product":                          row["product"],
                "functional_category":              row["functional_category"],
                "reverse_complemented_for_display": "Yes",
            })
    else:
        # Anchor at speC start coordinate
        anchor = speC_row["start"]
        for _, row in region.iterrows():
            d_start = (row["start"] - anchor) / 1000
            d_end = (row["end"] - anchor) / 1000
            rows.append({
                "genome_id":                        rep["genome_id"],
                "phylon":                           rep["phylon"],
                "paralog_pair":                     rep["paralog_pair"],
                "locus_tag":                        row["locus_tag"],
                "genomic_start":                    row["start"],
                "genomic_end":                      row["end"],
                "genomic_strand":                   row["strand"],
                "length_bp":                        row["end"] - row["start"] + 1,
                "display_position_kb":              f"{d_start:.3g} to {d_end:.3g}",
                "display_strand":                   row["strand"],
                "gene":                             row["gene"],
                "product":                          row["product"],
                "functional_category":              row["functional_category"],
                "reverse_complemented_for_display": "No",
            })

    return pd.DataFrame(rows)


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    tables = []
    for rep in REPS:
        print(f"Processing {rep['genome_id']} ({rep['phylon']}, {rep['paralog_pair']})...")
        t = build_genome_table(rep)
        tables.append(t)
        print(f"  {len(t)} features extracted")

    result = pd.concat(tables, ignore_index=True)

    # Verify chromosomal integration sites
    print(f"\nChromosomal integration sites (attR flanking genes):")
    for rep in REPS:
        sub = result[result["genome_id"] == rep["genome_id"]]
        chromo = sub[sub["functional_category"] == "Chromosomal"]
        flank_genes = chromo["product"].head(4).tolist()
        print(f"  {rep['genome_id']}: {', '.join(flank_genes)}")

    # Save
    out_path = OUTPUT_DIR / "TableS9_speC_mf2_synteny_annotations.csv"
    result.to_csv(out_path, index=False)
    print(f"\nTable saved: {out_path}")
    print(f"Shape: {result.shape}")


if __name__ == "__main__":
    main()
