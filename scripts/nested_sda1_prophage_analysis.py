#!/usr/bin/env python3
"""
Analysis of the nested Sda1-encoding prophage in ST28-emm1 strain 1314.4589.

This script performs three analyses on the double-prophage region:

1. SEQUENCE COMPARISON: Extracts the canonical sda1 (BENLNG_01457) and
   the paralog sda1_00717 (BENLNG_01492) coding sequences, compares them
   at nucleotide and protein level (BLASTn + direct translation).

2. PROPHAGE-WIDE BLASTp: Identifies all CDS in Prophage 1 (canonical,
   ~68 kb) and Prophage 2 (nested, ~28 kb), extracts their protein
   sequences, and runs BLASTp (P2 as query, P1 as database) to quantify
   the duplication relationship.

3. SYNTENY FIGURE: Generates Supplementary Figure S6 showing the nested
   prophage architecture with color-coded functional modules.

Dependencies:
    - biopython (pip install biopython)
    - ncbi-blast+ (apt install ncbi-blast+)
    - matplotlib, pandas

Input:
    1314_4589.gff3  (Bakta-annotated GFF3 with embedded FASTA)

Outputs:
    figS6_double_sda1_prophage.pdf
    figS6_double_sda1_prophage.png
    sda1_paralog_blastn.txt          (BLASTn: canonical vs paralog)
    prophage2_vs_prophage1_blastp.txt (BLASTp: all P2 genes vs P1)
    prophage2_vs_prophage1_blastp_summary.csv

Author: Siddharth M. Chauhan
"""

import re
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrow
from Bio.Seq import Seq

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
GFF_PATH   = Path("1314_4589.gff3")
OUTPUT_DIR = Path(".")

# Genomic window for the synteny figure (bp on contig_1)
WINDOW_LO = 1394000
WINDOW_HI = 1480000

# Prophage boundary locus tags
P1_ATT_L      = "BENLNG_01456"   # Paratox at 5' end of Prophage 1
P1_ATT_R      = "BENLNG_01545"   # Integrase at 3' end of Prophage 1
P2_START      = "BENLNG_01470"   # First gene of nested Prophage 2
P2_END        = "BENLNG_01504"   # Last gene of nested Prophage 2
SDA1_CANON_LT = "BENLNG_01457"   # Canonical sda1
SDA1_PARA_LT  = "BENLNG_01492"   # Paralog sda1_00717
P2_INTEGRASE  = "BENLNG_01490"   # Prophage 2 integrase

# Locus-tag mapping for the canonical 27-gene Sda1+ set in 1314.4589
CANONICAL_LTS = {
    "sda1": "BENLNG_01457", "irrE_00560": "BENLNG_01543",
    "groups_89": "BENLNG_01458", "groups_107": "BENLNG_01462",
    "groups_141": "BENLNG_01456", "groups_365": "BENLNG_01537",
    "groups_715": "BENLNG_01523", "groups_1105": "BENLNG_01461",
    "groups_1281": "BENLNG_01536", "groups_1462": "BENLNG_01524",
    "groups_1465": "BENLNG_01464", "groups_1725": "BENLNG_01522",
    "groups_2024": "BENLNG_01544", "groups_2028": "BENLNG_01525",
    "groups_2030": "BENLNG_01521", "groups_2032": "BENLNG_01460",
    "groups_2526": "BENLNG_01541", "groups_3349": "BENLNG_01539",
    "groups_3351": "BENLNG_01538", "groups_3362": "BENLNG_01498",
    "groups_4930": "BENLNG_01501", "groups_4935": "BENLNG_01500",
    "groups_4947": "BENLNG_01545", "groups_4948": "BENLNG_01540",
    "groups_4950": "BENLNG_01530",
}

# 27-gene functional module assignments
MODULE = {
    "groups_4947": "integration",   "irrE_00560":  "integration",
    "groups_1725": "dna_replication","groups_3351": "dna_replication",
    "groups_2028": "dna_replication",
    "groups_107":  "structural",    "groups_365":  "structural",
    "groups_715":  "structural",    "groups_1281": "structural",
    "groups_1462": "structural",    "groups_2024": "structural",
    "groups_2526": "structural",    "groups_3362": "structural",
    "groups_4935": "structural",    "groups_4950": "structural",
    "groups_1105": "lysis",         "groups_2032": "lysis",
    "groups_4948": "lysis",         "groups_89":   "lysis",
    "groups_141":  "paratox",       "groups_2030": "paratox",
    "sda1":        "cargo_sda1",
    "groups_1465": "cargo_other",   "groups_4930": "cargo_other",
    "groups_3323": "cargo_other",   "groups_3349": "cargo_other",
    "groups_343":  "cargo_other",
}

COLORS = {
    "integration":        "#9b59b6",
    "dna_replication":    "#3498db",
    "structural":         "#2ecc71",
    "paratox":            "#e67e22",
    "cargo_sda1":         "#e74c3c",
    "cargo_sda1_paralog": "#c0392b",
    "cargo_other":        "#f39c12",
    "lysis":              "#1abc9c",
    "other_phage":        "#bdc3c7",
    "host":               "#7f8c8d",
}

MODULE_LABELS = {
    "integration":        "Integration / lysogeny",
    "dna_replication":    "DNA modification / replication",
    "structural":         "Phage structural backbone",
    "paratox":            "Paratox cassette",
    "cargo_sda1":         "sda1 (cargo)",
    "cargo_sda1_paralog": "sda1 paralog (sda1_00717)",
    "cargo_other":        "Other cargo",
    "lysis":              "Lysis module (holins, amidase)",
    "other_phage":        "Other phage genes (not in 27-gene set)",
    "host":               "Host chromosome",
}


# ══════════════════════════════════════════════════════════════════════════════
# PARSERS
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
                "product":   attrs.get("product", ""),
            })
    return pd.DataFrame(feats)


def extract_fasta(path):
    """Extract embedded FASTA sequences from a GFF3 file."""
    in_fasta = False
    seq_lines, seq_id, seqs = [], None, {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("##FASTA"):
                in_fasta = True
                continue
            if in_fasta:
                if line.startswith(">"):
                    if seq_id:
                        seqs[seq_id] = "".join(seq_lines)
                    seq_id = line[1:].strip().split()[0]
                    seq_lines = []
                else:
                    seq_lines.append(line.strip())
        if seq_id:
            seqs[seq_id] = "".join(seq_lines)
    return seqs


def revcomp(seq):
    """Reverse complement a DNA sequence."""
    comp = str.maketrans("ATCGatcg", "TAGCtagc")
    return seq.translate(comp)[::-1]


def get_cds_nt(genome_seq, start, end, strand):
    """Extract CDS nucleotide sequence (GFF3 is 1-based, inclusive)."""
    nt = genome_seq[start - 1 : end]
    if strand == "-":
        nt = revcomp(nt)
    return nt


def translate_cds(nt_seq):
    """Translate a nucleotide CDS to protein (strip trailing stop)."""
    aa = str(Seq(nt_seq).translate())
    return aa.rstrip("*")


# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1: sda1 canonical vs paralog comparison
# ══════════════════════════════════════════════════════════════════════════════
def run_sda1_comparison(cds_df, genome_seq):
    """Compare canonical sda1 and sda1_00717 paralog at nt and protein level."""
    print("\n" + "=" * 72)
    print("ANALYSIS 1: sda1 canonical vs sda1_00717 paralog")
    print("=" * 72)

    canon = cds_df[cds_df["locus_tag"] == SDA1_CANON_LT].iloc[0]
    para  = cds_df[cds_df["locus_tag"] == SDA1_PARA_LT].iloc[0]

    canon_nt = get_cds_nt(genome_seq, canon["start"], canon["end"], canon["strand"])
    para_nt  = get_cds_nt(genome_seq, para["start"],  para["end"],  para["strand"])

    canon_aa = translate_cds(canon_nt)
    para_aa  = translate_cds(para_nt)

    # Direct nucleotide comparison
    nt_matches = sum(1 for a, b in zip(canon_nt, para_nt) if a == b)
    nt_total   = min(len(canon_nt), len(para_nt))
    nt_pident  = nt_matches / nt_total * 100

    # Protein comparison
    aa_matches = sum(1 for a, b in zip(canon_aa, para_aa) if a == b)
    aa_total   = min(len(canon_aa), len(para_aa))
    aa_pident  = aa_matches / aa_total * 100

    print(f"\n  Canonical sda1 ({SDA1_CANON_LT}):")
    print(f"    Position: {canon['start']}-{canon['end']} ({canon['strand']})")
    print(f"    CDS: {len(canon_nt)} bp -> {len(canon_aa)} aa")
    print(f"\n  Paralog sda1_00717 ({SDA1_PARA_LT}):")
    print(f"    Position: {para['start']}-{para['end']} ({para['strand']})")
    print(f"    CDS: {len(para_nt)} bp -> {len(para_aa)} aa")
    print(f"\n  Nucleotide identity: {nt_matches}/{nt_total} ({nt_pident:.1f}%)")
    print(f"  Protein identity:   {aa_matches}/{aa_total} ({aa_pident:.1f}%)")
    print(f"  Sequences identical? NT={canon_nt == para_nt}, AA={canon_aa == para_aa}")

    # Write FASTA files and run BLASTn for formal confirmation
    with tempfile.TemporaryDirectory() as tmpdir:
        canon_fa = Path(tmpdir) / "canon.fasta"
        para_fa  = Path(tmpdir) / "para.fasta"
        db_name  = Path(tmpdir) / "canon_db"
        out_file = OUTPUT_DIR / "sda1_paralog_blastn.txt"

        canon_fa.write_text(f">canonical_sda1_{SDA1_CANON_LT}\n{canon_nt}\n")
        para_fa.write_text(f">paralog_sda1_00717_{SDA1_PARA_LT}\n{para_nt}\n")

        subprocess.run(
            ["makeblastdb", "-in", str(canon_fa), "-dbtype", "nucl",
             "-out", str(db_name)],
            capture_output=True,
        )
        result = subprocess.run(
            ["blastn", "-query", str(para_fa), "-db", str(db_name),
             "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                        "qstart qend sstart send evalue bitscore",
             "-out", str(out_file)],
            capture_output=True,
        )
        blastn = pd.read_csv(
            out_file, sep="\t", header=None,
            names=["qseqid", "sseqid", "pident", "length", "mismatch",
                   "gapopen", "qstart", "qend", "sstart", "send",
                   "evalue", "bitscore"],
        )
        if len(blastn) > 0:
            r = blastn.iloc[0]
            print(f"\n  BLASTn confirmation:")
            print(f"    {r['pident']:.1f}% identity, {r['length']} bp aligned, "
                  f"{r['mismatch']} mismatches, E={r['evalue']}")

    return nt_pident, aa_pident


# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 2: Prophage-wide BLASTp
# ══════════════════════════════════════════════════════════════════════════════
def run_prophage_blastp(cds_df, genome_seq):
    """
    Extract all proteins from Prophage 1 (excluding the nested P2 region)
    and Prophage 2, then BLASTp P2 queries against P1 database.
    """
    print("\n" + "=" * 72)
    print("ANALYSIS 2: Prophage-wide BLASTp (P2 query vs P1 database)")
    print("=" * 72)

    # Define boundaries
    p1_lo = cds_df[cds_df["locus_tag"] == P1_ATT_L]["start"].iloc[0]
    p1_hi = cds_df[cds_df["locus_tag"] == P1_ATT_R]["end"].iloc[0]
    p2_lo = cds_df[cds_df["locus_tag"] == P2_START]["start"].iloc[0]
    p2_hi = cds_df[cds_df["locus_tag"] == P2_END]["end"].iloc[0]

    print(f"\n  Prophage 1: {p1_lo:,}-{p1_hi:,} ({(p1_hi - p1_lo) / 1000:.1f} kb)")
    print(f"  Prophage 2: {p2_lo:,}-{p2_hi:,} ({(p2_hi - p2_lo) / 1000:.1f} kb)")

    # Partition CDS
    p1_genes = cds_df[
        (cds_df["contig"] == "contig_1")
        & (cds_df["start"] >= p1_lo) & (cds_df["end"] <= p1_hi)
        & ~((cds_df["start"] >= p2_lo) & (cds_df["end"] <= p2_hi))
    ].sort_values("start").reset_index(drop=True)

    p2_genes = cds_df[
        (cds_df["contig"] == "contig_1")
        & (cds_df["start"] >= p2_lo) & (cds_df["end"] <= p2_hi)
    ].sort_values("start").reset_index(drop=True)

    print(f"  Prophage 1 CDS (excl. nested): {len(p1_genes)}")
    print(f"  Prophage 2 CDS: {len(p2_genes)}")

    # Extract proteins and run BLASTp
    with tempfile.TemporaryDirectory() as tmpdir:
        p1_fa = Path(tmpdir) / "p1.fasta"
        p2_fa = Path(tmpdir) / "p2.fasta"
        db_name = Path(tmpdir) / "p1_db"

        def write_proteins(genes, outpath, tag):
            with open(outpath, "w") as f:
                for _, r in genes.iterrows():
                    nt = get_cds_nt(genome_seq, r["start"], r["end"], r["strand"])
                    aa = translate_cds(nt)
                    product = r["product"].replace(" ", "_")[:40]
                    f.write(f">{r['locus_tag']}|{tag}|{product}\n{aa}\n")

        write_proteins(p1_genes, p1_fa, "P1")
        write_proteins(p2_genes, p2_fa, "P2")

        subprocess.run(
            ["makeblastdb", "-in", str(p1_fa), "-dbtype", "prot",
             "-out", str(db_name)],
            capture_output=True,
        )

        out_file = OUTPUT_DIR / "prophage2_vs_prophage1_blastp.txt"
        subprocess.run(
            ["blastp", "-query", str(p2_fa), "-db", str(db_name),
             "-outfmt", "6 qseqid sseqid pident length mismatch gapopen "
                        "qstart qend sstart send evalue bitscore qlen slen",
             "-evalue", "1e-5", "-max_target_seqs", "1",
             "-out", str(out_file)],
            capture_output=True,
        )

    # Parse results
    cols = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
            "qstart", "qend", "sstart", "send", "evalue", "bitscore",
            "qlen", "slen"]
    blastp = pd.read_csv(out_file, sep="\t", header=None, names=cols)

    blastp["P2_locus_tag"] = blastp["qseqid"].str.split("|").str[0]
    blastp["P1_locus_tag"] = blastp["sseqid"].str.split("|").str[0]
    blastp["P2_product"]   = blastp["qseqid"].str.split("|").str[2].str.replace("_", " ")
    blastp["P1_product"]   = blastp["sseqid"].str.split("|").str[2].str.replace("_", " ")
    blastp["qcov"]         = (blastp["qend"] - blastp["qstart"] + 1) / blastp["qlen"] * 100

    # Print summary table
    print(f"\n  {'P2 locus':18s}  {'P1 locus':18s}  {'%id':>6s}  {'aln':>4s}  "
          f"{'qlen':>5s}  {'slen':>5s}  {'qcov':>5s}  Product")
    print(f"  {'-' * 100}")
    for _, r in blastp.iterrows():
        print(f"  {r['P2_locus_tag']:18s}  {r['P1_locus_tag']:18s}  "
              f"{r['pident']:6.1f}  {r['length']:4d}  {r['qlen']:5d}  "
              f"{r['slen']:5d}  {r['qcov']:4.1f}%  {r['P2_product'][:45]}")

    # Summary statistics
    n_total = len(blastp)
    n_100   = len(blastp[blastp["pident"] == 100.0])
    n_ge986 = len(blastp[blastp["pident"] >= 98.6])
    p2_with_hit = blastp["P2_locus_tag"].nunique()

    print(f"\n  {'=' * 60}")
    print(f"  SUMMARY")
    print(f"  {'=' * 60}")
    print(f"  P2 genes with a P1 hit:    {p2_with_hit}/{len(p2_genes)}")
    print(f"  Paralog pairs found:       {n_total}")
    print(f"  At 100% protein identity:  {n_100}/{n_total}")
    print(f"  At ≥98.6% identity:        {n_ge986}/{n_total}")
    print(f"  Identity range:            {blastp['pident'].min():.1f}% – "
          f"{blastp['pident'].max():.1f}%")
    print(f"  Mean identity:             {blastp['pident'].mean():.1f}%")
    print(f"  Median identity:           {blastp['pident'].median():.1f}%")

    sub100 = blastp[blastp["pident"] < 100.0]
    if len(sub100) > 0:
        print(f"\n  Pairs below 100% ({len(sub100)}):")
        for _, r in sub100.iterrows():
            print(f"    {r['P2_locus_tag']} vs {r['P1_locus_tag']}: "
                  f"{r['pident']:.1f}%, {r['length']}aa — {r['P2_product'][:50]}")

    # Save CSV summary
    out_csv = OUTPUT_DIR / "prophage2_vs_prophage1_blastp_summary.csv"
    blastp[[
        "P2_locus_tag", "P1_locus_tag", "pident", "length",
        "mismatch", "qlen", "slen", "qcov", "evalue", "P2_product",
    ]].to_csv(out_csv, index=False)
    print(f"\n  Summary CSV saved: {out_csv}")

    return blastp


# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 3: Synteny figure
# ══════════════════════════════════════════════════════════════════════════════
def draw_arrow(ax, start, end, strand, y, height, color,
               edgecolor="black", lw=0.4):
    x0, x1 = start, end
    if strand == "-":
        x0, x1 = x1, x0
    arrow = FancyArrow(
        x0, y, x1 - x0, 0,
        width=height, head_width=height,
        head_length=min(abs(x1 - x0) * 0.4, 200),
        length_includes_head=True,
        facecolor=color, edgecolor=edgecolor, linewidth=lw,
    )
    ax.add_patch(arrow)


def build_figure_region(cds_df):
    """Build the annotated region DataFrame for the synteny figure."""
    region = cds_df[
        (cds_df["contig"] == "contig_1")
        & (cds_df["start"] >= WINDOW_LO)
        & (cds_df["end"] <= WINDOW_HI)
    ].sort_values("start").reset_index(drop=True)

    lt_to_gene = {v: k for k, v in CANONICAL_LTS.items()}
    region["gene_cluster"] = region["locus_tag"].map(lt_to_gene)
    region.loc[region["locus_tag"] == SDA1_PARA_LT, "gene_cluster"] = "sda1_paralog"

    def assign_module(row):
        g = row["gene_cluster"]
        if g == "sda1_paralog":
            return "cargo_sda1_paralog"
        if g in MODULE:
            return MODULE[g]
        prod = row["product"].lower()
        if any(k in prod for k in (
            "phage", "tail", "capsid", "portal", "terminase", "head",
            "morphogenesis", "siphovirus", "virion", "scaffold",
            "baseplate", "integrase", "paratox", "holin", "amidase",
            "hyaluron", "irre", "snf2", "methyltrans",
        )):
            return "other_phage"
        return "host"

    region["module"] = region.apply(assign_module, axis=1)
    return region


def make_figure(region, nt_id_str, aa_id_str, outpath_pdf, outpath_png):
    """Generate Supplementary Figure S6."""
    print("\n" + "=" * 72)
    print("ANALYSIS 3: Generating synteny figure")
    print("=" * 72)

    fig, ax = plt.subplots(figsize=(14, 4.8))

    anchor = region[region["locus_tag"] == SDA1_CANON_LT]
    x_anchor = (anchor["start"].iloc[0] + anchor["end"].iloc[0]) / 2

    x_lo = (region["start"] - x_anchor).min() - 1000
    x_hi = (region["end"] - x_anchor).max() + 1000

    arrow_h = 0.45
    y = 0.0
    ax.plot([x_lo, x_hi], [y, y], color="#888", lw=0.6, zorder=1)
    ax.text(x_lo - 1000, y, "1314.4589\n(ST28-emm1)",
            ha="right", va="center", fontsize=10, fontweight="bold")

    # Draw gene arrows
    for _, row in region.iterrows():
        s = row["start"] - x_anchor
        e = row["end"] - x_anchor
        color = COLORS[row["module"]]
        offset = arrow_h * 0.55 if row["strand"] == "+" else -arrow_h * 0.55
        draw_arrow(ax, s, e, row["strand"], y + offset, arrow_h * 0.9,
                   color, lw=0.3)

    # Annotate sda1 copies
    for lt_, label_ in [
        (SDA1_CANON_LT, "sda1\n(Prophage 1 cargo)"),
        (SDA1_PARA_LT,  "sda1 paralog\n(Prophage 2 cargo,\nsda1_00717)"),
    ]:
        r = region[region["locus_tag"] == lt_].iloc[0]
        mid = (r["start"] + r["end"]) / 2 - x_anchor
        ax.annotate(label_, xy=(mid, y + arrow_h),
                    xytext=(mid, y + arrow_h * 4.5),
                    ha="center", va="bottom", fontsize=9, fontweight="bold",
                    color="#7b241c",
                    arrowprops=dict(arrowstyle="-", color="#7b241c", lw=0.8))

    # Annotate integrases
    for lt_, label_ in [
        (P2_INTEGRASE, "Prophage 2 integrase\n(attL of Prophage 2)"),
        (CANONICAL_LTS["groups_4947"],
         "Prophage 1 integrase\n(attR of Prophage 1,\ngroups_4947)"),
    ]:
        match = region[region["locus_tag"] == lt_]
        if len(match) == 0:
            continue
        r = match.iloc[0]
        mid = (r["start"] + r["end"]) / 2 - x_anchor
        ax.annotate(label_, xy=(mid, y - arrow_h),
                    xytext=(mid, y - arrow_h * 5.5),
                    ha="center", va="top", fontsize=8.5, fontweight="bold",
                    color="#5b2c6f",
                    arrowprops=dict(arrowstyle="-", color="#5b2c6f", lw=0.8))

    # Prophage brackets
    p1_lo = region[region["locus_tag"] == P1_ATT_L]["start"].iloc[0] - x_anchor
    p1_hi = region[region["locus_tag"] == P1_ATT_R]["end"].iloc[0] - x_anchor
    p2_lo = region[region["locus_tag"] == P2_START]["start"].iloc[0] - x_anchor
    p2_hi = region[region["locus_tag"] == P2_END]["end"].iloc[0] - x_anchor

    yb1 = y + arrow_h * 13.5
    yb2 = y + arrow_h * 10.5

    ax.plot([p1_lo, p1_lo, p1_hi, p1_hi],
            [yb1 - 0.1, yb1, yb1, yb1 - 0.1], color="#7b241c", lw=1.2)
    ax.text((p1_lo + p1_hi) / 2, yb1 + 0.15,
            f"Prophage 1 (canonical Sda1-encoding prophage, "
            f"~{(p1_hi - p1_lo) / 1000:.0f} kb)",
            ha="center", va="bottom", fontsize=10, fontweight="bold",
            color="#7b241c")

    ax.plot([p2_lo, p2_lo, p2_hi, p2_hi],
            [yb2 - 0.1, yb2, yb2, yb2 - 0.1],
            color="#922b21", lw=1.0, linestyle="--")
    ax.text((p2_lo + p2_hi) / 2, yb2 + 0.15,
            f"Prophage 2 (nested Sda1-paralog prophage, "
            f"~{(p2_hi - p2_lo) / 1000:.0f} kb)",
            ha="center", va="bottom", fontsize=9, fontweight="bold",
            color="#922b21")

    # Axes
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(-arrow_h * 11, arrow_h * 18)
    xticks = list(range(int(x_lo / 1000) * 1000,
                        int(x_hi / 1000) * 1000 + 1, 5000))
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{t // 1000:+d}" if t != 0 else "0" for t in xticks])
    ax.set_xlabel("Genomic position relative to canonical sda1 start (kb)",
                  fontsize=10)
    ax.set_yticks([])
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.set_title(
        "A nested Sda1-encoding prophage in ST28-emm1 strain 1314.4589\n"
        "indicates recent intra-genomic duplication of the canonical "
        "Sda1 prophage",
        fontsize=11, pad=12,
    )

    # Identity footnote
    ax.text(
        (x_lo + x_hi) / 2, -arrow_h * 10.0,
        f"sda1 canonical vs sda1_00717 paralog: "
        f"{nt_id_str} nucleotide identity, {aa_id_str} protein identity",
        ha="center", va="top", fontsize=8.5, fontstyle="italic", color="#444",
    )

    # Legend
    legend_order = [
        "integration", "dna_replication", "structural", "paratox",
        "cargo_sda1", "cargo_sda1_paralog", "cargo_other", "lysis",
        "other_phage", "host",
    ]
    handles = [mpatches.Patch(color=COLORS[m], label=MODULE_LABELS[m])
               for m in legend_order]
    ax.legend(handles=handles, loc="lower center",
              bbox_to_anchor=(0.5, -0.4), ncol=3, fontsize=8.5,
              frameon=False)

    plt.tight_layout()
    plt.savefig(outpath_pdf, bbox_inches="tight")
    plt.savefig(outpath_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"\n  Saved: {outpath_pdf}")
    print(f"  Saved: {outpath_png}")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    print("Parsing GFF3 and genome sequence...")
    cds_df = parse_gff_cds(GFF_PATH)
    seqs = extract_fasta(GFF_PATH)
    genome_seq = seqs["contig_1"]
    print(f"  {len(cds_df)} CDS features, genome = {len(genome_seq):,} bp")

    # Analysis 1: sda1 comparison
    nt_id, aa_id = run_sda1_comparison(cds_df, genome_seq)

    # Analysis 2: Prophage-wide BLASTp
    blastp_df = run_prophage_blastp(cds_df, genome_seq)

    # Analysis 3: Synteny figure
    region = build_figure_region(cds_df)
    nt_str = f"100% (1173/1173 bp)" if nt_id == 100.0 else f"{nt_id:.1f}%"
    aa_str = f"100% (390/390 aa)" if aa_id == 100.0 else f"{aa_id:.1f}%"
    make_figure(
        region, nt_str, aa_str,
        OUTPUT_DIR / "figS6_double_sda1_prophage.pdf",
        OUTPUT_DIR / "figS6_double_sda1_prophage.png",
    )

    print("\n" + "=" * 72)
    print("ALL ANALYSES COMPLETE")
    print("=" * 72)


if __name__ == "__main__":
    main()
