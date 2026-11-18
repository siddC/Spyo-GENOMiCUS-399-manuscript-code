#!/usr/bin/env python3
"""
Supplementary Figure S6: Double Sda1-encoding prophage in 1314.4589.

Draws the nested prophage region in ST28-emm1 strain 1314.4589, which
carries two tandemly integrated Sda1-encoding prophages on the same
contig. The canonical prophage (Prophage 1, ~68 kb) contains the 25
in-locus members of the 27-gene Sda1+ shared-exclusive set. A second,
nested prophage (Prophage 2, ~28 kb) carries a paralog of sda1
(PanTA cluster sda1_00717), which is 100% identical at both
nucleotide (1173/1173 bp) and protein (390/390 aa) level to the
canonical copy.

Input:
    1314_4589.gff3  (Bakta-annotated GFF3 for strain 1314.4589)

Outputs:
    figS6_double_sda1_prophage.pdf
    figS6_double_sda1_prophage.png

Author: Siddharth M. Chauhan
"""

import re
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrow
import os
from pathlib import Path

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
GFF_PATH   = INPUT_DIR / "1314_4589.gff3"
# ---------------------------------------------------------------------------
# Input/output directories.  Defaults assume this script is run from the
# repository root with the Zenodo matrices unpacked into data/:
#     python scripts/figS_double_sda1_prophage.py
# Override either with an environment variable:
#     SPYO_DATA=/path/to/data SPYO_RESULTS=/path/to/out python scripts/figS_double_sda1_prophage.py
# ---------------------------------------------------------------------------
INPUT_DIR  = Path(os.environ.get("SPYO_DATA", "data"))
OUTPUT_DIR = Path(os.environ.get("SPYO_RESULTS", "results"))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Genomic window to display (bp coordinates on contig_1)
WINDOW_LO = 1394000
WINDOW_HI = 1480000

# 27-gene set module assignments
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
    "integration":          "#9b59b6",
    "dna_replication":      "#3498db",
    "structural":           "#2ecc71",
    "paratox":              "#e67e22",
    "cargo_sda1":           "#e74c3c",
    "cargo_sda1_paralog":   "#c0392b",
    "cargo_other":          "#f39c12",
    "lysis":                "#1abc9c",
    "other_phage":          "#bdc3c7",
    "host":                 "#7f8c8d",
}

MODULE_LABELS = {
    "integration":          "Integration / lysogeny",
    "dna_replication":      "DNA modification / replication",
    "structural":           "Phage structural backbone",
    "paratox":              "Paratox cassette",
    "cargo_sda1":           "sda1 (cargo)",
    "cargo_sda1_paralog":   "sda1 paralog (sda1_00717)",
    "cargo_other":          "Other cargo",
    "lysis":                "Lysis module (holins, amidase)",
    "other_phage":          "Other phage genes (not in 27-gene set)",
    "host":                 "Host chromosome",
}

# Locus-tag mapping for the 27-gene set in 1314.4589
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

# sda1_00717 paralog locus tag
SDA1_PARALOG_LT = "BENLNG_01492"

# Prophage boundary locus tags
P1_ATT_L = "BENLNG_01456"   # Paratox at 5' end of Prophage 1
P1_ATT_R = "BENLNG_01545"   # Integrase at 3' end of Prophage 1
P2_START = "BENLNG_01470"   # Start of nested Prophage 2
P2_END   = "BENLNG_01504"   # End of nested Prophage 2


# ══════════════════════════════════════════════════════════════════════════════
# GFF3 PARSER
# ══════════════════════════════════════════════════════════════════════════════
def parse_gff_cds(path):
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


# ══════════════════════════════════════════════════════════════════════════════
# BUILD REGION
# ══════════════════════════════════════════════════════════════════════════════
def build_region():
    df = parse_gff_cds(GFF_PATH)
    region = df[
        (df["contig"] == "contig_1")
        & (df["start"] >= WINDOW_LO)
        & (df["end"] <= WINDOW_HI)
    ].sort_values("start").reset_index(drop=True)

    lt_to_gene = {v: k for k, v in CANONICAL_LTS.items()}
    region["gene_cluster"] = region["locus_tag"].map(lt_to_gene)
    region.loc[region["locus_tag"] == SDA1_PARALOG_LT, "gene_cluster"] = "sda1_paralog"

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


# ══════════════════════════════════════════════════════════════════════════════
# MAIN FIGURE
# ══════════════════════════════════════════════════════════════════════════════
def make_figure(region, outpath_pdf, outpath_png):
    fig, ax = plt.subplots(figsize=(14, 4.8))

    anchor = region[region["locus_tag"] == CANONICAL_LTS["sda1"]]
    x_anchor = (anchor["start"].iloc[0] + anchor["end"].iloc[0]) / 2

    rel_s = region["start"] - x_anchor
    rel_e = region["end"] - x_anchor
    x_lo, x_hi = rel_s.min() - 1000, rel_e.max() + 1000

    arrow_h = 0.45
    y = 0.0
    ax.plot([x_lo, x_hi], [y, y], color="#888", lw=0.6, zorder=1)
    ax.text(x_lo - 1000, y, "1314.4589\n(ST28-emm1)",
            ha="right", va="center", fontsize=10, fontweight="bold")

    for _, row in region.iterrows():
        s = row["start"] - x_anchor
        e = row["end"] - x_anchor
        color = COLORS[row["module"]]
        offset = arrow_h * 0.55 if row["strand"] == "+" else -arrow_h * 0.55
        draw_arrow(ax, s, e, row["strand"], y + offset, arrow_h * 0.9,
                   color, lw=0.3)

    # Annotate sda1 copies
    for lt_, label_ in [
        (CANONICAL_LTS["sda1"], "sda1\n(Prophage 1 cargo)"),
        (SDA1_PARALOG_LT, "sda1 paralog\n(Prophage 2 cargo,\nsda1_00717)"),
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
        ("BENLNG_01490", "Prophage 2 integrase\n(attL of Prophage 2)"),
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

    # BLAST identity footnote
    ax.text(
        (x_lo + x_hi) / 2, -arrow_h * 10.0,
        "sda1 canonical vs sda1_00717 paralog: "
        "100% nucleotide identity (1173/1173 bp), "
        "100% protein identity (390/390 aa)",
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
    print(f"Saved: {outpath_pdf}")
    print(f"Saved: {outpath_png}")


if __name__ == "__main__":
    region = build_region()
    n_canonical = region["gene_cluster"].isin(MODULE).sum()
    n_paralog = (region["gene_cluster"] == "sda1_paralog").sum()
    print(f"1314.4589: {len(region)} CDS in window, "
          f"{n_canonical} canonical 27-gene set, {n_paralog} sda1 paralog")

    make_figure(
        region,
        OUTPUT_DIR / "figS6_double_sda1_prophage.pdf",
        OUTPUT_DIR / "figS6_double_sda1_prophage.png",
    )
