#!/usr/bin/env python3
"""
Figure 6: The speC/mf2 virulence cassette on three structurally distinct prophages.

Draws synteny of the speC/mf2 region in three representative genomes
(one per paralog pair), anchored at the speC start codon. The conserved
cargo cassette (speC + mf2, head-to-tail) is highlighted, with the
prophage backbone upstream and chromosomal flank downstream. All three
prophage backbones encode hyaluronidase (HylA) upstream of the cassette.

The GCF_041430965.1 prophage is integrated in the opposite orientation
relative to the chromosome; its coordinates are reverse-complemented
for direct visual comparison.

Representative genomes:
    1314.985            ST52-emm28,          speC + mf2
    1314.3072           ST39-emm4,           speC_01049 + mf2_01051
    GCF_041430965.1     emm-patternE-other,  speC_01050 + mf2_01052

Input files:
    1314_985.gff3
    1314_3072.gff3
    GCF_041430965_1.gff3

Outputs:
    fig6_speC_mf2_synteny.pdf
    fig6_speC_mf2_synteny.png

Author: Siddharth M. Chauhan
"""

import re
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrow, Rectangle
from pathlib import Path

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ══════════════════════════════════════════════════════════════════════════════
OUTPUT_DIR = Path(".")

COLORS = {
    "Phage structural":      "#7f8c8d",   # grey-brown
    "Hyaluronidase / HylA":  "#e67e22",   # orange
    "Phage lysis":           "#8d6e3f",   # brown
    "speC":                  "#c0392b",   # dark red
    "mf2":                   "#6b1a1a",   # very dark red
    "Paratox":               "#c9a84c",   # gold
    "Chromosomal":           "#95a5a6",   # grey
}

REPS = [
    {
        "label":         "1314.985\nST52-emm28\nspeC + mf2",
        "genome_id":     "1314.985",
        "gff":           "1314_985.gff3",
        "speC_lt":       "GGAFNI_00883",
        "mf2_lt":        "GGAFNI_00884",
        "reverse_comp":  False,
    },
    {
        "label":         "1314.3072\nST39-emm4\nspeC_01049\n+ mf2_01051",
        "genome_id":     "1314.3072",
        "gff":           "1314_3072.gff3",
        "speC_lt":       "FEFGFE_00607",
        "mf2_lt":        "FEFGFE_00608",
        "reverse_comp":  False,
    },
    {
        "label":         "GCF_041430965.1\nemm-patternE-other\nspeC_01050\n+ mf2_01052",
        "genome_id":     "GCF_041430965.1",
        "gff":           "GCF_041430965_1.gff3",
        "speC_lt":       "LABOOI_01147",
        "mf2_lt":        "LABOOI_01146",
        "reverse_comp":  True,
    },
]

FLANK_FEATURES = 10


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


def assign_category(row, speC_lt, mf2_lt):
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
# BUILD DISPLAY-READY REGION
# ══════════════════════════════════════════════════════════════════════════════
def build_region(rep):
    """
    Parse GFF3, extract window around speC, assign categories,
    compute display coordinates (reverse-complementing if needed).
    """
    df = parse_gff_cds(rep["gff"])
    df = df.sort_values("start").reset_index(drop=True)

    speC_idx = df[df["locus_tag"] == rep["speC_lt"]].index[0]
    speC_row = df.loc[speC_idx]

    lo = max(0, speC_idx - FLANK_FEATURES)
    hi = min(len(df) - 1, speC_idx + FLANK_FEATURES)
    region = df.iloc[lo : hi + 1].copy()

    region["category"] = region.apply(
        lambda r: assign_category(r, rep["speC_lt"], rep["mf2_lt"]), axis=1
    )

    if rep["reverse_comp"]:
        speC_end = speC_row["end"]
        region["disp_start"] = (speC_end - region["end"]) / 1000
        region["disp_end"]   = (speC_end - region["start"]) / 1000
        region["disp_strand"] = region["strand"].map({"+": "-", "-": "+"})
    else:
        anchor = speC_row["start"]
        region["disp_start"] = (region["start"] - anchor) / 1000
        region["disp_end"]   = (region["end"] - anchor) / 1000
        region["disp_strand"] = region["strand"]

    region = region.sort_values("disp_start").reset_index(drop=True)
    rep["region"] = region
    return rep


# ══════════════════════════════════════════════════════════════════════════════
# DRAWING
# ══════════════════════════════════════════════════════════════════════════════
def draw_arrow(ax, x0, x1, strand, y, height, color, edgecolor="black", lw=0.5):
    if strand == "-":
        x0, x1 = x1, x0
    arrow = FancyArrow(
        x0, y, x1 - x0, 0,
        width=height, head_width=height,
        head_length=min(abs(x1 - x0) * 0.3, 0.15),
        length_includes_head=True,
        facecolor=color, edgecolor=edgecolor, linewidth=lw,
    )
    ax.add_patch(arrow)


def make_figure(reps, outpath_pdf, outpath_png):
    fig, ax = plt.subplots(figsize=(12, 6))

    n_tracks = len(reps)
    track_spacing = 2.2
    arrow_h = 0.38

    # Global x extent
    all_lo = min(r["region"]["disp_start"].min() for r in reps)
    all_hi = max(r["region"]["disp_end"].max() for r in reps)
    x_lo, x_hi = all_lo - 0.5, all_hi + 0.5

    # Find cargo cassette extent for highlight box
    cargo_lo = min(
        r["region"][r["region"]["category"].isin(["speC", "mf2"])]["disp_start"].min()
        for r in reps
    )
    cargo_hi = max(
        r["region"][r["region"]["category"].isin(["speC", "mf2"])]["disp_end"].max()
        for r in reps
    )

    # Draw cargo cassette highlight (vertical band)
    y_bottom = -1.0
    y_top = (n_tracks - 1) * track_spacing + 1.2
    ax.add_patch(Rectangle(
        (cargo_lo - 0.05, y_bottom), cargo_hi - cargo_lo + 0.1, y_top - y_bottom,
        facecolor="#f5d5d0", edgecolor="none", alpha=0.5, zorder=0,
    ))

    # attR dashed line
    attr_x = cargo_hi + 0.15
    ax.axvline(x=attr_x, color="#888", linestyle=":", linewidth=0.8, zorder=0)
    ax.text(attr_x + 0.05, y_top - 0.1, "attR", fontsize=8, fontstyle="italic",
            color="#666", va="top")

    # Section labels at top
    phage_x = (x_lo + cargo_lo) / 2
    cargo_x = (cargo_lo + cargo_hi) / 2
    chromo_x = (cargo_hi + x_hi) / 2
    label_y = y_top + 0.2
    for x, txt, style in [
        (phage_x, "Prophage backbone", "italic"),
        (cargo_x, "Conserved cargo cassette", "italic"),
        (chromo_x, "Chromosomal flank", "italic"),
    ]:
        ax.text(x, label_y, txt, ha="center", va="bottom", fontsize=9,
                fontstyle=style, color="#555")

    # Label speC and mf2 above the cargo highlight
    ax.text(cargo_x, y_top - 0.35, "speC    mf2", ha="center", va="bottom",
            fontsize=9, fontweight="bold", color="#7b241c")

    # Draw tracks
    for i, rep in enumerate(reps):
        y = (n_tracks - 1 - i) * track_spacing

        # Backbone line
        ax.plot([x_lo, x_hi], [y, y], color="#aaa", lw=0.5, zorder=1)

        # Track label (left)
        ax.text(x_lo - 0.3, y, rep["label"], ha="right", va="center",
                fontsize=8, fontweight="bold", linespacing=1.3)

        # Draw genes
        for _, row in rep["region"].iterrows():
            color = COLORS.get(row["category"], "#bbb")
            offset = arrow_h * 0.55 if row["disp_strand"] == "+" else -arrow_h * 0.55
            draw_arrow(ax, row["disp_start"], row["disp_end"],
                       row["disp_strand"], y + offset, arrow_h * 0.85,
                       color, lw=0.4)

    # Axes
    ax.set_xlim(x_lo - 0.3, x_hi + 0.3)
    ax.set_ylim(y_bottom - 0.5, label_y + 0.5)
    xticks = list(range(int(np.floor(x_lo)), int(np.ceil(x_hi)) + 1, 3))
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{t:+d}" if t != 0 else "0" for t in xticks])
    ax.set_xlabel("Position relative to speC start (kb)", fontsize=10)
    ax.set_yticks([])
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)

    ax.set_title(
        "Synteny diagram for speC/mf2 gene cassette",
        fontsize=12, fontweight="bold", pad=15,
    )

    # Legend
    legend_order = [
        "Phage structural", "Hyaluronidase / HylA", "Phage lysis",
        "speC", "mf2", "Paratox", "Chromosomal",
    ]
    handles = [mpatches.Patch(color=COLORS[c], label=c) for c in legend_order]
    ax.legend(
        handles=handles, loc="lower center",
        bbox_to_anchor=(0.5, -0.22), ncol=4, fontsize=8,
        frameon=False,
    )

    # Footnote about reverse complement
    ax.text(
        x_lo, y_bottom - 0.3,
        "*The GCF_041430965.1 prophage is integrated in the opposite "
        "orientation relative to the chromosome;\n"
        " its coordinates have been reverse-complemented for visual comparison.",
        ha="left", va="top", fontsize=7, fontstyle="italic", color="#777",
    )

    plt.tight_layout()
    plt.savefig(outpath_pdf, bbox_inches="tight")
    plt.savefig(outpath_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outpath_pdf}")
    print(f"Saved: {outpath_png}")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    built = [build_region(r) for r in REPS]

    for r in built:
        n = len(r["region"])
        span = r["region"]["disp_end"].max() - r["region"]["disp_start"].min()
        label = r["genome_id"]
        rc = " (reverse-complemented)" if r["reverse_comp"] else ""
        print(f"{label}: {n} CDS, {span:.1f} kb display span{rc}")

    make_figure(
        built,
        OUTPUT_DIR / "fig6_speC_mf2_synteny.pdf",
        OUTPUT_DIR / "fig6_speC_mf2_synteny.png",
    )
