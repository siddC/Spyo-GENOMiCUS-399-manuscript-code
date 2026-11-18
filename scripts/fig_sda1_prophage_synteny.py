#!/usr/bin/env python3
"""
Sda1 prophage cross-phylon synteny plot (Fig. 5, Panel 2).

Draws the canonical Sda1-encoding prophage region in three representative
genomes — one per Sda1+ phylon — anchored at the sda1 start codon.
Locates and color-codes 25 of the 27 Sda1+ shared-exclusive genes; the
remaining 2 (groups_343, groups_3323) sit at distant chromosomal positions
and are noted in a footnote.

Representative genomes (highest A_norm affinity for their phylon):
    ST28-emm1:   GCF_051530335.1
    ST36:        GCF_041429935.1
    ST63-emm77:  GCF_041074825.1

Input files:
    GCF_051530335_1.gff3
    GCF_041429935_1.gff3
    GCF_041074825_1.gff3

Outputs:
    fig5_sda1_prophage_synteny.pdf
    fig5_sda1_prophage_synteny.png

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
# ---------------------------------------------------------------------------
# Input/output directories.  Defaults assume this script is run from the
# repository root with the Zenodo matrices unpacked into data/:
#     python scripts/fig_sda1_prophage_synteny.py
# Override either with an environment variable:
#     SPYO_DATA=/path/to/data SPYO_RESULTS=/path/to/out python scripts/fig_sda1_prophage_synteny.py
# ---------------------------------------------------------------------------
INPUT_DIR  = Path(os.environ.get("SPYO_DATA", "data"))
OUTPUT_DIR = Path(os.environ.get("SPYO_RESULTS", "results"))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Module assignment for the 27 Sda1+ shared-exclusive genes
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
    "integration":     "#9b59b6",
    "dna_replication": "#3498db",
    "structural":      "#2ecc71",
    "paratox":         "#e67e22",
    "cargo_sda1":      "#e74c3c",
    "cargo_other":     "#f39c12",
    "lysis":           "#1abc9c",
    "other_phage":     "#bdc3c7",
    "host":            "#7f8c8d",
}

MODULE_LABELS = {
    "integration":     "Integration / lysogeny",
    "dna_replication": "DNA modification / replication",
    "structural":      "Phage structural backbone",
    "paratox":         "Paratox cassette",
    "cargo_sda1":      "sda1 (cargo)",
    "cargo_other":     "Other cargo (HylA, HylB, GlpK, NDPK, secreted)",
    "lysis":           "Lysis module (holins, amidase)",
    "other_phage":     "Other phage genes (not in 27-gene set)",
    "host":            "Host chromosome",
}

# Short labels for annotated genes in the figure
LABEL = {
    "sda1": "sda1", "groups_141": "paratox", "groups_2030": "paratox",
    "groups_89": "amidase", "groups_1105": "holin", "groups_2032": "holin",
    "groups_4948": "holin", "groups_4947": "int", "irrE_00560": "irrE",
    "groups_1725": "meth", "groups_3351": "DNA-bind",
    "groups_2028": "helicase", "groups_1465": "hyaluronidase",
    "groups_4930": "glcnacase", "groups_3349": "glycerate-K",
    "groups_343": "NDPK", "groups_3323": "secr.",
}

# Representative genomes and their locus-tag mappings (from
# gene_pa_Sda1_exclusive_gene_set.csv)
REPS = [
    {
        "label":    "ST28-emm1\n(GCF_051530335.1)",
        "gff_path": INPUT_DIR / "GCF_051530335_1.gff3",
        "mapping": {
            "sda1": "JAOCMF_01460", "irrE_00560": "JAOCMF_01511",
            "groups_89": "JAOCMF_01461", "groups_107": "JAOCMF_01465",
            "groups_141": "JAOCMF_01459", "groups_343": "JAOCMF_00822",
            "groups_365": "JAOCMF_01505", "groups_715": "JAOCMF_01491",
            "groups_1105": "JAOCMF_01464", "groups_1281": "JAOCMF_01504",
            "groups_1462": "JAOCMF_01492", "groups_1465": "JAOCMF_01467",
            "groups_1725": "JAOCMF_01490", "groups_2024": "JAOCMF_01512",
            "groups_2028": "JAOCMF_01493", "groups_2030": "JAOCMF_01489",
            "groups_2032": "JAOCMF_01463", "groups_2526": "JAOCMF_01509",
            "groups_3323": "JAOCMF_00191", "groups_3349": "JAOCMF_01507",
            "groups_3351": "JAOCMF_01506", "groups_3362": "JAOCMF_01466",
            "groups_4930": "JAOCMF_01469", "groups_4935": "JAOCMF_01468",
            "groups_4947": "JAOCMF_01513", "groups_4948": "JAOCMF_01508",
            "groups_4950": "JAOCMF_01498",
        },
    },
    {
        "label":    "ST36\n(GCF_041429935.1)",
        "gff_path": INPUT_DIR / "GCF_041429935_1.gff3",
        "mapping": {
            "sda1": "ELMNOM_01350", "irrE_00560": "ELMNOM_01401",
            "groups_89": "ELMNOM_01351", "groups_107": "ELMNOM_01355",
            "groups_141": "ELMNOM_01349", "groups_343": "ELMNOM_00836",
            "groups_365": "ELMNOM_01395", "groups_715": "ELMNOM_01381",
            "groups_1105": "ELMNOM_01354", "groups_1281": "ELMNOM_01394",
            "groups_1462": "ELMNOM_01382", "groups_1465": "ELMNOM_01357",
            "groups_1725": "ELMNOM_01380", "groups_2024": "ELMNOM_01402",
            "groups_2028": "ELMNOM_01383", "groups_2030": "ELMNOM_01379",
            "groups_2032": "ELMNOM_01353", "groups_2526": "ELMNOM_01399",
            "groups_3323": "ELMNOM_00187", "groups_3349": "ELMNOM_01397",
            "groups_3351": "ELMNOM_01396", "groups_3362": "ELMNOM_01356",
            "groups_4930": "ELMNOM_01359", "groups_4935": "ELMNOM_01358",
            "groups_4947": "ELMNOM_01403", "groups_4948": "ELMNOM_01398",
            "groups_4950": "ELMNOM_01388",
        },
    },
    {
        "label":    "ST63-emm77\n(GCF_041074825.1)",
        "gff_path": INPUT_DIR / "GCF_041074825_1.gff3",
        "mapping": {
            "sda1": "GPDFCO_01541", "irrE_00560": "GPDFCO_01592",
            "groups_89": "GPDFCO_01542", "groups_107": "GPDFCO_01546",
            "groups_141": "GPDFCO_01540", "groups_343": "GPDFCO_00868",
            "groups_365": "GPDFCO_01586", "groups_715": "GPDFCO_01572",
            "groups_1105": "GPDFCO_01545", "groups_1281": "GPDFCO_01585",
            "groups_1462": "GPDFCO_01573", "groups_1465": "GPDFCO_01548",
            "groups_1725": "GPDFCO_01571", "groups_2024": "GPDFCO_01593",
            "groups_2028": "GPDFCO_01574", "groups_2030": "GPDFCO_01570",
            "groups_2032": "GPDFCO_01544", "groups_2526": "GPDFCO_01590",
            "groups_3323": "GPDFCO_00187", "groups_3349": "GPDFCO_01588",
            "groups_3351": "GPDFCO_01587", "groups_3362": "GPDFCO_01547",
            "groups_4930": "GPDFCO_01550", "groups_4935": "GPDFCO_01549",
            "groups_4947": "GPDFCO_01594", "groups_4948": "GPDFCO_01589",
            "groups_4950": "GPDFCO_01579",
        },
    },
]


# ══════════════════════════════════════════════════════════════════════════════
# GFF3 PARSER
# ══════════════════════════════════════════════════════════════════════════════
def parse_gff_cds(path):
    """Parse CDS features from a GFF3 file."""
    feats = []
    with open(path) as fh:
        for line in fh:
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


# ══════════════════════════════════════════════════════════════════════════════
# BUILD PROPHAGE REGION FOR EACH REPRESENTATIVE
# ══════════════════════════════════════════════════════════════════════════════
def build_region(rep):
    """
    For one representative genome: parse the GFF3, locate the in-prophage
    members of the 27-gene set (excluding 2 distant genes), determine the
    prophage extent, expand ±2 kb for flanking context, and annotate each
    CDS with its module assignment and short label.
    """
    df = parse_gff_cds(rep["gff_path"])
    distant = {"groups_343", "groups_3323"}
    in_proph_lts = [rep["mapping"][g] for g in rep["mapping"] if g not in distant]
    sub = df[df["locus_tag"].isin(in_proph_lts)]
    contig = sub["contig"].iloc[0]
    lo, hi = sub["start"].min(), sub["end"].max()
    lo_show, hi_show = lo - 2000, hi + 2000
    region = df[
        (df["contig"] == contig)
        & (df["start"] >= lo_show)
        & (df["end"] <= hi_show)
    ].sort_values("start").reset_index(drop=True)

    lt_to_gene = {v: k for k, v in rep["mapping"].items()}
    region["gene_cluster"] = region["locus_tag"].map(lt_to_gene)

    def assign_module(row):
        g = row["gene_cluster"]
        if g in MODULE:
            return MODULE[g]
        prod = row["product"].lower()
        if any(k in prod for k in (
            "phage", "tail", "capsid", "portal", "terminase", "head",
            "morphogenesis", "siphovirus", "virion", "scaffold", "baseplate",
        )):
            return "other_phage"
        return "host"

    region["module"] = region.apply(assign_module, axis=1)
    region["short_label"] = region["gene_cluster"].map(LABEL)
    rep["region"] = region
    rep["lo"] = lo
    rep["hi"] = hi
    return rep


# ══════════════════════════════════════════════════════════════════════════════
# DRAWING HELPERS
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


# ══════════════════════════════════════════════════════════════════════════════
# MAIN FIGURE
# ══════════════════════════════════════════════════════════════════════════════
def make_figure(reps, outpath_pdf, outpath_png):
    fig, ax = plt.subplots(figsize=(14, 7))
    n_tracks = len(reps)
    arrow_h = 0.45
    track_spacing = 2.4

    # Anchor each rep at sda1 midpoint
    for rep in reps:
        sda1_lt = rep["mapping"]["sda1"]
        sda1_row = rep["region"][rep["region"]["locus_tag"] == sda1_lt].iloc[0]
        rep["anchor"] = (sda1_row["start"] + sda1_row["end"]) / 2

    # Global x extent
    all_lo, all_hi = [], []
    for rep in reps:
        all_lo.append((rep["region"]["start"] - rep["anchor"]).min())
        all_hi.append((rep["region"]["end"] - rep["anchor"]).max())
    x_lo, x_hi = min(all_lo) - 500, max(all_hi) + 500

    # Draw tracks
    for i, rep in enumerate(reps):
        y = (n_tracks - 1 - i) * track_spacing
        ax.plot([x_lo, x_hi], [y, y], color="#888", lw=0.6, zorder=1)
        ax.text(x_lo - 700, y, rep["label"], ha="right", va="center",
                fontsize=10, fontweight="bold")

        # Collision-aware label placement (3 rows)
        labeled = rep["region"][
            rep["region"]["short_label"].notna()
        ].copy().sort_values("start").reset_index(drop=True)
        label_row = []
        last_x = {0: -1e12, 1: -1e12, 2: -1e12}
        MIN_SEP = 1800
        for _, r in labeled.iterrows():
            x_mid = (r["start"] + r["end"]) / 2 - rep["anchor"]
            chosen = next((k for k in (0, 1, 2)
                           if x_mid - last_x[k] >= MIN_SEP), 0)
            label_row.append(chosen)
            last_x[chosen] = x_mid
        labeled["label_row"] = label_row
        lr_lookup = dict(zip(labeled["locus_tag"], labeled["label_row"]))

        for _, row in rep["region"].iterrows():
            s = row["start"] - rep["anchor"]
            e = row["end"] - rep["anchor"]
            color = COLORS[row["module"]]
            offset = arrow_h * 0.55 if row["strand"] == "+" else -arrow_h * 0.55
            draw_arrow(ax, s, e, row["strand"], y + offset,
                       arrow_h * 0.9, color, lw=0.3)
            if pd.notna(row["short_label"]):
                ri = lr_lookup[row["locus_tag"]]
                if row["strand"] == "+":
                    ax.text((s + e) / 2, y + 0.85 + 0.30 * ri,
                            row["short_label"], ha="center", va="bottom",
                            fontsize=6.5)
                else:
                    ax.text((s + e) / 2, y - 0.85 - 0.30 * ri,
                            row["short_label"], ha="center", va="top",
                            fontsize=6.5)

    # Axes
    ax.set_xlim(x_lo, x_hi)
    ax.set_ylim(-1.7, (n_tracks - 1) * track_spacing + 1.7)
    xticks = list(range(int(x_lo / 1000) * 1000,
                        int(x_hi / 1000) * 1000 + 1, 5000))
    ax.set_xticks(xticks)
    ax.set_xticklabels([f"{t // 1000:+d}" if t != 0 else "0" for t in xticks])
    ax.set_xlabel("Genomic position relative to sda1 start (kb)", fontsize=10)
    ax.set_yticks([])
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.set_title(
        "Synteny of the canonical Sda1-encoding prophage across the "
        "three Sda1+ phylons of Spyo-GENOMiCUS-400",
        fontsize=11, pad=12,
    )

    # Legend
    legend_order = [
        "integration", "dna_replication", "structural", "paratox",
        "cargo_sda1", "cargo_other", "lysis", "other_phage", "host",
    ]
    handles = [mpatches.Patch(color=COLORS[m], label=MODULE_LABELS[m])
               for m in legend_order]
    ax.legend(handles=handles, loc="lower center",
              bbox_to_anchor=(0.5, -0.32), ncol=3, fontsize=8.5,
              frameon=False)

    # Footnote
    ax.text(
        x_lo, -1.5,
        "* groups_343 (NDP kinase) and groups_3323 (secreted protein) — "
        "the remaining 2 of 27 Sda1+ shared-exclusive genes — sit at "
        "distant chromosomal positions and are not shown.",
        ha="left", va="top", fontsize=7.5, fontstyle="italic", color="#555",
    )

    plt.tight_layout()
    plt.savefig(outpath_pdf, bbox_inches="tight")
    plt.savefig(outpath_png, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {outpath_pdf}")
    print(f"Saved: {outpath_png}")


if __name__ == "__main__":
    REPS_built = [build_region(r) for r in REPS]
    for rep in REPS_built:
        n = len(rep["region"])
        n27 = rep["region"]["gene_cluster"].isin(MODULE).sum()
        span = (rep["hi"] - rep["lo"]) / 1000
        label = rep["label"].replace("\n", " ")
        print(f"{label}: {n} CDS in window, {n27}/25 in-prophage genes, "
              f"~{span:.1f} kb span")

    make_figure(
        REPS_built,
        OUTPUT_DIR / "fig5_sda1_prophage_synteny.pdf",
        OUTPUT_DIR / "fig5_sda1_prophage_synteny.png",
    )
