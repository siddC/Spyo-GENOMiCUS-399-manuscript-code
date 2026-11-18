#!/usr/bin/env python3
"""
Verification of the Results section: "NMF recovers emm-pattern regulons
and emm-switching events"

Recomputes:
  1. Composite phylon composition (STs, emm types)
  2. Exclusive genes for emm-patternD-other and emm-patternE-other
  3. emmtyper validation (emm-cluster/pattern concordance)
  4. ST36 emm-switching event (emm12/emm82 composition)
  5. Concordance metrics (ARI/NMI for 7- and 9-phylon comparisons)
  6. Per-phylon purity

Input files:
    - L_bin_panta_cdhit_Spyo.csv
    - A_bin_panta_cdhit_Spyo.csv
    - car_metadata_Sypo.csv
    - gene_annot_panta_Spyo.csv

Author: Siddharth M. Chauhan
"""

import pandas as pd
import numpy as np
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import os
from pathlib import Path

# ---------------------------------------------------------------------------
# Input/output directories.  Defaults assume this script is run from the
# repository root with the Zenodo matrices unpacked into data/:
#     python scripts/regulon_emmswitching_verification.py
# Override either with an environment variable:
#     SPYO_DATA=/path/to/data SPYO_RESULTS=/path/to/out python scripts/regulon_emmswitching_verification.py
# ---------------------------------------------------------------------------
INPUT_DIR  = Path(os.environ.get("SPYO_DATA", "data"))
OUTPUT_DIR = Path(os.environ.get("SPYO_RESULTS", "results"))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

L = pd.read_csv(INPUT_DIR / "L_bin_panta_cdhit_Spyo.csv", index_col=0)
A = pd.read_csv(INPUT_DIR / "A_bin_panta_cdhit_Spyo.csv", index_col=0)
meta = pd.read_csv(INPUT_DIR / "car_metadata_Sypo.csv", low_memory=False)
gene_annot = pd.read_csv(INPUT_DIR / "gene_annot_panta_Spyo.csv", low_memory=False)
annot_map = dict(zip(gene_annot["Gene"], gene_annot["Annotation"]))

TAX_PHYLONS = [p for p in A.index if p not in ["SpeK-prophage", "SpeC-prophage"]]
LINEAGE_PHYLONS = [
    "ST28-emm1", "ST36", "ST63-emm77", "ST52-emm28",
    "ST39-emm4", "emm3/6/5", "ST101-emm89",
]


def classify_pattern(cluster):
    """Classify an emm cluster into A-C, D, E, or Clade Y."""
    if pd.isna(cluster):
        return "Unknown"
    c = str(cluster)
    if c.startswith("D"):
        return "D-family"
    if c.startswith("E"):
        return "E-family"
    if c.startswith("A-C"):
        return "A-C family"
    if "clade Y" in c or "Single" in c:
        return "Clade Y"
    if c == "-":
        return "Unknown"
    return f"Other ({c})"


def base_emm(val):
    """Convert raw emm subtype (emm_12.0) to base type (emm12)."""
    s = str(val)
    if s in ("-", "nan"):
        return None
    return s.replace("emm_", "emm").split(".")[0]


# ══════════════════════════════════════════════════════════════════════════════
# 1. COMPOSITE PHYLON COMPOSITION
# ══════════════════════════════════════════════════════════════════════════════
print("=" * 72)
print("1. COMPOSITE PHYLON COMPOSITION")
print("=" * 72)

for p in ["emm-clusterD-other", "emm-clusterE-other"]:
    strains = list(A.columns[A.loc[p] == 1])
    sub = meta[meta["genome_id"].astype(str).isin(strains)]
    print(f"\n{p} (n={len(strains)}):")
    print(f"  Unique STs: {sub['mlst'].nunique()}")
    print(f"  Unique emm types: {sub['Predicted emm-type'].nunique()}")

# ══════════════════════════════════════════════════════════════════════════════
# 2. EXCLUSIVE GENES
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'=' * 72}")
print("2. EXCLUSIVE GENES")
print("=" * 72)

for p in ["emm-clusterD-other", "emm-clusterE-other"]:
    others = [q for q in L.columns if q != p]
    exclusive = L.index[(L[p] == 1) & (L[others] == 0).all(axis=1)]
    print(f"\n{p}: {len(exclusive)} exclusive genes")
    for g in exclusive:
        a = str(annot_map.get(g, "")).split(",")[0][:65]
        print(f"  {g:20s} {a}")

# ══════════════════════════════════════════════════════════════════════════════
# 3. EMMTYPER VALIDATION
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'=' * 72}")
print("3. EMMTYPER VALIDATION")
print("=" * 72)

for p in ["emm-clusterD-other", "emm-clusterE-other"]:
    strains = list(A.columns[A.loc[p] == 1])
    sub = meta[meta["genome_id"].astype(str).isin(strains)].copy()
    n_total = len(sub)

    # Full distribution (all strains)
    sub["pattern"] = sub["emm cluster"].apply(classify_pattern)
    pattern_counts = sub["pattern"].value_counts()
    print(f"\n{p} (n={n_total} total):")
    for pat, count in pattern_counts.items():
        print(f"  {pat}: {count}")

    # Exclude uninterpretable emmtyper assignments from denominator
    # (emm cluster = "-" or compound/ambiguous entries starting with "(")
    evaluable = sub[
        ~sub["emm cluster"].astype(str).isin(["-", "nan"])
        & ~sub["emm cluster"].astype(str).str.startswith("(")
    ]
    n_eval = len(evaluable)
    eval_patterns = evaluable["pattern"].value_counts()

    if "D" in p:
        d_y = sum(c for pat, c in eval_patterns.items() if pat in ["D-family", "Clade Y"])
        print(f"  → D-family + Clade Y: {d_y}/{n_eval} ({d_y / n_eval * 100:.0f}%) "
              f"[{n_total - n_eval} uninterpretable excluded]")
    else:
        e = eval_patterns.get("E-family", 0)
        print(f"  → E-family: {e}/{n_eval} ({e / n_eval * 100:.0f}%) "
              f"[{n_total - n_eval} uninterpretable excluded]")

# ══════════════════════════════════════════════════════════════════════════════
# 4. ST36 EMM-SWITCHING EVENT
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'=' * 72}")
print("4. ST36 EMM-SWITCHING EVENT")
print("=" * 72)

st36 = list(A.columns[A.loc["ST36"] == 1])
sub = meta[meta["genome_id"].astype(str).isin(st36)]
print(f"\nST36 phylon: {len(st36)} strains")
print(f"\nemm types: {sub['Predicted emm-type'].value_counts().to_dict()}")
print(f"STs: {sub['mlst'].value_counts().to_dict()}")

emm12 = sub[sub["Predicted emm-type"].astype(str).str.contains("emm_12", na=False)]
emm82 = sub[sub["Predicted emm-type"].astype(str).str.contains("emm_82", na=False)]
print(f"\nemm12: {len(emm12)} strains (STs: {emm12['mlst'].value_counts().to_dict()})")
print(f"emm82: {len(emm82)} strains (STs: {emm82['mlst'].value_counts().to_dict()})")

# ══════════════════════════════════════════════════════════════════════════════
# 5. CONCORDANCE METRICS
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'=' * 72}")
print("5. CONCORDANCE METRICS (ARI / NMI)")
print("=" * 72)

# Build phylon labels
phylon_labels = {}
for p in TAX_PHYLONS:
    for s in A.columns[A.loc[p] == 1]:
        if s not in phylon_labels:
            phylon_labels[s] = p
assigned_meta = meta[meta["genome_id"].astype(str).isin(phylon_labels)].copy()
assigned_meta["phylon"] = assigned_meta["genome_id"].astype(str).map(phylon_labels)

# Derive base emm type (emm_12.0 → emm12, emm_82.7 → emm82)
assigned_meta["emm_base"] = assigned_meta["Predicted emm-type"].apply(base_emm)

# 7-phylon (lineage-specific only)
lin_meta = assigned_meta[assigned_meta["phylon"].isin(LINEAGE_PHYLONS)]

# vs MLST (all strains)
truth = lin_meta["mlst"].astype(str).values
pred = lin_meta["phylon"].values
ari = adjusted_rand_score(truth, pred)
nmi = normalized_mutual_info_score(truth, pred)
print(f"\n7-phylon vs MLST: ARI={ari:.3f}, NMI={nmi:.3f} (n={len(lin_meta)})")

# vs emm type (exclude untyped)
lin_emm = lin_meta[lin_meta["emm_base"].notna()]
truth = lin_emm["emm_base"].values
pred = lin_emm["phylon"].values
ari = adjusted_rand_score(truth, pred)
nmi = normalized_mutual_info_score(truth, pred)
print(f"7-phylon vs emm type (base): ARI={ari:.3f}, NMI={nmi:.3f} (n={len(lin_emm)})")

# 9-phylon (all taxonomic)
truth = assigned_meta["mlst"].astype(str).values
pred = assigned_meta["phylon"].values
ari = adjusted_rand_score(truth, pred)
nmi = normalized_mutual_info_score(truth, pred)
print(f"\n9-phylon vs MLST: ARI={ari:.3f}, NMI={nmi:.3f} (n={len(assigned_meta)})")

all_emm = assigned_meta[assigned_meta["emm_base"].notna()]
truth = all_emm["emm_base"].values
pred = all_emm["phylon"].values
ari = adjusted_rand_score(truth, pred)
nmi = normalized_mutual_info_score(truth, pred)
print(f"9-phylon vs emm type (base): ARI={ari:.3f}, NMI={nmi:.3f} (n={len(all_emm)})")

# ══════════════════════════════════════════════════════════════════════════════
# 6. PER-PHYLON PURITY
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'=' * 72}")
print("6. PER-PHYLON PURITY")
print("=" * 72)

for p in LINEAGE_PHYLONS:
    strains = list(A.columns[A.loc[p] == 1])
    sub = meta[meta["genome_id"].astype(str).isin(strains)].copy()
    n = len(sub)
    st_pur = sub["mlst"].value_counts().iloc[0] / n * 100
    dom_st = sub["mlst"].value_counts().index[0]

    sub["emm_base"] = sub["Predicted emm-type"].apply(base_emm)

    if p == "emm3/6/5":
        combined = sub["emm_base"].isin(["emm3", "emm5", "emm6"])
        emm_pur = combined.sum() / n * 100
        emm_str = f"emm3+5+6={combined.sum()}/{n}"
    else:
        emm_dom = sub["emm_base"].value_counts()
        emm_pur = emm_dom.iloc[0] / n * 100 if len(emm_dom) > 0 else 0
        emm_str = f"{emm_dom.index[0]}={emm_dom.iloc[0]}/{n}"

    print(f"  {p:20s}: ST{dom_st} {st_pur:.0f}%, {emm_str} {emm_pur:.0f}%")

print(f"\n{'=' * 72}")
print("VERIFICATION COMPLETE")
print("=" * 72)
