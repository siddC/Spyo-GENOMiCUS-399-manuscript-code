#!/usr/bin/env bash
# copy_to_mash.sh
#
# Copies all genome .fna files listed in unified_genome_metadata.csv into
# mash_genomes/, looking in:
#   BV-BRC: /mnt/craig/pan_phylon/S_pyogenes/raw/genomes/{genome_id}/{genome_id}.fna
#   NCBI:   /mnt/craig/pan_phylon/S_pyogenes/raw/temp_genomes/{accession}.fna

set -euo pipefail

METADATA="/mnt/craig/pan_phylon/S_pyogenes/metadata/unified_genome_metadata.csv"
BVBRC_DIR="/mnt/craig/pan_phylon/S_pyogenes/raw/genomes"
NCBI_DIR="/mnt/craig/pan_phylon/S_pyogenes/raw/temp_genomes"
OUT_DIR="/mnt/craig/pan_phylon/S_pyogenes/raw/mash_genomes"

mkdir -p "$OUT_DIR"

python3 - "$METADATA" "$BVBRC_DIR" "$NCBI_DIR" "$OUT_DIR" << 'PYEOF'
import csv, sys, shutil
from pathlib import Path

metadata, bvbrc_dir, ncbi_dir, out_dir = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

n_ok = n_fail = 0
with open(metadata, newline="") as f:
    for row in csv.DictReader(f):
        source  = row.get("source", "").strip()
        acc     = row.get("assembly_accession_resolved", "").strip()
        gid     = row.get("genome_id", "").strip()
        file_id = gid if source == "BV-BRC" and gid else acc

        if source == "BV-BRC":
            src = Path(bvbrc_dir) / file_id / f"{file_id}.fna"
        else:
            src = Path(ncbi_dir) / f"{file_id}.fna"

        dest = Path(out_dir) / f"{file_id}.fna"
        if src.exists():
            shutil.copy2(src, dest)
            n_ok += 1
        else:
            print(f"MISSING: {src}")
            n_fail += 1

print(f"\nDone.\n  Copied  : {n_ok}\n  Missing : {n_fail}\n  Output  : {out_dir}")
PYEOF