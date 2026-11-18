#!/usr/bin/env python3
"""
run_bakta.py
============
Annotate all S. pyogenes genomes from the mash-scrubbed metadata CSV using Bakta.

Usage (from within the bakta-deepsig conda env):
    # Dry run — prints commands without executing
    python run_bakta.py --dry-run

    # Full run with 8 threads per genome
    python run_bakta.py --threads 8

    # Resume after a crash (skips genomes with existing output)
    python run_bakta.py --threads 8

    # Run with GNU parallel (recommended for throughput)
    python run_bakta.py --dry-run --parallel-mode > bakta_commands.txt
    parallel -j 4 < bakta_commands.txt
"""

import argparse
import subprocess
import sys
import os
import glob
import pandas as pd
from pathlib import Path

# ══════════════════════════════════════════════════════════════════════════════
# CONFIGURATION — edit these paths as needed
# ══════════════════════════════════════════════════════════════════════════════
METADATA_CSV = "/mnt/craig/pan_phylon/S_pyogenes/metadata/mash_scrubbed_genome_metadata.csv"
GENOME_DIR   = "/mnt/craig/pan_phylon/S_pyogenes/raw/mash_genomes"
OUTPUT_ROOT  = "/mnt/craig/pan_phylon/S_pyogenes/processed/bakta"
BAKTA_DB     = "/mnt/craig/DATABASES/bakta/db/db"

# Bakta flags
GENUS   = "Streptococcus"
SPECIES = "pyogenes"
GRAM    = "+"          # S. pyogenes is Gram-positive
TRANSLATION_TABLE = 11
# ══════════════════════════════════════════════════════════════════════════════


def find_genome_fna(genome_id: str, genome_dir: str) -> str | None:
    """
    Try to locate the .fna file for a given genome_id.
    Checks common naming conventions from BV-BRC downloads:
      1. {genome_id}.fna
      2. {genome_id}.fasta
      3. {genome_id}.fa
      4. {genome_id}.fna.gz  (compressed)
    """
    for ext in [".fna", ".fasta", ".fa", ".fna.gz", ".fasta.gz", ".fa.gz"]:
        candidate = os.path.join(genome_dir, f"{genome_id}{ext}")
        if os.path.isfile(candidate):
            return candidate
    return None


def build_bakta_cmd(
    fna_path: str,
    output_dir: str,
    prefix: str,
    strain: str | None,
    threads: int,
) -> list[str]:
    """Construct the bakta command list."""
    cmd = [
        "bakta",
        "--db", BAKTA_DB,
        "--genus", GENUS,
        "--species", SPECIES,
        "--gram", GRAM,
        "--translation-table", str(TRANSLATION_TABLE),
        "--complete",          # these are all complete genomes
        "--output", output_dir,
        "--prefix", prefix,
        "--threads", str(threads),
        "--force",
        "--verbose",
    ]
    if strain and str(strain) != "nan":
        cmd.extend(["--strain", str(strain)])

    cmd.append(fna_path)
    return cmd


def main():
    parser = argparse.ArgumentParser(description="Run Bakta on all S. pyogenes genomes")
    parser.add_argument("--threads", type=int, default=4,
                        help="Threads per bakta job (default: 4)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print commands without executing")
    parser.add_argument("--parallel-mode", action="store_true",
                        help="Output one command per line for GNU parallel")
    parser.add_argument("--metadata", default=METADATA_CSV,
                        help="Path to metadata CSV")
    parser.add_argument("--genome-dir", default=GENOME_DIR,
                        help="Directory containing .fna files")
    parser.add_argument("--output-root", default=OUTPUT_ROOT,
                        help="Root output directory for bakta results")
    args = parser.parse_args()

    # ── Load metadata ──
    df = pd.read_csv(args.metadata, low_memory=False)
    print(f"[INFO] Loaded {len(df)} genomes from metadata", file=sys.stderr)

    os.makedirs(args.output_root, exist_ok=True)

    missing = []
    skipped = []
    queued = []

    for _, row in df.iterrows():
        genome_id = str(row["genome_id"])
        strain = row.get("strain", None)
        prefix = genome_id  # use genome_id as output prefix

        # Locate the fna file
        fna_path = find_genome_fna(genome_id, args.genome_dir)
        if fna_path is None:
            missing.append(genome_id)
            continue

        # Output directory per genome
        out_dir = os.path.join(args.output_root, genome_id)

        # Skip if already annotated (check for .gff3 output)
        gff_glob = glob.glob(os.path.join(out_dir, "*.gff3"))
        if gff_glob and not args.dry_run:
            skipped.append(genome_id)
            continue

        cmd = build_bakta_cmd(fna_path, out_dir, prefix, strain, args.threads)
        queued.append((genome_id, cmd))

    # ── Report ──
    print(f"[INFO] Found:   {len(df) - len(missing)}/{len(df)}", file=sys.stderr)
    print(f"[INFO] Missing: {len(missing)}", file=sys.stderr)
    print(f"[INFO] Skipped: {len(skipped)} (already annotated)", file=sys.stderr)
    print(f"[INFO] To run:  {len(queued)}", file=sys.stderr)

    if missing:
        print(f"\n[WARN] Could not find .fna files for:", file=sys.stderr)
        for gid in missing[:20]:
            print(f"       {gid}", file=sys.stderr)
        if len(missing) > 20:
            print(f"       ... and {len(missing)-20} more", file=sys.stderr)

    # ── Execute or print ──
    if args.parallel_mode:
        # Output shell commands for GNU parallel
        for genome_id, cmd in queued:
            # Make sure output dir exists
            out_dir = os.path.join(args.output_root, genome_id)
            shell_cmd = f"mkdir -p {out_dir} && " + " ".join(f"'{c}'" for c in cmd)
            print(shell_cmd)
        return

    if args.dry_run:
        print(f"\n[DRY RUN] Would execute {len(queued)} bakta jobs:\n",
              file=sys.stderr)
        for genome_id, cmd in queued[:5]:
            print(" ".join(cmd))
        if len(queued) > 5:
            print(f"... and {len(queued)-5} more")
        return

    # Sequential execution with progress tracking
    total = len(queued)
    failed = []
    for i, (genome_id, cmd) in enumerate(queued, 1):
        out_dir = os.path.join(args.output_root, genome_id)
        os.makedirs(out_dir, exist_ok=True)

        print(f"\n[{i}/{total}] Annotating {genome_id} ...", file=sys.stderr)
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=1800,  # 30 min timeout per genome
            )
            if result.returncode != 0:
                print(f"[ERROR] bakta failed for {genome_id}", file=sys.stderr)
                print(result.stderr[-500:], file=sys.stderr)
                failed.append(genome_id)
            else:
                print(f"[OK] {genome_id} done", file=sys.stderr)
        except subprocess.TimeoutExpired:
            print(f"[TIMEOUT] {genome_id} exceeded 30 min", file=sys.stderr)
            failed.append(genome_id)
        except Exception as e:
            print(f"[ERROR] {genome_id}: {e}", file=sys.stderr)
            failed.append(genome_id)

    # ── Final report ──
    print(f"\n{'='*60}", file=sys.stderr)
    print(f"BAKTA ANNOTATION COMPLETE", file=sys.stderr)
    print(f"  Succeeded: {total - len(failed)}/{total}", file=sys.stderr)
    print(f"  Failed:    {len(failed)}", file=sys.stderr)
    if failed:
        print(f"  Failed IDs: {failed}", file=sys.stderr)
    print(f"{'='*60}", file=sys.stderr)


if __name__ == "__main__":
    main()

