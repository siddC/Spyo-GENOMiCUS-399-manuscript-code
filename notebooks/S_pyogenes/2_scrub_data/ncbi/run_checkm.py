#!/usr/bin/env python3
"""
run_checkm.py

Runs CheckM lineage_wf on a directory of .fna genome files, parses the
results into a TSV, and optionally filters by completeness/contamination.

CheckM must be installed and available in the active conda environment:
    conda activate checkm-env
    python run_checkm.py <genome_dir> <output_dir>

Usage:
    python run_checkm.py GENOME_DIR OUTPUT_DIR [options]

Example:
    python run_checkm.py \\
        /mnt/craig/pan_phylon/S_pyogenes/raw/temp_genomes \\
        /mnt/craig/pan_phylon/S_pyogenes/raw/checkm_output \\
        --threads 16 \\
        --min-completeness 97 \\
        --max-contamination 3
"""

import argparse
import csv
import shutil
import subprocess
import sys
from pathlib import Path


# ── Run CheckM ─────────────────────────────────────────────────────────────────

def run_checkm_lineage_wf(
    genome_dir: Path,
    output_dir: Path,
    threads: int,
    extension: str,
) -> int:
    """
    Run: checkm lineage_wf <genome_dir> <output_dir> -x <ext> -t <threads>
    Returns subprocess returncode.
    """
    cmd = [
        "checkm", "lineage_wf",
        str(genome_dir),
        str(output_dir),
        "-x", extension,
        "-t", str(threads),
        "--pplacer_threads", str(min(threads, 16)),  # pplacer caps at 16 usefully
    ]
    print(f"Running: {' '.join(cmd)}\n")
    result = subprocess.run(cmd)
    return result.returncode


def run_checkm_qa(
    output_dir: Path,
    qa_tsv: Path,
    threads: int,
) -> int:
    """
    Run: checkm qa <lineage.ms> <output_dir> -o 2 --tab_table -f <qa_tsv>

    -o 2  →  extended output table (includes completeness, contamination,
             strain heterogeneity, genome size, # contigs, N50, etc.)
    """
    lineage_ms = output_dir / "lineage.ms"
    cmd = [
        "checkm", "qa",
        str(lineage_ms),
        str(output_dir),
        "-o", "2",                   # extended output format
        "--tab_table",
        "-f", str(qa_tsv),
        "-t", str(threads),
    ]
    print(f"Running: {' '.join(cmd)}\n")
    result = subprocess.run(cmd)
    return result.returncode


# ── Parse + filter QA table ────────────────────────────────────────────────────

# CheckM qa -o 2 column names (tab-separated)
# Bin Id | Marker lineage | # genomes | # markers | # marker sets |
# 0 | 1 | 2 | 3 | 4 | 5+ | Completeness | Contamination | Strain heterogeneity |
# Genome size (bp) | # ambiguous bases | # scaffolds | # contigs |
# N50 (scaffolds) | N50 (contigs) | Mean scaffold length (bp) |
# Mean contig length (bp) | Longest scaffold (bp) | Longest contig (bp) |
# GC | GC std (scaffolds > 1kbp) | Coding density | Translation table |
# # predicted genes | 0 | 1 | 2 | 3 | 4 | 5+

COMPLETENESS_COL  = "Completeness"
CONTAMINATION_COL = "Contamination"


def parse_and_filter(
    qa_tsv: Path,
    filtered_tsv: Path,
    min_completeness: float,
    max_contamination: float,
) -> tuple[int, int]:
    """
    Read CheckM QA TSV, add a PASS/FAIL column, write filtered TSV.
    Returns (n_total, n_passed).
    """
    n_total  = 0
    n_passed = 0

    with open(qa_tsv, newline="", encoding="utf-8") as fin, \
         open(filtered_tsv, "w", newline="", encoding="utf-8") as fout:

        reader = csv.DictReader(fin, delimiter="\t")
        fieldnames = (reader.fieldnames or []) + ["checkm_pass"]
        writer = csv.DictWriter(fout, fieldnames=fieldnames,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        for row in reader:
            n_total += 1
            try:
                completeness  = float(row[COMPLETENESS_COL])
                contamination = float(row[CONTAMINATION_COL])
            except (KeyError, ValueError):
                row["checkm_pass"] = "MISSING"
                writer.writerow(row)
                continue

            passes = (
                completeness  >= min_completeness and
                contamination <= max_contamination
            )
            row["checkm_pass"] = "PASS" if passes else "FAIL"
            writer.writerow(row)
            if passes:
                n_passed += 1

    return n_total, n_passed


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("genome_dir",
                        help="Directory containing .fna genome files")
    parser.add_argument("output_dir",
                        help="CheckM output directory (will be created if absent)")
    parser.add_argument("--threads", type=int, default=16,
                        help="Number of threads for CheckM (default: 16)")
    parser.add_argument("--extension", default="fna",
                        help="Genome file extension (default: fna)")
    parser.add_argument("--min-completeness", type=float, default=97.0,
                        help="Minimum completeness %% to pass (default: 97.0)")
    parser.add_argument("--max-contamination", type=float, default=3.0,
                        help="Maximum contamination %% to pass (default: 3.0)")
    parser.add_argument("--skip-lineage-wf", action="store_true",
                        help="Skip lineage_wf and go straight to qa "
                             "(use if lineage_wf already completed)")
    args = parser.parse_args()

    if shutil.which("checkm") is None:
        sys.exit(
            "ERROR: checkm not found in PATH.\n"
            "Make sure you have activated the correct conda environment:\n"
            "    conda activate checkm-env"
        )

    genome_dir = Path(args.genome_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    qa_tsv       = output_dir / "checkm_qa.tsv"
    filtered_tsv = output_dir / "checkm_qa_filtered.tsv"

    # Count input genomes
    fna_files = sorted(genome_dir.glob(f"*.{args.extension}"))
    print(f"Genome files found : {len(fna_files)}")
    print(f"Output directory   : {output_dir}")
    print(f"Threads            : {args.threads}\n")

    # ── Step 1: lineage_wf ────────────────────────────────────────────────────
    if not args.skip_lineage_wf:
        rc = run_checkm_lineage_wf(genome_dir, output_dir, args.threads, args.extension)
        if rc != 0:
            sys.exit(f"ERROR: checkm lineage_wf failed with exit code {rc}")
    else:
        print("Skipping lineage_wf (--skip-lineage-wf set)\n")

    # ── Step 2: qa → TSV ─────────────────────────────────────────────────────
    rc = run_checkm_qa(output_dir, qa_tsv, args.threads)
    if rc != 0:
        sys.exit(f"ERROR: checkm qa failed with exit code {rc}")

    # ── Step 3: parse + filter ────────────────────────────────────────────────
    print(f"\nParsing {qa_tsv} ...")
    n_total, n_passed = parse_and_filter(
        qa_tsv, filtered_tsv,
        min_completeness=args.min_completeness,
        max_contamination=args.max_contamination,
    )

    n_failed = n_total - n_passed
    print(f"\nResults:")
    print(f"  Total genomes  : {n_total}")
    print(f"  PASS           : {n_passed}  (completeness >= {args.min_completeness}%, contamination <= {args.max_contamination}%)")
    print(f"  FAIL           : {n_failed}")
    print(f"\nFull results  → {qa_tsv}")
    print(f"Filtered TSV  → {filtered_tsv}")


if __name__ == "__main__":
    main()

