#!/usr/bin/env bash
# run_mash.sh
#
# Runs Mash sketch + dist on all .fna files in the mash_genomes directory.
#
# Steps:
#   1. mash sketch  — sketch all genomes into a single .msh file
#   2. mash dist    — compute all-vs-all distances
#   3. Output a TSV with columns: query, reference, distance, p-value, hashes
#
# Usage:
#   ./run_mash.sh [--threads N]
#
# Output:
#   /mnt/craig/pan_phylon/S_pyogenes/raw/mash/S_pyogenes.msh
#   /mnt/craig/pan_phylon/S_pyogenes/raw/mash/S_pyogenes_distances.tsv

set -euo pipefail

GENOME_DIR="/mnt/craig/pan_phylon/S_pyogenes/raw/mash_genomes"
OUT_DIR="/mnt/craig/pan_phylon/S_pyogenes/raw/mash"
PREFIX="S_pyogenes"
THREADS=16
KMER=21       # Mash default; good for species-level discrimination
SKETCH=1000   # Mash default sketch size
SKIP_SKETCH=0

# Parse flags
while [[ $# -gt 0 ]]; do
    case "$1" in
        --threads)      THREADS="$2";    shift 2 ;;
        --skip-sketch)  SKIP_SKETCH=1;   shift   ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
done

mkdir -p "$OUT_DIR"

SKETCH_FILE="$OUT_DIR/${PREFIX}"   # mash sketch adds .msh automatically
DIST_FILE="/mnt/craig/pan_phylon/S_pyogenes/raw/${PREFIX}_distances.tsv"
VALID_GENOMES_FILE="$OUT_DIR/valid_genome_list.txt"

echo "Genome dir  : $GENOME_DIR"
echo "Output dir  : $OUT_DIR"
echo "Threads     : $THREADS"
echo "k-mer size  : $KMER"
echo "Sketch size : $SKETCH"
echo "Genomes     : $(ls "$GENOME_DIR"/*.fna | wc -l)"
echo ""

# ── Pre-filter: remove empty / zero-byte FNA files ───────────────────────────
echo "── Pre-filter: checking for empty FNA files ──────"
n_empty=0
for f in "$GENOME_DIR"/*.fna; do
    if [[ ! -s "$f" ]]; then
        echo "  EMPTY (removing): $f"
        rm "$f"
        (( n_empty++ )) || true
    fi
done
echo "  Removed $n_empty empty file(s)"
echo "  Valid genomes: $(ls "$GENOME_DIR"/*.fna | wc -l)"
echo ""

# ── Step 1: Sketch ────────────────────────────────────────────────────────────
if [[ "$SKIP_SKETCH" -eq 1 ]]; then
    echo "── Step 1: Sketching (SKIPPED) ───────────────────"
    echo "  Using existing: ${SKETCH_FILE}.msh"
else
    echo "── Step 1: Sketching ─────────────────────────────"
    mash sketch \
        -o "$SKETCH_FILE" \
        -k "$KMER" \
        -s "$SKETCH" \
        -p "$THREADS" \
        "$GENOME_DIR"/*.fna
    echo "  Sketch written to: ${SKETCH_FILE}.msh"
fi
echo ""

# ── Step 2: Distance ──────────────────────────────────────────────────────────
echo "── Step 2: All-vs-all distances ──────────────────"
mash dist \
    -p "$THREADS" \
    "${SKETCH_FILE}.msh" \
    "${SKETCH_FILE}.msh" \
    > "$DIST_FILE"

echo "  Distances written to: $DIST_FILE"
echo "  Rows: $(wc -l < "$DIST_FILE")"
echo ""
echo "Done."

