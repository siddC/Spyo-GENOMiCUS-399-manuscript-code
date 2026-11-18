#!/usr/bin/env python3
"""
filter_ncbi_metadata.py

Filters the NCBI metadata TSV produced by fetch_ncbi_metadata.py using
assembly quality criteria appropriate for pangenome analysis.

Filters applied:
  1. Assembly completeness : num_contigs == 1  AND  contig_n50 >= ref_n50 * n50_fraction
  2. CheckM quality        : completeness >= min_completeness  AND  contamination <= max_contamination
                             (genomes with missing CheckM data are SKIPPED by default;
                              use --require-checkm to exclude them instead)
  3. GC content            : gc_min <= gc_percent <= gc_max

Usage:
    python filter_ncbi_metadata.py [options] INPUT_TSV

Examples:
    # S. pyogenes defaults (reference N50 ~1.7 Mb)
    python filter_ncbi_metadata.py ncbi_metadata_1314.tsv

    # Exclude genomes with missing CheckM data
    python filter_ncbi_metadata.py ncbi_metadata_1314.tsv --require-checkm

    # Override reference N50 for a different species
    python filter_ncbi_metadata.py ncbi_metadata_1280.tsv --ref-n50 2.8e6 --gc-min 32 --gc-max 36
"""

import argparse
import csv
import math
import sys
from pathlib import Path


# ── Parsing helpers ────────────────────────────────────────────────────────────

def parse_float(val: str) -> float | None:
    """Return float or None if missing/unparseable."""
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


def parse_int(val: str) -> int | None:
    """Return int or None if missing/unparseable."""
    try:
        return int(val)
    except (ValueError, TypeError):
        return None


# ── Filtering logic ────────────────────────────────────────────────────────────

def passes_filters(
    row: dict,
    ref_n50: float,
    n50_fraction: float,
    min_completeness: float,
    max_contamination: float,
    gc_min: float,
    gc_max: float,
    require_checkm: bool,
) -> tuple[bool, list[str]]:
    """Returns (passed, reasons_failed). reasons_failed is empty if passed."""
    failures = []

    # ── Filter 1: Assembly completeness ───────────────────────────────────────
    # contig_l50 == 1 means a single contig covers >=50% of the assembly.
    # This correctly passes strains with a chromosome + plasmid(s), because
    # the chromosome dominates the length. num_contigs == 1 would be WRONG:
    # it would incorrectly exclude all plasmid-bearing strains.
    n50_threshold = ref_n50 * n50_fraction

    contig_l50 = parse_int(row.get("contig_l50", ""))
    if contig_l50 is None:
        failures.append("contig_l50 missing (re-run fetch_ncbi_metadata.py to populate)")
    elif contig_l50 != 1:
        failures.append(f"contig_l50={contig_l50} (expected 1)")

    contig_n50 = parse_int(row.get("contig_n50", ""))
    if contig_n50 is None:
        failures.append("contig_n50 missing")
    elif contig_n50 < n50_threshold:
        failures.append(
            f"contig_n50={contig_n50:,} < {n50_threshold:,.0f} "
            f"({n50_fraction*100:.0f}% of ref {ref_n50:.2e})"
        )

    # ── Filter 2: CheckM quality ───────────────────────────────────────────────
    completeness  = parse_float(row.get("checkm_completeness", ""))
    contamination = parse_float(row.get("checkm_contamination", ""))

    if completeness is None or contamination is None:
        if require_checkm:
            failures.append("checkm_data missing")
        # else: silently pass — missing CheckM is treated as unknown, not failed
    else:
        if completeness < min_completeness:
            failures.append(
                f"checkm_completeness={completeness:.2f}% < {min_completeness}%"
            )
        if contamination > max_contamination:
            failures.append(
                f"checkm_contamination={contamination:.2f}% > {max_contamination}%"
            )

    # ── Filter 3: GC content ───────────────────────────────────────────────────
    gc = parse_float(row.get("gc_percent", ""))
    if gc is None:
        failures.append("gc_percent missing")
    elif not (gc_min <= gc <= gc_max):
        failures.append(f"gc_percent={gc:.1f}% outside [{gc_min}, {gc_max}]")

    return (len(failures) == 0, failures)


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input_tsv",
                        help="Input TSV from fetch_ncbi_metadata.py")
    parser.add_argument("--outfile", default=None,
                        help="Output TSV path (default: <input>_filtered.tsv)")
    parser.add_argument("--ref-n50", type=float, default=1.7e6,
                        help="Reference genome N50 in bp (default: 1.7e6 for S. pyogenes)")
    parser.add_argument("--n50-fraction", type=float, default=0.85,
                        help="Minimum N50 as fraction of ref_n50 (default: 0.85)")
    parser.add_argument("--min-completeness", type=float, default=97.0,
                        help="Minimum CheckM completeness %% (default: 97.0)")
    parser.add_argument("--max-contamination", type=float, default=3.0,
                        help="Maximum CheckM contamination %% (default: 3.0)")
    parser.add_argument("--gc-min", type=float, default=30.0,
                        help="Minimum GC content %% (default: 30.0)")
    parser.add_argument("--gc-max", type=float, default=40.0,
                        help="Maximum GC content %% (default: 40.0)")
    parser.add_argument("--require-checkm", action="store_true",
                        help="Exclude genomes with missing CheckM data (default: pass them through)")
    parser.add_argument("--verbose", action="store_true",
                        help="Print failure reasons for every excluded genome")
    args = parser.parse_args()

    inpath  = Path(args.input_tsv)
    outpath = Path(args.outfile) if args.outfile \
              else inpath.parent / (inpath.stem + "_filtered.tsv")

    n50_threshold = args.ref_n50 * args.n50_fraction

    print(f"Input  : {inpath}")
    print(f"Output : {outpath}")
    print(f"\nFilter thresholds:")
    print(f"  contig_l50           == 1")
    print(f"  contig_n50           >= {n50_threshold:,.0f} bp  ({args.n50_fraction*100:.0f}% of ref {args.ref_n50:.2e})")
    print(f"  checkm_completeness  >= {args.min_completeness}%  (missing = {'FAIL' if args.require_checkm else 'PASS'})")
    print(f"  checkm_contamination <= {args.max_contamination}%  (missing = {'FAIL' if args.require_checkm else 'PASS'})")
    print(f"  gc_percent           in [{args.gc_min}, {args.gc_max}]")
    print()

    failure_counts: dict[str, int] = {}
    checkm_missing = 0
    n_total = 0
    n_passed = 0

    with open(inpath, newline="", encoding="utf-8") as fin, \
         open(outpath, "w", newline="", encoding="utf-8") as fout:

        reader = csv.DictReader(fin, delimiter="\t")
        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames,
                                delimiter="\t", extrasaction="ignore")
        writer.writeheader()

        for row in reader:
            n_total += 1

            # Track how many are missing CheckM regardless of filter outcome
            if not row.get("checkm_completeness") or not row.get("checkm_contamination"):
                checkm_missing += 1

            passed, failures = passes_filters(
                row,
                ref_n50=args.ref_n50,
                n50_fraction=args.n50_fraction,
                min_completeness=args.min_completeness,
                max_contamination=args.max_contamination,
                gc_min=args.gc_min,
                gc_max=args.gc_max,
                require_checkm=args.require_checkm,
            )

            if passed:
                writer.writerow(row)
                n_passed += 1
            else:
                for reason in failures:
                    key = reason.split("=")[0].split(" ")[0]
                    failure_counts[key] = failure_counts.get(key, 0) + 1
                if args.verbose:
                    print(f"  FAIL {row.get('accession','?')}: {'; '.join(failures)}")

    n_failed = n_total - n_passed
    print(f"CheckM data present : {n_total - checkm_missing} / {n_total} genomes")
    print(f"Results             : {n_passed} / {n_total} passed ({n_failed} excluded)\n")

    if failure_counts:
        print("Exclusion reasons (genomes can fail multiple filters):")
        for reason, count in sorted(failure_counts.items(), key=lambda x: -x[1]):
            print(f"  {reason:<30} {count}")

    print(f"\nFiltered TSV saved to: {outpath}")


if __name__ == "__main__":
    main()

