#!/usr/bin/env python3
"""
fetch_ncbi_metadata.py

Fetches genome metadata for a given NCBI taxon from the NCBI Datasets API
and writes a flat TSV suitable for prefiltration (CheckM, GC, length, CDS).

Usage:
    python fetch_ncbi_metadata.py [--taxon TAXON_ID] [--assembly-level LEVEL] [--outfile FILE]

Examples:
    # All complete S. pyogenes genomes (default)
    python fetch_ncbi_metadata.py

    # All assembly levels (for broader survey)
    python fetch_ncbi_metadata.py --assembly-level all

    # Different taxon
    python fetch_ncbi_metadata.py --taxon 1280 --outfile saureus_metadata.tsv
"""

import argparse
import csv
import json
import sys
import time
import urllib.request
import urllib.parse
from typing import Any, Iterator

# ── Config ────────────────────────────────────────────────────────────────────

NCBI_API_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/taxon"
PAGE_SIZE = 200  # max allowed by NCBI API
RETRY_LIMIT = 3
RETRY_DELAY = 5  # seconds

# Fields to extract into TSV columns
# (column_name, dot-path into the JSON report object)
FIELDS = [
    # Identity
    ("accession",                   "accession"),
    ("paired_accession",            "paired_accession"),       # GCF <-> GCA counterpart
    ("organism_name",               "organism.organism_name"),
    ("tax_id",                      "organism.tax_id"),
    ("strain",                      "organism.infraspecific_names.strain"),

    # Assembly info
    ("assembly_level",              "assembly_info.assembly_level"),
    ("assembly_status",             "assembly_info.assembly_status"),
    ("assembly_name",               "assembly_info.assembly_name"),
    ("release_date",                "assembly_info.release_date"),
    ("sequencing_tech",             "assembly_info.sequencing_tech"),
    ("bioproject",                  "assembly_info.bioproject_accession"),
    ("biosample",                   "assembly_info.biosample.accession"),

    # Assembly stats  ← your prefiltration fields
    ("genome_length",               "assembly_stats.total_sequence_length"),
    ("num_contigs",                 "assembly_stats.number_of_contigs"),
    ("num_scaffolds",               "assembly_stats.number_of_scaffolds"),
    ("contig_n50",                  "assembly_stats.contig_n50"),
    ("contig_l50",                  "assembly_stats.contig_l50"),
    ("gc_percent",                  "assembly_stats.gc_percent"),
    ("gc_count",                    "assembly_stats.gc_count"),
    ("num_chromosomes",             "assembly_stats.total_number_of_chromosomes"),

    # Annotation stats  ← CDS count for length-vs-CDS filter
    ("total_genes",                 "annotation_info.stats.gene_counts.total"),
    ("protein_coding_genes",        "annotation_info.stats.gene_counts.protein_coding"),
    ("non_coding_genes",            "annotation_info.stats.gene_counts.non_coding"),

    # CheckM  ← completeness / contamination filters
    ("checkm_completeness",         "checkm_info.completeness"),
    ("checkm_contamination",        "checkm_info.contamination"),
    ("checkm_completeness_pct",     "checkm_info.completeness_percentile"),
    ("checkm_version",              "checkm_info.checkm_version"),
    ("checkm_marker_set",           "checkm_info.checkm_marker_set"),
    ("checkm_marker_set_rank",      "checkm_info.checkm_marker_set_rank"),

    # ANI taxonomy check  ← bonus: catch misidentified strains
    ("ani_tax_check_status",        "average_nucleotide_identity.taxonomy_check_status"),
    ("ani_match_status",            "average_nucleotide_identity.match_status"),
    ("ani_best_match_assembly",     "average_nucleotide_identity.best_ani_match.assembly"),
    ("ani_best_match_organism",     "average_nucleotide_identity.best_ani_match.organism_name"),
    ("ani_best_ani",                "average_nucleotide_identity.best_ani_match.ani"),

    # BioSample metadata  ← isolation source, geography, host
    ("collection_date",             "assembly_info.biosample.collection_date"),
    ("geo_loc_name",                "assembly_info.biosample.geo_loc_name"),
    ("host",                        "assembly_info.biosample.host"),
    ("isolation_source",            "assembly_info.biosample.isolation_source"),
]


# ── JSON path extractor ───────────────────────────────────────────────────────

def get_nested(obj: Any, path: str, default: str = "") -> str:
    """Traverse a dot-separated path through a nested dict/list."""
    parts = path.split(".")
    for part in parts:
        if isinstance(obj, dict):
            obj = obj.get(part)
        else:
            return default
        if obj is None:
            return default
    return str(obj) if obj is not None else default


# ── NCBI paginated fetch ──────────────────────────────────────────────────────

def build_url(taxon: str, assembly_level: str, page_token: str = "") -> str:
    params: dict[str, str] = {"page_size": str(PAGE_SIZE)}
    if assembly_level != "all":
        params["filters.assembly_level"] = assembly_level
    if page_token:
        params["page_token"] = page_token
    query = urllib.parse.urlencode(params)
    return f"{NCBI_API_BASE}/{taxon}/dataset_report?{query}"


def fetch_page(url: str) -> dict:
    for attempt in range(1, RETRY_LIMIT + 1):
        try:
            req = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.load(resp)
        except Exception as e:
            print(f"  Attempt {attempt}/{RETRY_LIMIT} failed: {e}", file=sys.stderr)
            if attempt < RETRY_LIMIT:
                time.sleep(RETRY_DELAY)
    raise RuntimeError(f"Failed to fetch {url} after {RETRY_LIMIT} attempts")


def iter_reports(taxon: str, assembly_level: str) -> Iterator[dict]:
    """Yield every genome report, handling pagination automatically."""
    page_token = ""
    page_num = 0
    total = None

    while True:
        url = build_url(taxon, assembly_level, page_token)
        data = fetch_page(url)

        if total is None:
            total = data.get("total_count", "?")
            print(f"  Total genomes reported by NCBI: {total}", file=sys.stderr)

        reports = data.get("reports", [])
        page_num += 1
        print(f"  Page {page_num}: {len(reports)} records fetched", file=sys.stderr)

        yield from reports

        page_token = data.get("next_page_token", "")
        if not page_token:
            break


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--taxon", default="1314",
                        help="NCBI taxon ID (default: 1314 = S. pyogenes)")
    parser.add_argument("--assembly-level", default="complete_genome",
                        choices=["complete_genome", "chromosome", "scaffold", "contig", "all"],
                        help="Assembly level filter (default: complete_genome)")
    parser.add_argument("--outfile", default=None,
                        help="Output TSV path (default: ncbi_metadata_<taxon>.tsv)")
    args = parser.parse_args()

    outfile = args.outfile or f"ncbi_metadata_{args.taxon}.tsv"
    col_names = [f[0] for f in FIELDS]

    print(f"Fetching NCBI metadata for taxon {args.taxon} "
          f"(assembly_level={args.assembly_level}) ...", file=sys.stderr)

    n_written = 0
    with open(outfile, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=col_names, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()

        for report in iter_reports(args.taxon, args.assembly_level):
            row = {col: get_nested(report, path) for col, path in FIELDS}
            writer.writerow(row)
            n_written += 1

    print(f"\nDone. Wrote {n_written} records to {outfile}", file=sys.stderr)


if __name__ == "__main__":
    main()

