#!/usr/bin/env python3
"""
download_ncbi_genomes.py

Downloads .fna.gz files from NCBI FTP for all accessions in a filtered
metadata TSV (produced by filter_ncbi_metadata.py), decompresses them,
and saves each as {accession}.fna in the output directory.

NCBI FTP path structure:
  https://ftp.ncbi.nlm.nih.gov/genomes/all/{db}/{NNN}/{NNN}/{NNN}/
      {accession}_{assembly_name}/{accession}_{assembly_name}_genomic.fna.gz

Usage:
    python download_ncbi_genomes.py INPUT_TSV OUTPUT_DIR

Example:
    python download_ncbi_genomes.py \\
        ncbi_metadata_S_pyogenes_filtered.tsv \\
        /mnt/craig/pan_phylon/S_pyogenes/raw/temp_genomes
"""

import argparse
import csv
import gzip
import re
import shutil
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from tqdm.auto import tqdm


# ── URL construction ───────────────────────────────────────────────────────────

def accession_to_ftp_url(accession: str, assembly_name: str) -> str:
    """
    Build the NCBI HTTPS-FTP URL for a genomic FASTA file.

    Example:
        GCA_900475035.1, 41965_F01
        -> https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/900/475/035/
               GCA_900475035.1_41965_F01/GCA_900475035.1_41965_F01_genomic.fna.gz
    """
    # Strip version suffix to get numeric part: GCA_900475035.1 -> 900475035
    match = re.match(r"(GC[AF])_(\d+)\.\d+", accession)
    if not match:
        raise ValueError(f"Cannot parse accession: {accession!r}")

    db, digits = match.group(1), match.group(2)
    digits_padded = digits.zfill(9)                          # zero-pad to 9 digits
    path_chunks = "/".join([digits_padded[i:i+3] for i in (0, 3, 6)])  # NNN/NNN/NNN

    asm_dir  = f"{accession}_{assembly_name}"
    filename = f"{asm_dir}_genomic.fna.gz"

    return (
        f"https://ftp.ncbi.nlm.nih.gov/genomes/all/{db}/"
        f"{path_chunks}/{asm_dir}/{filename}"
    )


# ── Download + decompress ──────────────────────────────────────────────────────

# ── Download + decompress ──────────────────────────────────────────────────────

MAX_RETRIES = 4
RETRY_BACKOFF = [5, 15, 30, 60]  # seconds between retries

def download_and_decompress(url: str, out_fna: Path, accession: str) -> bool:
    """
    Stream-download .fna.gz from NCBI and decompress directly to out_fna.
    No intermediate .gz file. Retries with backoff on failure.
    """
    import time
    import urllib.request

    headers = {"User-Agent": "Mozilla/5.0"}  # NCBI rejects blank user-agents

    for attempt in range(MAX_RETRIES + 1):
        if attempt > 0:
            wait = RETRY_BACKOFF[min(attempt - 1, len(RETRY_BACKOFF) - 1)]
            tqdm.write(f"  RETRY {attempt}/{MAX_RETRIES} for {accession} (waiting {wait}s)")
            time.sleep(wait)

        out_fna.unlink(missing_ok=True)

        try:
            req = urllib.request.Request(url, headers=headers)
            with urllib.request.urlopen(req, timeout=120) as resp, \
                 gzip.open(resp, "rb") as gz, \
                 open(out_fna, "wb") as fout:
                shutil.copyfileobj(gz, fout)

            if out_fna.exists() and out_fna.stat().st_size > 0:
                return True
            out_fna.unlink(missing_ok=True)

        except Exception as e:
            out_fna.unlink(missing_ok=True)
            tqdm.write(f"  attempt {attempt+1} failed for {accession}: {e}")

    return False


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input_tsv",
                        help="Filtered metadata TSV (from filter_ncbi_metadata.py)")
    parser.add_argument("output_dir",
                        help="Directory to save .fna files")
    parser.add_argument("--redownload", action="store_true",
                        help="Re-download files that already exist locally")
    parser.add_argument("--workers", type=int, default=4,
                        help="Number of parallel download threads (default: 4 — NCBI rate-limits aggressive clients)")
    args = parser.parse_args()

    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Read filtered TSV
    with open(args.input_tsv, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))

    total = len(rows)
    print(f"Genomes to download : {total}")
    print(f"Output directory    : {outdir}\n")

    n_ok      = 0
    n_skipped = 0
    n_failed  = 0
    failed_accessions = []

    # Build work list (skip already-downloaded)
    work = []
    for row in rows:
        accession     = row.get("accession", "").strip()
        assembly_name = row.get("assembly_name", "").strip()
        if not accession or not assembly_name:
            tqdm.write(f"SKIP — missing accession or assembly_name: {row}")
            n_failed += 1
            continue
        out_fna = outdir / f"{accession}.fna"
        if out_fna.exists() and not args.redownload:
            n_skipped += 1
            continue
        try:
            url = accession_to_ftp_url(accession, assembly_name)
        except ValueError as e:
            tqdm.write(f"ERROR  {e}")
            n_failed += 1
            failed_accessions.append(accession)
            continue
        work.append((accession, url, out_fna))

    tqdm.write(f"To download: {len(work)}  |  Already exists: {n_skipped}  |  Workers: {args.workers}\n")

    def _download_one(item):
        accession, url, out_fna = item
        success = download_and_decompress(url, out_fna, accession)
        return accession, out_fna, success

    with tqdm(total=len(work), desc="Downloading genomes", unit="genome") as pbar:
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = {pool.submit(_download_one, item): item for item in work}
            for future in as_completed(futures):
                accession, out_fna, success = future.result()
                if success:
                    size_mb = out_fna.stat().st_size / 1e6
                    tqdm.write(f"OK     {accession}  ({size_mb:.1f} MB)")
                    n_ok += 1
                else:
                    tqdm.write(f"FAILED {accession}")
                    n_failed += 1
                    failed_accessions.append(accession)
                pbar.set_postfix(ok=n_ok, fail=n_failed)
                pbar.update(1)

    # Summary
    print(f"\n{'─'*50}")
    print(f"Downloaded : {n_ok}")
    print(f"Skipped    : {n_skipped}  (already existed)")
    print(f"Failed     : {n_failed}")

    if failed_accessions:
        fail_log = outdir / "failed_accessions.txt"
        fail_log.write_text("\n".join(failed_accessions) + "\n")
        print(f"\nFailed accessions written to: {fail_log}")


if __name__ == "__main__":
    main()

