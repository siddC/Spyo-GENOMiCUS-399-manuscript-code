#!/usr/bin/env python3
"""
download_bvbrc_metadata.py

Downloads genome_summary, genome_metadata, and genome_lineage from the
BV-BRC FTPS server (ftps://ftp.bvbrc.org/RELEASE_NOTES/) using wget.

Usage:
    python download_bvbrc_metadata.py [--outdir OUTPUT_DIR]

Example:
    python download_bvbrc_metadata.py --outdir /mnt/craig/pan_phylon/bv_brc_files
"""

import argparse
import shutil
import subprocess
import sys
from pathlib import Path


# ── Config ────────────────────────────────────────────────────────────────────

FTP_BASE = "ftps://ftp.bvbrc.org/RELEASE_NOTES"

FILES = [
    "genome_summary",
    "genome_metadata",
    "genome_lineage",
]


# ── Helpers ───────────────────────────────────────────────────────────────────

def check_wget() -> None:
    if shutil.which("wget") is None:
        sys.exit(
            "ERROR: wget not found. Install it with:\n"
            "  conda install -c conda-forge wget\n"
            "  # or: sudo apt install wget"
        )


def download_file(url: str, dest: Path) -> None:
    """
    Download a single file via wget with FTPS (explicit TLS).
    --continue resumes partial downloads if re-run.
    """
    cmd = [
        "wget",
        "--ftp-user=anonymous",
        "--ftp-password=guest",
        "--secure-protocol=auto",
        "--no-check-certificate",
        "--continue",           # resume if partially downloaded
        "--show-progress",
        "-O", str(dest),
        url,
    ]
    print(f"-> {dest.name}")
    print(f"   {url}")
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print(f"   ERROR: wget exited with code {result.returncode}")
    else:
        size_mb = dest.stat().st_size / 1e6
        print(f"   Saved: {dest}  ({size_mb:.1f} MB)")
    print()


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--outdir",
        default="bvbrc_metadata",
        help="Directory to save downloaded files (default: bvbrc_metadata/)",
    )
    args = parser.parse_args()

    check_wget()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Downloading {len(FILES)} files to {outdir}/\n")

    for fname in FILES:
        url  = f"{FTP_BASE}/{fname}"
        dest = outdir / fname
        download_file(url, dest)

    print("Done.")


if __name__ == "__main__":
    main()