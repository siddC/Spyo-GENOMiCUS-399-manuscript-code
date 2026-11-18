#!/usr/bin/env python3
"""
merge_genome_metadata.py

Deduplicates and merges BV-BRC and NCBI genome collections into a single
unified metadata table for S. pyogenes pangenome analysis.

Steps:
  1. Load BV-BRC metadata + summary, filter by CheckM thresholds
  2. Load NCBI CheckM results, keep PASS only, deduplicate GCA/GCF pairs
  3. Deduplicate between BV-BRC and NCBI via assembly accession
  4. Merge BV-BRC lineage onto BV-BRC genomes
  5. Write unified TSV with source column

Usage:
    python merge_genome_metadata.py [options]
"""

import argparse
import re
import sys

import pandas as pd


# ── Defaults ───────────────────────────────────────────────────────────────────

DEFAULTS = dict(
    bvbrc_metadata   = "/mnt/craig/pan_phylon/S_pyogenes/metadata/scrubbed_species_metadata.csv",
    bvbrc_summary    = "/mnt/craig/pan_phylon/S_pyogenes/metadata/scrubbed_species_summary.csv",
    bvbrc_lineage    = "/mnt/craig/pan_phylon/bv_brc_files/genome_lineage_04_02_26",
    ncbi_checkm      = "/mnt/craig/pan_phylon/S_pyogenes/raw/checkm_output/checkm_qa_filtered.tsv",
    ncbi_metadata    = "/home/sidd/jupyterlab/S_pyogenes/1_obtain_data/ncbi/ncbi_metadata_S_pyogenes_filtered.tsv",
    outfile          = "/mnt/craig/pan_phylon/S_pyogenes/metadata/unified_genome_metadata.csv",
    bvbrc_genome_dir = "/mnt/craig/pan_phylon/S_pyogenes/raw/genomes",
    ncbi_genome_dir  = "/mnt/craig/pan_phylon/S_pyogenes/raw/temp_genomes",
    min_completeness = 97.0,
    max_contamination = 3.0,
)

LINEAGE_COLS = ["kingdom", "phylum", "class", "order", "family", "genus", "species",
                "taxon_lineage_ids", "taxon_lineage_names"]


# ── Helpers ────────────────────────────────────────────────────────────────────

def log(msg: str) -> None:
    print(msg, flush=True)


def extract_accession(val: str) -> str | None:
    """
    Pull a GCA_/GCF_ accession out of a BV-BRC accession string.
    BV-BRC stores these as e.g. 'GCA_000006785.2' or comma-separated lists.
    Returns the first GCA match (strips version), or None.
    """
    if not isinstance(val, str):
        return None
    match = re.search(r"(GC[AF]_\d+\.\d+)", val)
    return match.group(1) if match else None


def strip_version(acc: str) -> str:
    """GCA_000006785.2 -> GCA_000006785"""
    return re.sub(r"\.\d+$", "", acc)


def accession_base(acc: str) -> str:
    """GCA_000006785.2 -> 000006785  (digits only, for GCA/GCF cross-matching)"""
    m = re.search(r"GC[AF]_(\d+)", acc)
    return m.group(1) if m else acc


# ── Step 1: Load + filter BV-BRC ──────────────────────────────────────────────

def load_bvbrc(meta_path, summary_path, min_comp, max_cont):
    log("\n── BV-BRC ────────────────────────────────────────")

    meta    = pd.read_csv(meta_path, index_col=0, low_memory=False)
    summary = pd.read_csv(summary_path, index_col=0, low_memory=False)
    log(f"  Metadata : {len(meta)} genomes × {meta.shape[1]} columns")
    log(f"  Summary  : {len(summary)} genomes × {summary.shape[1]} columns")

    # Merge summary checkm columns onto metadata (genome_id is the shared key)
    checkm_cols = ["genome_id", "checkm_completeness", "checkm_contamination",
                   "genome_length", "gc_content", "contig_l50", "contig_n50",
                   "chromosomes", "plasmids", "contigs", "patric_cds", "refseq_cds"]
    available   = [c for c in checkm_cols if c in summary.columns]
    meta = meta.merge(summary[available], on="genome_id", how="left", suffixes=("", "_summary"))

    # CheckM filter — require data to be present AND pass thresholds.
    # NaNs are dropped (not passed through) so every genome in the final
    # pangenome can be audited as having explicitly passed CheckM QC.
    before    = len(meta)
    n_missing = (meta["checkm_completeness"].isna() | meta["checkm_contamination"].isna()).sum()
    if n_missing > 0:
        log(f"  WARNING: {n_missing} BV-BRC genome(s) have missing CheckM data — dropping")
    meta = meta.dropna(subset=["checkm_completeness", "checkm_contamination"])
    meta = meta[
        (meta["checkm_completeness"] >= min_comp) &
        (meta["checkm_contamination"] <= max_cont)
    ]
    log(f"  After CheckM filter (≥{min_comp}% complete, ≤{max_cont}% contam): "
        f"{len(meta)} / {before}  ({n_missing} dropped for missing data)")

    # Extract a normalised assembly accession for dedup
    # BV-BRC stores GCA in genbank_accessions, GCF in refseq_accessions
    def get_acc(row):
        for col in ("genbank_accessions", "refseq_accessions", "assembly_accession"):
            acc = extract_accession(str(row.get(col, "")))
            if acc:
                return acc
        return None

    meta["_acc"] = meta.apply(get_acc, axis=1)
    meta["_acc_base"] = meta["_acc"].apply(
        lambda x: accession_base(x) if isinstance(x, str) else None
    )
    meta["source"] = "BV-BRC"
    log(f"  Accessions resolved: {meta['_acc'].notna().sum()} / {len(meta)}")
    return meta


# ── Step 2: Load + filter NCBI CheckM ─────────────────────────────────────────

def load_ncbi_checkm(checkm_path, ncbi_metadata_path, min_comp, max_cont):
    log("\n── NCBI CheckM ───────────────────────────────────")

    df = pd.read_csv(checkm_path, sep="\t", low_memory=False)
    log(f"  Loaded   : {len(df)} rows")

    # Keep only PASS (already filtered, but be explicit)
    df = df[df["checkm_pass"] == "PASS"].copy()
    log(f"  PASS     : {len(df)}")

    # Rename key columns to snake_case for the unified table
    df = df.rename(columns={
        "Bin Id"              : "accession",
        "Completeness"        : "checkm_completeness",
        "Contamination"       : "checkm_contamination",
        "Strain heterogeneity": "checkm_strain_heterogeneity",
        "Genome size (bp)"    : "genome_length",
        "# contigs"           : "contigs",
        "# scaffolds"         : "scaffolds",
        "N50 (contigs)"       : "contig_n50",
        "N50 (scaffolds)"     : "scaffold_n50",
        "GC"                  : "gc_content",
        "# predicted genes"   : "predicted_genes",
        "Marker lineage"      : "checkm_marker_lineage",
        "# markers"           : "checkm_markers",
        "# marker sets"       : "checkm_marker_sets",
    })

    df["_acc"]      = df["accession"]
    df["_acc_base"] = df["accession"].apply(accession_base)

    # ── Deduplicate GCA/GCF pairs within NCBI ─────────────────────────────────
    # Both GCA and GCF for the same genome share the same digit block.
    # Prefer GCF (RefSeq, curated annotation) over GCA when both present.
    df["_is_gcf"] = df["accession"].str.startswith("GCF")
    df = df.sort_values(["_acc_base", "_is_gcf"], ascending=[True, False])
    df = df.drop_duplicates(subset="_acc_base", keep="first")
    log(f"  After GCA/GCF dedup : {len(df)}")

    # ── Join NCBI assembly metadata ───────────────────────────────────────────
    # Brings in assembly_level (-> genome_status), strain, biosample, geo,
    # isolation_source, host, collection_date, and ANI taxonomy check fields.
    ncbi_meta = pd.read_csv(ncbi_metadata_path, sep="\t", low_memory=False)
    log(f"  NCBI metadata rows: {len(ncbi_meta)}")

    # Normalise accession for joining — match on base digits (GCA/GCF agnostic)
    ncbi_meta["_acc_base"] = ncbi_meta["accession"].apply(accession_base)

    # Deduplicate ncbi_meta on _acc_base (same GCA/GCF logic: prefer GCF)
    ncbi_meta["_is_gcf"] = ncbi_meta["accession"].str.startswith("GCF")
    ncbi_meta = ncbi_meta.sort_values(["_acc_base", "_is_gcf"], ascending=[True, False])
    ncbi_meta = ncbi_meta.drop_duplicates(subset="_acc_base", keep="first")
    ncbi_meta = ncbi_meta.drop(columns=["_is_gcf"])

    # Rename NCBI metadata columns to match BV-BRC naming conventions
    ncbi_meta_rename = {
        "assembly_level"      : "genome_status",
        "biosample"           : "biosample_accession",
        "bioproject"          : "bioproject_accession",
        "host"                : "host_name",
        "geo_loc_name"        : "geographic_location",
        "gc_percent"          : "gc_content",
        "num_contigs"         : "contigs",
        "num_chromosomes"     : "chromosomes",
        "tax_id"              : "taxon_id",
        "sequencing_tech"     : "sequencing_platform",
        "paired_accession"    : "refseq_accessions",
        # protein_coding_genes intentionally excluded — patric_cds is already
        # populated from the CheckM df (renamed from predicted_genes)
    }
    ncbi_meta = ncbi_meta.rename(columns={
        k: v for k, v in ncbi_meta_rename.items() if k in ncbi_meta.columns
    })

    # Bring in all columns not already present in the checkm df,
    # checking AFTER renaming to avoid post-rename clashes.
    already_in_checkm = set(df.columns) - {"_acc_base"}
    meta_cols = ["_acc_base"] + [
        c for c in ncbi_meta.columns
        if c not in already_in_checkm and c != "_acc_base"
    ]

    df = df.merge(ncbi_meta[meta_cols], on="_acc_base", how="left")

    # Rename assembly_level -> genome_status to match BV-BRC convention
    if "assembly_level" in df.columns:
        df = df.rename(columns={"assembly_level": "genome_status"})

    log(f"  genome_status populated: {df['genome_status'].notna().sum()} / {len(df)}")

    df["source"] = "NCBI"
    return df


# ── Step 3: Deduplicate between BV-BRC and NCBI ───────────────────────────────

def dedup_between(bvbrc: pd.DataFrame, ncbi: pd.DataFrame):
    log("\n── Deduplication ─────────────────────────────────")

    bvbrc_bases = set(bvbrc["_acc_base"].dropna())
    ncbi["_in_bvbrc"] = ncbi["_acc_base"].isin(bvbrc_bases)

    ncbi_only = ncbi[~ncbi["_in_bvbrc"]].copy()
    overlap   = ncbi["_in_bvbrc"].sum()

    log(f"  BV-BRC genomes              : {len(bvbrc)}")
    log(f"  NCBI genomes (post-dedup)   : {len(ncbi)}")
    log(f"  Overlap (in both)           : {overlap}  → keeping BV-BRC version")
    log(f"  NCBI-only (novel)           : {len(ncbi_only)}")
    log(f"  Total after merge           : {len(bvbrc) + len(ncbi_only)}")

    return bvbrc, ncbi_only


# ── Step 4: Attach lineage to BV-BRC genomes ──────────────────────────────────

def attach_lineage(bvbrc: pd.DataFrame, lineage_path: str) -> pd.DataFrame:
    log("\n── Lineage ───────────────────────────────────────")

    lineage = pd.read_csv(lineage_path, sep="\t", low_memory=False)
    log(f"  Loaded lineage: {len(lineage)} rows")

    # Filter to S. pyogenes only before merging (speeds up merge significantly)
    spyo = lineage[lineage["species"].str.contains("pyogenes", na=False)]
    log(f"  S. pyogenes rows: {len(spyo)}")

    # Cast genome_id to str on both sides — BV-BRC reads it as float (1010840.4),
    # lineage file stores it as string; mismatched dtypes silently produce 0 matches.
    bvbrc["genome_id"] = bvbrc["genome_id"].astype(str)
    spyo = spyo.copy()
    spyo["genome_id"]  = spyo["genome_id"].astype(str)

    # Deduplicate lineage on genome_id — some genomes appear multiple times
    # (e.g. different taxon entries), which would inflate row count after merge.
    spyo = spyo.drop_duplicates(subset="genome_id", keep="first")
    log(f"  S. pyogenes rows (after dedup): {len(spyo)}")

    # Note: 'kingdom' is NaN for all bacterial entries in BV-BRC lineage;
    # use 'species' as the success check instead.
    cols = ["genome_id"] + [c for c in LINEAGE_COLS if c in spyo.columns]
    bvbrc = bvbrc.merge(spyo[cols], on="genome_id", how="left")
    log(f"  Lineage attached to {bvbrc['species'].notna().sum()} / {len(bvbrc)} BV-BRC genomes")
    return bvbrc



# ── Step 5: Validate FNA files exist on disk ───────────────────────────────────

def validate_fna_on_disk(
    bvbrc: pd.DataFrame,
    ncbi_only: pd.DataFrame,
    bvbrc_genome_dir: str,
    ncbi_genome_dir: str,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Drop any genome whose .fna file cannot be found on disk.
    BV-BRC: {bvbrc_genome_dir}/{genome_id}/{genome_id}.fna
    NCBI:   {ncbi_genome_dir}/{assembly_accession_resolved}.fna

    Phantom entries (deprecated PATRIC IDs present in metadata but removed
    from the FTP) are logged and excluded so the final pangenome only contains
    genomes with verified sequence data on disk.
    """
    from pathlib import Path

    log("\n── FNA validation ────────────────────────────────")

    bvbrc_dir = Path(bvbrc_genome_dir)
    ncbi_dir  = Path(ncbi_genome_dir)

    def bvbrc_exists(row) -> bool:
        gid = str(row.get("genome_id", "")).strip()
        return gid and (bvbrc_dir / gid / f"{gid}.fna").exists()

    def ncbi_exists(row) -> bool:
        acc = str(row.get("_acc", "")).strip()
        return acc and (ncbi_dir / f"{acc}.fna").exists()

    bvbrc_mask = bvbrc.apply(bvbrc_exists, axis=1)
    ncbi_mask  = ncbi_only.apply(ncbi_exists, axis=1)

    bvbrc_missing = bvbrc[~bvbrc_mask]
    ncbi_missing  = ncbi_only[~ncbi_mask]

    if len(bvbrc_missing) > 0:
        log(f"  BV-BRC phantom entries dropped ({len(bvbrc_missing)}):")
        for _, row in bvbrc_missing.iterrows():
            log(f"    {row.get('genome_id', '?')}  ({row.get('genome_name', '?')})")

    if len(ncbi_missing) > 0:
        log(f"  NCBI missing FNA dropped ({len(ncbi_missing)}):")
        for _, row in ncbi_missing.iterrows():
            log(f"    {row.get('_acc', '?')}")

    bvbrc_clean    = bvbrc[bvbrc_mask].copy()
    ncbi_only_clean = ncbi_only[ncbi_mask].copy()

    log(f"  BV-BRC : {len(bvbrc_clean)} retained  ({len(bvbrc_missing)} dropped)")
    log(f"  NCBI   : {len(ncbi_only_clean)} retained  ({len(ncbi_missing)} dropped)")

    return bvbrc_clean, ncbi_only_clean


# ── Step 6: Build unified table ────────────────────────────────────────────────

def build_unified(bvbrc: pd.DataFrame, ncbi_only: pd.DataFrame) -> pd.DataFrame:
    log("\n── Building unified table ────────────────────────")

    # ── BV-BRC: keep ALL columns, just rename internal helpers ────────────────
    # Drop internal dedup helper columns, rename _acc to assembly_accession_resolved
    bvbrc_out = bvbrc.copy()
    bvbrc_out = bvbrc_out.rename(columns={"_acc": "assembly_accession_resolved"})
    bvbrc_out = bvbrc_out.drop(columns=[c for c in ("_acc_base", "_is_gcf") if c in bvbrc_out.columns])

    # ── NCBI: rename columns to match BV-BRC naming conventions ──────────────
    ncbi_rename = {
        "accession"                  : "assembly_accession_resolved",
        "checkm_strain_heterogeneity": "checkm_strain_heterogeneity",
        "checkm_marker_lineage"      : "checkm_marker_lineage",
        "predicted_genes"            : "patric_cds",
    }
    ncbi_out = ncbi_only.copy()
    ncbi_out = ncbi_out.rename(columns={k: v for k, v in ncbi_rename.items() if k in ncbi_out.columns})
    ncbi_out = ncbi_out.drop(columns=[c for c in ("_acc", "_acc_base", "_is_gcf", "_in_bvbrc") if c in ncbi_out.columns])

    # ── Concatenate ───────────────────────────────────────────────────────────
    # Safety dedup — should never trigger, but guards against future regressions
    bvbrc_dupes = bvbrc_out.columns[bvbrc_out.columns.duplicated()].tolist()
    ncbi_dupes  = ncbi_out.columns[ncbi_out.columns.duplicated()].tolist()
    if bvbrc_dupes:
        log(f"  WARNING: duplicate columns in BV-BRC (dropping): {bvbrc_dupes}")
        bvbrc_out = bvbrc_out.loc[:, ~bvbrc_out.columns.duplicated()]
    if ncbi_dupes:
        log(f"  WARNING: duplicate columns in NCBI (dropping): {ncbi_dupes}")
        ncbi_out = ncbi_out.loc[:, ~ncbi_out.columns.duplicated()]

    unified = pd.concat([bvbrc_out, ncbi_out], ignore_index=True, sort=False)

    # NCBI rows have no genome_id (BV-BRC internal ID) — use the assembly
    # accession as the genome_id so every row has a unique, meaningful identifier.
    if "genome_id" in unified.columns:
        unified.loc[unified["source"] == "NCBI", "genome_id"] = \
            unified.loc[unified["source"] == "NCBI", "assembly_accession_resolved"]

    # ── Deduplicate on genome_id ──────────────────────────────────────────────
    # Some BV-BRC genomes appear under multiple assembly accessions (e.g. both
    # GCA and GCF) resulting in the same genome_id after merging. Keep one row
    # per genome_id, preferring GCF over GCA where both exist.
    n_before = len(unified)
    unified["_is_gcf"] = unified["assembly_accession_resolved"].str.startswith("GCF").fillna(False)
    unified = unified.sort_values("_is_gcf", ascending=False)  # GCF rows first
    unified = unified.drop_duplicates(subset="genome_id", keep="first")
    unified = unified.drop(columns=["_is_gcf"])
    n_dropped = n_before - len(unified)
    if n_dropped > 0:
        log(f"  Dropped {n_dropped} duplicate genome_id rows (kept GCF over GCA where ambiguous)")
    priority_cols = ["assembly_accession_resolved", "source", "genome_id",
                     "genome_name", "organism_name", "taxon_id", "strain"]
    front = [c for c in priority_cols if c in unified.columns]
    rest  = sorted([c for c in unified.columns if c not in front])
    unified = unified[front + rest]

    log(f"  Unified table: {len(unified)} genomes × {unified.shape[1]} columns")
    log(f"  BV-BRC : {(unified['source'] == 'BV-BRC').sum()}")
    log(f"  NCBI   : {(unified['source'] == 'NCBI').sum()}")
    return unified


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--bvbrc-metadata",   default=DEFAULTS["bvbrc_metadata"])
    parser.add_argument("--bvbrc-summary",    default=DEFAULTS["bvbrc_summary"])
    parser.add_argument("--bvbrc-lineage",    default=DEFAULTS["bvbrc_lineage"])
    parser.add_argument("--ncbi-checkm",      default=DEFAULTS["ncbi_checkm"])
    parser.add_argument("--ncbi-metadata",    default=DEFAULTS["ncbi_metadata"],
                        help="Filtered NCBI assembly metadata TSV (from fetch_ncbi_metadata.py)")
    parser.add_argument("--outfile",          default=DEFAULTS["outfile"])
    parser.add_argument("--bvbrc-genome-dir", default=DEFAULTS["bvbrc_genome_dir"],
                        help="Directory containing BV-BRC genomes (subdirs per genome_id)")
    parser.add_argument("--ncbi-genome-dir",  default=DEFAULTS["ncbi_genome_dir"],
                        help="Directory containing NCBI .fna files (flat)")
    parser.add_argument("--min-completeness", type=float, default=DEFAULTS["min_completeness"])
    parser.add_argument("--max-contamination",type=float, default=DEFAULTS["max_contamination"])
    args = parser.parse_args()

    log(f"CheckM thresholds: completeness ≥ {args.min_completeness}%, "
        f"contamination ≤ {args.max_contamination}%")

    bvbrc    = load_bvbrc(args.bvbrc_metadata, args.bvbrc_summary,
                          args.min_completeness, args.max_contamination)
    ncbi     = load_ncbi_checkm(args.ncbi_checkm, args.ncbi_metadata,
                                args.min_completeness, args.max_contamination)
    bvbrc, ncbi_only = dedup_between(bvbrc, ncbi)
    bvbrc    = attach_lineage(bvbrc, args.bvbrc_lineage)
    bvbrc, ncbi_only = validate_fna_on_disk(
        bvbrc, ncbi_only,
        args.bvbrc_genome_dir, args.ncbi_genome_dir,
    )
    unified  = build_unified(bvbrc, ncbi_only)

    unified.to_csv(args.outfile, sep=",", index=False)
    log(f"\n✓ Unified metadata written to: {args.outfile}")


if __name__ == "__main__":
    main()

