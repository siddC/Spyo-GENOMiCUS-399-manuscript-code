#!/bin/bash
# copy_gff3.sh
# Copy all Bakta GFF3 files into a single flat directory for PanTA

SRC="/mnt/craig/pan_phylon/S_pyogenes/processed/bakta"
DST="/mnt/craig/pan_phylon/S_pyogenes/processed/gff3_files"

rm -rf "$DST"
mkdir -p "$DST"

count=0
for gff in "$SRC"/*/*.gff3; do
    cp "$gff" "$DST/"
    ((count++))
done

echo "Copied $count GFF3 files to $DST"

