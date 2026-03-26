#!/usr/bin/env bash
set -euo pipefail

UNFILTERED_DIR="/Users/cossa/Desktop/projects/manas_heart/results/unfiltered"
OUTPUT_DIR="/Users/cossa/Desktop/projects/manas_heart/results/filtered.1"

mkdir -p "$OUTPUT_DIR"

for vcf in "$UNFILTERED_DIR"/*_unfiltered.vcf.gz; do
    chunk=$(basename "$vcf" _unfiltered.vcf.gz)
    echo "Processing $chunk ..."
    /Users/cossa/miniforge3/envs/MiTo/bin/python filter_mutect2.py \
        "$vcf" \
        --chunk "$chunk" \
        -o "${OUTPUT_DIR}/${chunk}_filtered.tsv"
done

echo "Done."
