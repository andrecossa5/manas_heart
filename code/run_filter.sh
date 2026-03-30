#!/usr/bin/env bash
set -euo pipefail

VCF_DIR="${1:-/lustre/scratch125/cellgen/behjati/ac87/manas_heart/results/vcfs}"
OUTPUT_DIR="${2:-/lustre/scratch125/cellgen/behjati/ac87/manas_heart/results/muts.1}"
BIN_DIR="/lustre/scratch125/cellgen/behjati/ac87/manas_heart/code/manas_heart/bin"
PYTHON="/lustre/scratch126/cellgen/behjati/ac87/miniforge3/envs/mutect/bin/python"

mkdir -p "$OUTPUT_DIR"

# Run filter_mutect2.py on each unfiltered VCF
for vcf in "$VCF_DIR"/*_unfiltered.vcf.gz; do
    chunk=$(basename "$vcf" _unfiltered.vcf.gz)
    echo "Processing $chunk ..."
    "$PYTHON" "$BIN_DIR/filter_mutect2.py" \
        "$vcf" \
        --chunk "$chunk" \
        -o "${OUTPUT_DIR}/${chunk}_filtered.tsv"
done

# Gather all TSVs into ALL_FILTERED.tsv.gz
echo "Gathering mutation tables..."
"$PYTHON" "$BIN_DIR/gather_tables.py" \
    --mode muts \
    --output "${OUTPUT_DIR}/ALL_FILTERED.tsv.gz" \
    "${OUTPUT_DIR}"/*.tsv

# Gather all stats into ALL_FILTERED.stats
echo "Gathering stats..."
"$PYTHON" "$BIN_DIR/gather_tables.py" \
    --mode stats \
    --output "${OUTPUT_DIR}/ALL_FILTERED.stats" \
    "${OUTPUT_DIR}"/*.stats

echo "Done."
