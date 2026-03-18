#!/usr/bin/env bash
##############################################################################
# create_symlinks_placenta.sh
#   For every sample listed in placenta_samples.csv, create symbolic links in
#   $TARGET_DIR pointing at the newest *.sample.dupmarked.bam file and its
#   matching *.bam.bai index.  Works for IDs with or without "_loNNNN" suffix.
##############################################################################
set -euo pipefail

# ──────────────────────────  paths ──────────────────────────
TARGET_DIR="/lustre/scratch125/cellgen/behjati/ac87/manas_heart/work/placenta"
SAMPLE_CSV="$(dirname "$0")/placenta_samples.csv"

mkdir -p "$TARGET_DIR"

# ─────────────────────────  projects ────────────────────────
projects=(2855 3254)

# ───────────  irods-CGP roots visible under /nfs  ───────────
mapfile -t irods_roots < <(find /nfs -maxdepth 1 -type d -name 'irods-cgp*' -print)

# ───────────  load unique sample IDs from the CSV (column 1, skip header)  ────
if [[ ! -r $SAMPLE_CSV ]]; then
  echo "ERROR: cannot read $SAMPLE_CSV" >&2
  exit 1
fi
mapfile -t samples < <(
  tail -n +2 "$SAMPLE_CSV" | cut -d',' -f1 | sort -u
)
printf 'Loaded %d samples\n' "${#samples[@]}"

# ─────────────────────────  helper function  ────────────────
link_latest () {
  local sample=$1
  local base=${sample%%_lo*}    # PDxxxxx   (without lane suffix)
  local alt=$sample             # PDxxxxx_loNNNN (if suffix present)
  local hits=()

  # search all irods roots / projects / candidate sample dirs
  for root in "${irods_roots[@]}"; do
    for proj in "${projects[@]}"; do
      for d in "$root/intproj/$proj/sample/$base" \
               "$root/intproj/$proj/sample/$alt"; do
        [[ -d $d ]] || continue
        shopt -s globstar nullglob
        hits+=( "$d"/**/"$sample".v*.sample.dupmarked.bam )
        hits+=( "$d"/**/"$sample".sample.dupmarked.bam )
        shopt -u globstar nullglob
      done
    done
  done

  if (( ${#hits[@]} == 0 )); then
    echo "⚠  $sample not found" >&2
    return
  fi

  # ─── choose the highest version, ignoring the directory path ───
  local newest
  newest=$(
    printf '%s\n' "${hits[@]}" \
    | awk -F/ '{print $NF, $0}' \
    | sort -k1,1V \
    | tail -n1 \
    | cut -d' ' -f2-
  )

  ln -sf "$newest" "$TARGET_DIR/$(basename "$newest")"

  # link the .bai if it exists
  local bai="${newest}.bai"
  [[ -f $bai ]] && ln -sf "$bai" "$TARGET_DIR/$(basename "$bai")"

  echo "linked $(basename "$newest")"
}

export TARGET_DIR irods_roots projects
for s in "${samples[@]}"; do
  link_latest "$s"
done
