#!/usr/bin/env bash
##############################################################################
# create_links_blood.sh
#   For every sample listed in the Sample_ID column of $CSV, create symbolic
#   links in $TARGET_FOLDER pointing at the newest *.sample.dupmarked.bam
#   file and, if it exists, its matching *.bam.bai index.
#
##############################################################################
set -euo pipefail

# ──────────────────────────  paths ──────────────────────────
CSV="$(dirname "$0")/../data/input/blood_samples.csv"
TARGET_FOLDER="/lustre/scratch125/cellgen/behjati/ac87/manas_heart/work/bams/blood"

mkdir -p "$TARGET_FOLDER"

# ─────────────────────────  projects ────────────────────────
projects=(3187)

# ───────────  irods-CGP roots visible under /nfs  ───────────
mapfile -t irods_roots < <(find /nfs -maxdepth 1 -type d -name 'irods-cgp*' -print)

# ───────────  load sample IDs from the CSV Sample column  ───
if [[ ! -r $CSV ]]; then
  echo "ERROR: cannot read $CSV" >&2
  exit 1
fi
mapfile -t samples < <(
  awk -F',' 'NR>1 && $1!="" {print $1}' "$CSV" | sort -u
)
printf 'Loaded %d samples\n' "${#samples[@]}"

# ─────────────────────────  helper function  ────────────────
link_latest () {
  local sample=$1
  local base=${sample%%_lo*}    # PDxxxxx   (without lane suffix)
  local alt=$sample             # PDxxxxx_loNNNN (if suffix present)
  local hits=()

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

  ln -sf "$newest" "$TARGET_FOLDER/$(basename "$newest")"

  local bai="${newest}.bai"
  [[ -f $bai ]] && ln -sf "$bai" "$TARGET_FOLDER/$(basename "$bai")"

  echo "linked $(basename "$newest")"
}

export TARGET_FOLDER irods_roots projects
for s in "${samples[@]}"; do
  link_latest "$s"
done
