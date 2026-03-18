#!/usr/bin/env bash
##############################################################################
# create_symlinks_bai.sh
#   Link the newest *.sample.dupmarked.bam.bai for every sample listed
#   in all_samples.txt.  Works for IDs with or without "_loNNNN" suffix.
##############################################################################
set -euo pipefail

# ──────────────────────────  paths ──────────────────────────
TARGET_DIR="/lustre/scratch126/cellgen/behjati/md39/foetal_project/run1/Tree"
SAMPLE_LIST="$TARGET_DIR/all_samples.txt"

mkdir -p "$TARGET_DIR"

# ─────────────────────────  projects ────────────────────────
projects=(2855 3254 3186 3187 3152 3137)

# ───────────  irods-CGP roots visible under /nfs  ───────────
mapfile -t irods_roots < <(find /nfs -maxdepth 1 -type d -name 'irods-cgp*' -print)

# ───────────  load unique sample IDs from the text file  ────
if [[ ! -r $SAMPLE_LIST ]]; then
  echo "ERROR: cannot read $SAMPLE_LIST" >&2
  exit 1
fi
mapfile -t samples < <(
  grep -Eo 'PD[0-9]{5,6}[A-Za-z]{0,2}(_lo[0-9]+)?' "$SAMPLE_LIST" | sort -u
)
printf 'Loaded %d samples\n' "${#samples[@]}"

# ─────────────────────────  helper function  ────────────────
link_latest_bai () {
  local sample=$1
  local base=${sample%%_lo*}   # PDxxxxx
  local alt=$sample            # PDxxxxx_loNNNN (if present)
  local hits=()

  for root in "${irods_roots[@]}"; do
    for proj in "${projects[@]}"; do
      for d in "$root/intproj/$proj/sample/$base" "$root/intproj/$proj/sample/$alt"; do
        [[ -d $d ]] || continue
        shopt -s globstar nullglob
        # look for both versioned and un-versioned .bam.bai
        hits+=( "$d"/**/"$sample".v*.sample.dupmarked.bam.bai )
        hits+=( "$d"/**/"$sample".sample.dupmarked.bam.bai )
        shopt -u globstar
      done
    done
  done

  if (( ${#hits[@]} == 0 )); then
    echo "⚠  $sample  .bai not found" >&2
    return
  fi

  # ─── pick the highest version number ──────────────────────
  # 1. prepend basename (e.g. PD53943lt.v2.sample.dupmarked.bam.bai)
  # 2. sort -V only on that field
  local newest_bai
  newest_bai=$(
    printf '%s\n' "${hits[@]}" \
    | awk -F/ '{print $NF, $0}' \
    | sort -k1,1V \
    | tail -n1 \
    | cut -d' ' -f2-
  )

  ln -sf "$newest_bai" "$TARGET_DIR/$(basename "$newest_bai")"
  echo "linked $(basename "$newest_bai")"
}

export TARGET_DIR irods_roots projects
for s in "${samples[@]}"; do
  link_latest_bai "$s"
done
