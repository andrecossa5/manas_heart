#!/usr/bin/env bash
##############################################################################
# create_symlinks.sh
#   For every sample listed in $SAMPLE_LIST, create symbolic links in
#   $TARGET_DIR pointing at the newest *.sample.dupmarked.bam file and, if it
#   exists, its matching *.bam.bai index.  Works for IDs with or without the
#   “_loNNNN” lane suffix.
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
link_latest () {
  local sample=$1
  local base=${sample%%_lo*}    # PDxxxxx   (without lane suffix)
  local alt=$sample            # PDxxxxx_loNNNN (if suffix present)
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
