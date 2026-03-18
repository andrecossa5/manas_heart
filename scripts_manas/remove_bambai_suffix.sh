#!/usr/bin/env bash
##############################################################################
# Remove_BAMBAI_version_suffix.sh
#   In $TARGET_DIR, rename
#     PDxxxxx(.loNNNN)?.vN.sample.dupmarked.bam     → PDxxxxx(.loNNNN)?.sample.dupmarked.bam
#     …and the same for *.bam.bai, forcing overwrites
##############################################################################
set -euo pipefail

TARGET_DIR="/lustre/scratch126/cellgen/behjati/md39/foetal_project/run1/Tree"
cd "$TARGET_DIR"

echo "Renaming .bam.bai first…"
for f in *.v*.sample.dupmarked.bam.bai; do
  new=$(printf '%s\n' "$f" | sed -E 's/\.v[0-9]+\.sample\.dupmarked\.bam\.bai$/.sample.dupmarked.bam.bai/')
  [ "$f" = "$new" ] && continue
  echo "  $f → $new"
  mv -f "$f" "$new"
done

echo "Now renaming .bam…"
for f in *.v*.sample.dupmarked.bam; do
  new=$(printf '%s\n' "$f" | sed -E 's/\.v[0-9]+\.sample\.dupmarked\.bam$/.sample.dupmarked.bam/')
  [ "$f" = "$new" ] && continue
  echo "  $f → $new"
  mv -f "$f" "$new"
done

echo "Done."
