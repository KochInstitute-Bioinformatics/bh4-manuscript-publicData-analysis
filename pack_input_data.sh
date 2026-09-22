#!/bin/bash
# Collects every gitignored data file under datasets/ (prep-produced/
# ready-to-score .rds checkpoints, plus each dataset's raw data/ or
# from_geo/ directory) into input_data/, mirroring the exact
# datasets/<name>/... path each file already has -- then tars it up as
# input_data.tar.gz.
#
# Run this whenever a dataset's data changes, to refresh the tarball that
# gets uploaded somewhere durable (Zenodo/Dropbox/lab share -- GitHub
# itself only ever gets this script, never the tarball). A fresh clone
# gets the data back with restore_input_data.sh -- see that script and
# README's "Large files / data availability" section.
set -euo pipefail
cd "$(dirname "$0")"

DEST=input_data
rm -rf "$DEST"
mkdir -p "$DEST"

echo "Collecting .rds checkpoints..."
find datasets -type f -name '*.rds' -print0 | while IFS= read -r -d '' f; do
  mkdir -p "$DEST/$(dirname "$f")"
  cp "$f" "$DEST/$f"
done

echo "Collecting raw data/ and from_geo/ directories..."
for sub in data from_geo; do
  find datasets -mindepth 2 -maxdepth 2 -type d -name "$sub" -print0 | while IFS= read -r -d '' d; do
    mkdir -p "$DEST/$(dirname "$d")"
    cp -r "$d" "$DEST/$d"
  done
done

echo "Creating input_data.tar.gz..."
tar -czf input_data.tar.gz "$DEST"

echo "Done: $(du -sh "$DEST" | cut -f1) in $DEST/, tarball $(du -sh input_data.tar.gz | cut -f1)"
