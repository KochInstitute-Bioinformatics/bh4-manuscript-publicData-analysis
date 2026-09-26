#!/bin/bash
# Restores the data files packed by pack_input_data.sh into their correct
# datasets/<name>/... locations.
#
# Quick start (downloads ~3.6 GB, needs ~8 GB free for tarball + extraction):
#
#   bash restore_input_data.sh
#
# The script will fetch input_data.tar.gz from Google Drive if it is not
# already present, unpack it, and copy input_data/datasets/ onto datasets/.
#
# To download manually instead (either command works):
#
#   curl -L -o input_data.tar.gz \
#     "https://drive.usercontent.google.com/download?id=1x3JLTMNisNe02R4xyrH1we21wqEEObnL&export=download&confirm=t"
#
#   wget -c -O input_data.tar.gz \
#     "https://drive.usercontent.google.com/download?id=1x3JLTMNisNe02R4xyrH1we21wqEEObnL&export=download&confirm=t"
#
#   tar -xzf input_data.tar.gz
#
# Copies input_data/datasets/ straight onto datasets/, merging into the
# existing directory tree -- no per-file list to maintain here, so a
# future input_data.tar.gz with more/different files just works without
# editing this script.
set -euo pipefail
cd "$(dirname "$0")"

FILE_ID=1x3JLTMNisNe02R4xyrH1we21wqEEObnL
TARBALL=input_data.tar.gz
URL="https://drive.usercontent.google.com/download?id=${FILE_ID}&export=download&confirm=t"

# 1. Download the tarball if we have neither it nor the unpacked tree.
if [ ! -d input_data/datasets ] && [ ! -f "$TARBALL" ]; then
  echo "Downloading $TARBALL (~3.6 GB) ..."
  if command -v curl >/dev/null 2>&1; then
    curl -L -C - -o "$TARBALL" "$URL"
  elif command -v wget >/dev/null 2>&1; then
    wget -c -O "$TARBALL" "$URL"
  else
    echo "Need curl or wget to download $TARBALL." >&2
    exit 1
  fi
fi

# 2. Unpack it if needed. Guard against Google handing back an HTML page.
if [ ! -d input_data/datasets ]; then
  if ! gzip -t "$TARBALL" 2>/dev/null; then
    echo "$TARBALL is not a valid gzip file (likely a Google Drive error page)." >&2
    echo "Delete it and retry, or download manually -- see header of this script." >&2
    exit 1
  fi
  echo "Unpacking $TARBALL ..."
  tar -xzf "$TARBALL"
fi

if [ ! -d input_data/datasets ]; then
  echo "input_data/datasets/ not found after unpacking $TARBALL." >&2
  exit 1
fi

# 3. Merge into datasets/.
n_files=$(find input_data/datasets -type f | wc -l)
mkdir -p datasets
cp -R input_data/datasets/. datasets/

echo "Restored $n_files file(s) into datasets/."
