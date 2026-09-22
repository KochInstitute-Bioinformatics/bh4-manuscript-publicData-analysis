#!/bin/bash
# Restores the data files packed by pack_input_data.sh into their correct
# datasets/<name>/... locations. Run this after cloning the repo and
# unpacking input_data.tar.gz (obtained separately -- see README's "Large
# files / data availability" section for where to get it):
#
#   tar -xzf input_data.tar.gz
#   ./restore_input_data.sh
#
# Copies input_data/datasets/ straight onto datasets/, merging into the
# existing directory tree -- no per-file list to maintain here, so a
# future input_data.tar.gz with more/different files just works without
# editing this script.
set -euo pipefail
cd "$(dirname "$0")"

if [ ! -d input_data/datasets ]; then
  echo "input_data/datasets/ not found." >&2
  echo "Download input_data.tar.gz first, then: tar -xzf input_data.tar.gz" >&2
  exit 1
fi

n_files=$(find input_data/datasets -type f | wc -l)
cp -R input_data/datasets/. datasets/

echo "Restored $n_files file(s) into datasets/."
