#!/usr/bin/env bash
set -euo pipefail

WHEEL="$1"
DEST_DIR="$2"

echo "cibw_repair.sh running on $(uname -s) for $WHEEL"

# ---- CRITICAL GUARD ----
# cibuildwheel runs repair-wheel-command on *all* platforms.
# We must skip non-macOS platforms.
if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "Non-macOS platform; skipping custom repair."
  exit 0
fi
# ------------------------

# macOS-only logic below
python -m pip install -U delocate
delocate-wheel -v -w "$DEST_DIR" "$WHEEL"

# (your install_name_tool fixups go here later)

