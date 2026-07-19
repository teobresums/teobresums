#!/usr/bin/env bash
# Freeze reference waveform/dynamics outputs for all profiling configs using the
# optimized (-O2) build, plus an md5 manifest. Bit-identical optimizations (e.g.
# the redundant-recomputation cache, geometric buffer growth) are verified by
# regenerating outputs and diffing against these references / the manifest.
#
# Usage:
#   Utils/Profiling/freeze_baselines.sh [DEST_DIR]
# DEST_DIR defaults to Utils/Profiling/baselines/<git-rev>/ .
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CDIR="$REPO/C"
PARDIR="$REPO/Utils/Profiling/parfiles"
REV="$(git -C "$REPO" rev-parse --short HEAD)"
DEST="${1:-$REPO/Utils/Profiling/baselines/$REV}"

echo "[baselines] building opt ..."
make -C "$CDIR" opt >/dev/null

mkdir -p "$DEST"
MANIFEST="$DEST/manifest.md5"
: > "$MANIFEST"

for p in "$PARDIR"/*.par; do
  name="$(basename "$p" .par)"
  odir="$(grep '^output_dir' "$p" | sed 's/.*"\.\/\(.*\)\/".*/\1/')"
  echo "[baselines] $name"
  ( cd "$CDIR" && rm -rf "$odir" && ./TEOBResumS.x -p "$p" >/dev/null 2>&1 )
  # Copy data files (exclude params.txt: it records the run's output_dir path).
  mkdir -p "$DEST/$name"
  for f in "$CDIR/$odir"/*; do
    bn="$(basename "$f")"
    [ "$bn" = "params.txt" ] && continue
    cp "$f" "$DEST/$name/$bn"
    ( cd "$DEST" && md5sum "$name/$bn" >> "$MANIFEST" )
  done
  ( cd "$CDIR" && rm -rf "$odir" )
done

echo "[baselines] wrote $DEST ($(wc -l < "$MANIFEST") files hashed)"
