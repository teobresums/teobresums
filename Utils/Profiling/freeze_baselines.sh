#!/usr/bin/env bash
# Runs the given parfiles with the opt build and stores their output and an md5 manifest, to check
# later that a change leaves the output the same.
#
# Usage:
#   Utils/Profiling/freeze_baselines.sh DEST_DIR PARFILE [PARFILE ...]
set -euo pipefail

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
CDIR="$REPO/C"
DEST="$1"; shift

echo "[baselines] building opt ..."
make -C "$CDIR" opt >/dev/null

mkdir -p "$DEST"
MANIFEST="$DEST/manifest.md5"
: > "$MANIFEST"

for p in "$@"; do
  p="$(cd "$(dirname "$p")" && pwd)/$(basename "$p")"
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
