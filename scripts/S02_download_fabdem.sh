#!/usr/bin/env bash
# Download FABDEM v1.2 (Hawker et al. 2022) — bare-earth 30 m DEM derived
# from Copernicus GLO-30 with vegetation/buildings removed via ML.
#
# Distribution: University of Bristol data.bris — 10°×10° zip archives by
# SW corner. For the NileDeltaBase bbox we need 4 archives:
#   N20E020-N30E030, N20E030-N30E040, N30E020-N40E030, N30E030-N40E040
#
# Each archive ~1–2 GB. After extraction we keep only the 1°×1° tiles
# that intersect our bbox.

set -euo pipefail

ROOT="${NILEDELTABASE_ROOT:-/data4/EGY/NileDeltaBase}"
DEST="$ROOT/data/raw/fabdem"
STAGING="$DEST/_staging"
LOG="$ROOT/logs/S02_fabdem.log"

mkdir -p "$DEST" "$STAGING"
mkdir -p "$(dirname "$LOG")"

# data.bris base (FABDEM V1.2, DOI 10.5523/bris.s5hqmjcdj8yo2ibzi9b4ew3sn)
BASE_URL="https://data.bris.ac.uk/datasets/s5hqmjcdj8yo2ibzi9b4ew3sn"

echo "=== NileDeltaBase S02: FABDEM v1.2 download ===" | tee "$LOG"
echo "Destination: $DEST" | tee -a "$LOG"
echo "Started: $(date -u +%FT%TZ)" | tee -a "$LOG"

# Archive names (each contains 100 1°×1° tiles)
ARCHIVES=(
  "N20E020-N30E030_FABDEM_V1-2.zip"
  "N20E030-N30E040_FABDEM_V1-2.zip"
  "N30E020-N40E030_FABDEM_V1-2.zip"
  "N30E030-N40E040_FABDEM_V1-2.zip"
)

# Tiles to keep (1°×1°, named by SW corner)
KEEP=()
for LAT in 29 30 31; do
  for LON in 029 030 031 032 033 034 035; do
    KEEP+=("N${LAT}E${LON}_FABDEM_V1-2.tif")
  done
done

for ARC in "${ARCHIVES[@]}"; do
  OUT="$STAGING/$ARC"
  URL="$BASE_URL/$ARC"

  if [[ -f "$OUT" && -s "$OUT" ]]; then
    echo "[skip] $ARC — already downloaded ($(du -h "$OUT" | cut -f1))" | tee -a "$LOG"
  else
    echo "[get ] $URL" | tee -a "$LOG"
    if ! curl -fsS --retry 3 --retry-delay 5 -o "$OUT.part" "$URL"; then
      rm -f "$OUT.part"
      echo "[fail] $ARC — download failed" | tee -a "$LOG"
      continue
    fi
    mv "$OUT.part" "$OUT"
    echo "       -> $(du -h "$OUT" | cut -f1)" | tee -a "$LOG"
  fi

  echo "[unzip] extracting needed tiles from $ARC" | tee -a "$LOG"
  for TILE in "${KEEP[@]}"; do
    unzip -jo "$OUT" "$TILE" -d "$DEST" 2>/dev/null && \
      echo "         got $TILE" | tee -a "$LOG" || true
  done
done

echo "" | tee -a "$LOG"
echo "Tiles in $DEST:" | tee -a "$LOG"
ls "$DEST"/*.tif 2>/dev/null | tee -a "$LOG" || echo "(none)" | tee -a "$LOG"

# Build per-source VRT
if ls "$DEST"/*.tif >/dev/null 2>&1; then
  ls "$DEST"/*.tif > "$ROOT/data/filelists/fabdem.txt"
  gdalbuildvrt -overwrite \
    -input_file_list "$ROOT/data/filelists/fabdem.txt" \
    "$ROOT/data/vrt/fabdem_nile.vrt" 2>&1 | tee -a "$LOG"
  echo "VRT: $ROOT/data/vrt/fabdem_nile.vrt" | tee -a "$LOG"
fi

echo "Finished: $(date -u +%FT%TZ)" | tee -a "$LOG"
