#!/usr/bin/env bash
# Download Copernicus GLO-30 tiles covering the NileDeltaBase bbox
# (29.0–35.1°E, 29.3–31.8°N). Public AWS Open Data bucket, no auth.
#
# Tiles are 1°×1° COGs (~60 MB each). 21 tiles → ~1.3 GB.
# Heights are orthometric (EGM2008).

set -euo pipefail

ROOT="${NILEDELTABASE_ROOT:-/data4/EGY/NileDeltaBase}"
DEST="$ROOT/data/raw/cop_glo30"
LOG="$ROOT/logs/S01_cop_glo30.log"

mkdir -p "$DEST"
mkdir -p "$(dirname "$LOG")"

BASE_URL="https://copernicus-dem-30m.s3.amazonaws.com"

echo "=== NileDeltaBase S01: Copernicus GLO-30 download ===" | tee "$LOG"
echo "Destination: $DEST" | tee -a "$LOG"
echo "Started: $(date -u +%FT%TZ)" | tee -a "$LOG"

N_OK=0; N_FAIL=0
for LAT in N29 N30 N31; do
  for LON in E029 E030 E031 E032 E033 E034 E035; do
    TILE="Copernicus_DSM_COG_10_${LAT}_00_${LON}_00_DEM"
    URL="$BASE_URL/$TILE/$TILE.tif"
    OUT="$DEST/$TILE.tif"

    if [[ -f "$OUT" && -s "$OUT" ]]; then
      echo "[skip] $TILE — already present ($(du -h "$OUT" | cut -f1))" | tee -a "$LOG"
      N_OK=$((N_OK + 1))
      continue
    fi

    echo "[get ] $URL" | tee -a "$LOG"
    if curl -fsS --retry 3 --retry-delay 5 -o "$OUT.part" "$URL"; then
      mv "$OUT.part" "$OUT"
      echo "       -> $(du -h "$OUT" | cut -f1)" | tee -a "$LOG"
      N_OK=$((N_OK + 1))
    else
      rm -f "$OUT.part"
      echo "[fail] $TILE" | tee -a "$LOG"
      N_FAIL=$((N_FAIL + 1))
    fi
  done
done

echo "" | tee -a "$LOG"
echo "Done: $N_OK ok, $N_FAIL failed" | tee -a "$LOG"
echo "Finished: $(date -u +%FT%TZ)" | tee -a "$LOG"

# Build the per-source VRT
ls "$DEST"/*.tif > "$ROOT/data/filelists/cop_glo30.txt"
gdalbuildvrt -overwrite \
  -input_file_list "$ROOT/data/filelists/cop_glo30.txt" \
  "$ROOT/data/vrt/cop_glo30_nile.vrt" 2>&1 | tee -a "$LOG"

echo "VRT: $ROOT/data/vrt/cop_glo30_nile.vrt" | tee -a "$LOG"
