#!/usr/bin/env bash
# ESA WorldCover v200 (2021) — global 10 m land cover from Sentinel-1/2.
# 11 classes:
#   10 tree cover           20 shrubland        30 grassland
#   40 cropland             50 built-up         60 bare/sparse
#   70 snow/ice             80 perm. water      90 herbaceous wetland
#   95 mangroves            100 moss/lichen
#
# Tiles are 3°×3° GeoTIFFs (Byte, 10 m / ~1 arcsec).
# AWS Open Data:  s3://esa-worldcover/v200/2021/map/<tile>.tif
# Naming:         ESA_WorldCover_10m_2021_v200_<N..E..>.tif where N.. is SW lat
#                 (multiples of 3, e.g. N30), E.. is SW lon (multiples of 3, e.g. E027).
#
# For bbox 28.0–32.5°E, 30.0–32.5°N we need one 3×3 tile: N30E027 (covers 27–30°E)
# plus N30E030 (covers 30–33°E).

set -euo pipefail

ROOT="${NILEDELTABASE_ROOT:-/data4/EGY/NileDeltaBase}"
DEST="$ROOT/data/raw/worldcover"
LOG="$ROOT/logs/S05c_worldcover.log"
mkdir -p "$DEST" "$(dirname "$LOG")"

BASE_URL="https://esa-worldcover.s3.eu-central-1.amazonaws.com/v200/2021/map"

echo "=== NileDeltaBase S05c: ESA WorldCover v200 (2021) ===" | tee "$LOG"
echo "Started: $(date -u +%FT%TZ)" | tee -a "$LOG"

TILES=(
  "ESA_WorldCover_10m_2021_v200_N30E027_Map.tif"
  "ESA_WorldCover_10m_2021_v200_N30E030_Map.tif"
)

for TILE in "${TILES[@]}"; do
  OUT="$DEST/$TILE"
  URL="$BASE_URL/$TILE"
  if [[ -f "$OUT" && -s "$OUT" ]]; then
    echo "[skip] $TILE — already present ($(du -h "$OUT" | cut -f1))" | tee -a "$LOG"
    continue
  fi
  echo "[get ] $URL" | tee -a "$LOG"
  curl -fsS --retry 3 --retry-delay 5 -o "$OUT.part" "$URL" && mv "$OUT.part" "$OUT"
  echo "       -> $(du -h "$OUT" | cut -f1)" | tee -a "$LOG"
done

# VRT over the downloaded tiles
ls "$DEST"/*.tif > "$ROOT/data/filelists/worldcover.txt"
gdalbuildvrt -overwrite \
  -input_file_list "$ROOT/data/filelists/worldcover.txt" \
  "$ROOT/data/vrt/worldcover_nile.vrt" 2>&1 | tail -3 | tee -a "$LOG"

# Clip to the extended display bbox at the display resolution.  Use nearest-
# neighbour so class codes are preserved (Byte).
XMIN=28.0; YMIN=30.0; XMAX=32.5; YMAX=32.5
RES=0.003
CLIP="$ROOT/data/derived/display/nile_worldcover_display_ext.tif"
echo "[warp] WorldCover clip at ${RES}° (nearest)" | tee -a "$LOG"
gdalwarp -overwrite \
  -t_srs EPSG:4326 \
  -te $XMIN $YMIN $XMAX $YMAX \
  -tr $RES $RES \
  -r near -wo NUM_THREADS=ALL_CPUS \
  -of GTiff -ot Byte -dstnodata 0 \
  -co TILED=YES -co COMPRESS=DEFLATE \
  "$ROOT/data/vrt/worldcover_nile.vrt" "$CLIP" 2>&1 | tail -3 | tee -a "$LOG"

echo "Clipped: $CLIP" | tee -a "$LOG"
ls -lh "$CLIP" | tee -a "$LOG"
echo "Finished: $(date -u +%FT%TZ)" | tee -a "$LOG"
