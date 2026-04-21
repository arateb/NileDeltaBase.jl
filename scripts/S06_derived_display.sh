#!/usr/bin/env bash
# Build display-resolution (0.005° ≈ 500 m) products from fused_30m:
#   elevation.tif       bilinear-resampled DEM
#   hillshade.tif       gdaldem hillshade (multidirectional)
#   slope.tif           gdaldem slope
#   below_zero.tif      mask of elevation < 0 (uint8, 0/1)
# These are the inputs for the Julia plotting script.

set -euo pipefail

ROOT="${NILEDELTABASE_ROOT:-/data4/EGY/NileDeltaBase}"
VRT="$ROOT/data/vrt/nile_dem_fused_30m.vrt"
DISP="$ROOT/data/derived/display"
LOG="$ROOT/logs/S06_derived.log"

mkdir -p "$DISP"
mkdir -p "$(dirname "$LOG")"

echo "=== NileDeltaBase S06: derived display products ===" | tee "$LOG"
echo "Started: $(date -u +%FT%TZ)" | tee -a "$LOG"

[[ -f "$VRT" ]] || { echo "Missing $VRT — run S04 first" | tee -a "$LOG"; exit 1; }

XMIN=29.0; YMIN=29.3; XMAX=35.1; YMAX=31.8
RES=0.003   # ~330 m display grid → 2033 × 834 pixels

DEM="$DISP/nile_elevation_display.tif"
HLS="$DISP/nile_hillshade_display.tif"
SLP="$DISP/nile_slope_display.tif"
BZR="$DISP/nile_below_zero_display.tif"

echo "[warp] display DEM at ${RES}° (~330 m)" | tee -a "$LOG"
gdalwarp -overwrite \
  -t_srs EPSG:4326 \
  -te $XMIN $YMIN $XMAX $YMAX \
  -tr $RES $RES \
  -r bilinear -wo NUM_THREADS=ALL_CPUS \
  -of GTiff -ot Float32 -dstnodata -9999 \
  -co TILED=YES -co COMPRESS=DEFLATE \
  "$VRT" "$DEM" 2>&1 | tail -5 | tee -a "$LOG"

echo "[hillshade] multidirectional" | tee -a "$LOG"
gdaldem hillshade -multidirectional \
  -co TILED=YES -co COMPRESS=DEFLATE \
  "$DEM" "$HLS" 2>&1 | tail -3 | tee -a "$LOG"

echo "[slope]" | tee -a "$LOG"
gdaldem slope -s 111120 \
  -co TILED=YES -co COMPRESS=DEFLATE \
  "$DEM" "$SLP" 2>&1 | tail -3 | tee -a "$LOG"

echo "[below-zero mask]" | tee -a "$LOG"
gdal_calc.py \
  -A "$DEM" \
  --outfile="$BZR" \
  --calc='where((A < 0) & (A > -100), 1, 0)' \
  --type=Byte --NoDataValue=0 \
  --co='TILED=YES' --co='COMPRESS=DEFLATE' \
  --overwrite 2>&1 | tail -3 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "Outputs:" | tee -a "$LOG"
ls -lh "$DISP"/ | tee -a "$LOG"

echo "Finished: $(date -u +%FT%TZ)" | tee -a "$LOG"
