#!/usr/bin/env bash
# Extended-bbox display products (Dabaa → Suez Canal, 28.0–32.5°E, 30.0–32.5°N).
# Produces the same set as S06 but at the wider extent, with "_ext" suffix:
#   nile_elevation_display_ext.tif    bilinear-resampled DEM
#   nile_hillshade_display_ext.tif    gdaldem hillshade (multidirectional)
#   nile_slope_display_ext.tif        gdaldem slope
#   nile_landmask_display_ext.tif     rasterized GSHHG L1
#
# Used by S09 for the extended basemap + zoom panels (Dabaa, Alexandria,
# Cairo, Rosetta, Damietta, Port Said, Ismailia, Manzala, Burullus, etc.).

set -euo pipefail

ROOT="${NILEDELTABASE_ROOT:-/data4/EGY/NileDeltaBase}"
VRT="$ROOT/data/vrt/nile_dem_fused_30m.vrt"
DISP="$ROOT/data/derived/display"
LOG="$ROOT/logs/S06b_derived_ext.log"

mkdir -p "$DISP"
mkdir -p "$(dirname "$LOG")"

echo "=== NileDeltaBase S06b: extended-bbox display products ===" | tee "$LOG"
echo "Started: $(date -u +%FT%TZ)" | tee -a "$LOG"

[[ -f "$VRT" ]] || { echo "Missing $VRT — run S04 first" | tee -a "$LOG"; exit 1; }

# Extended bbox: west to Dabaa (~28.5°E) + buffer, east unchanged.
XMIN=28.0; YMIN=30.0; XMAX=32.5; YMAX=32.5
RES=0.003   # ~330 m display grid

DEM="$DISP/nile_elevation_display_ext.tif"
HLS="$DISP/nile_hillshade_display_ext.tif"
SLP="$DISP/nile_slope_display_ext.tif"
LMSK="$DISP/nile_landmask_display_ext.tif"
GSHHG_L1="$ROOT/data/raw/shoreline/gshhg_GSHHS_f_L1_ext.shp"

echo "[warp] display DEM at ${RES}° (~330 m)" | tee -a "$LOG"
gdalwarp -overwrite \
  -t_srs EPSG:4326 \
  -te $XMIN $YMIN $XMAX $YMAX \
  -tr $RES $RES \
  -r bilinear -wo NUM_THREADS=ALL_CPUS \
  -of GTiff -ot Float32 -dstnodata -9999 \
  -co TILED=YES -co COMPRESS=DEFLATE \
  "$VRT" "$DEM" 2>&1 | tail -3 | tee -a "$LOG"

echo "[hillshade] multidirectional, -s 111120, z=4" | tee -a "$LOG"
gdaldem hillshade -multidirectional -s 111120 -z 4 \
  -co TILED=YES -co COMPRESS=DEFLATE \
  "$DEM" "$HLS" 2>&1 | tail -3 | tee -a "$LOG"

echo "[slope]" | tee -a "$LOG"
gdaldem slope -s 111120 \
  -co TILED=YES -co COMPRESS=DEFLATE \
  "$DEM" "$SLP" 2>&1 | tail -3 | tee -a "$LOG"

echo "[landmask] rasterize GSHHG L1 land polygon at ${RES}°" | tee -a "$LOG"
if [[ -f "$GSHHG_L1" ]]; then
  gdal_rasterize -burn 1 -init 0 -at \
    -te $XMIN $YMIN $XMAX $YMAX \
    -tr $RES $RES \
    -ot Byte -a_nodata 255 \
    -l "$(basename "$GSHHG_L1" .shp)" \
    -co TILED=YES -co COMPRESS=DEFLATE \
    "$GSHHG_L1" "$LMSK" 2>&1 | tail -3 | tee -a "$LOG"
else
  echo "SKIP: $GSHHG_L1 missing — run S05b_download_shoreline.jl (ext bbox)" | tee -a "$LOG"
fi

echo "" | tee -a "$LOG"
echo "Outputs:" | tee -a "$LOG"
ls -lh "$DISP"/*_ext.tif 2>/dev/null | tee -a "$LOG"
echo "Finished: $(date -u +%FT%TZ)" | tee -a "$LOG"
