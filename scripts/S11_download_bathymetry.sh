#!/usr/bin/env bash
# GMRT (Global Multi-Resolution Topography) bathymetry + topography grid
# for the Nile Delta bbox.  GMRT fuses multibeam swaths with SRTM/GEBCO
# fallbacks — far better near-shore detail than raw GEBCO in the eastern
# Mediterranean (Nile cone, Herodotus basin, Rosetta/Damietta canyons).
#
# Free HTTPS subset API, no auth, GeoTIFF (Int16 metres).
#   https://www.gmrt.org/services/GridServer?...
#
# Resolution: "max" ≈ 100 m offshore where multibeam exists, SRTM
# resampled ~30 m on land.  Bathymetry is our only source <0 m.

set -euo pipefail

ROOT="${NILEDELTABASE_ROOT:-/data4/EGY/NileDeltaBase}"
DEST="$ROOT/data/raw/bathymetry"
LOG="$ROOT/logs/S11_bathymetry.log"
mkdir -p "$DEST" "$(dirname "$LOG")"

echo "=== NileDeltaBase S11: GMRT bathymetry ===" | tee "$LOG"
echo "Started: $(date -u +%FT%TZ)" | tee -a "$LOG"

# Extended bbox — includes Mediterranean shelf north of the Delta and
# Herodotus basin offshore of Rosetta/Damietta.  Pad 0.5° north for
# offshore context.
XMIN=28.0; YMIN=30.0; XMAX=35.5; YMAX=33.0
OUT="$DEST/gmrt_nile.tif"

if [[ -f "$OUT" && -s "$OUT" ]]; then
  echo "[skip] $OUT already present ($(du -h "$OUT" | cut -f1))" | tee -a "$LOG"
else
  URL="https://www.gmrt.org/services/GridServer?west=${XMIN}&east=${XMAX}&south=${YMIN}&north=${YMAX}&format=geotiff&resolution=max&layer=topo"
  echo "[get ] $URL" | tee -a "$LOG"
  curl -fsS --retry 3 --retry-delay 5 -o "$OUT.part" "$URL" && mv "$OUT.part" "$OUT"
  echo "       -> $(du -h "$OUT" | cut -f1)" | tee -a "$LOG"
fi

echo "" | tee -a "$LOG"
echo "[gdalinfo] summary" | tee -a "$LOG"
gdalinfo -stats "$OUT" 2>&1 | grep -E "Size|Pixel Size|Minimum|Maximum|Mean|NoData" | tee -a "$LOG"

# Build VRT + ocean-only clip (bathy < 0) at display grids.
VRTDIR="$ROOT/data/vrt"
DISPDIR="$ROOT/data/derived/display"
mkdir -p "$VRTDIR" "$DISPDIR"

gdalbuildvrt -overwrite "$VRTDIR/gmrt_nile.vrt" "$OUT" 2>&1 | tail -3 | tee -a "$LOG"

# Ocean-only bathymetry (<0 m) at the extended display grid (0.003°).
DEPTH_EXT="$DISPDIR/nile_bathymetry_display_ext.tif"
echo "[warp] bathy clip — ext bbox" | tee -a "$LOG"
gdalwarp -overwrite \
  -t_srs EPSG:4326 \
  -te 28.0 30.0 32.5 32.5 \
  -tr 0.003 0.003 \
  -r bilinear -wo NUM_THREADS=ALL_CPUS \
  -of GTiff -ot Float32 -dstnodata -9999 \
  -co TILED=YES -co COMPRESS=DEFLATE \
  "$VRTDIR/gmrt_nile.vrt" "$DEPTH_EXT" 2>&1 | tail -3 | tee -a "$LOG"

echo "Finished: $(date -u +%FT%TZ)" | tee -a "$LOG"
