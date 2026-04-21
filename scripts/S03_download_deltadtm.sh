#!/usr/bin/env bash
# Download DeltaDTM v1.1 (Pronk et al. 2024, Scientific Data) — bare-earth
# 30 m coastal DEM corrected with ICESat-2 ATL08 + GEDI L2A and non-terrain
# morphological filtering. Best free product for flat coastal plains (<30 m).
#
# Distribution: 4TU.ResearchData DOI:10.4121/21997565.v4
# Cloud-Optimized GeoTIFF, zipped by continent (Africa for Nile Delta).

set -euo pipefail

ROOT="${NILEDELTABASE_ROOT:-/data4/EGY/NileDeltaBase}"
DEST="$ROOT/data/raw/deltadtm"
STAGING="$DEST/_staging"
LOG="$ROOT/logs/S03_deltadtm.log"

mkdir -p "$DEST" "$STAGING"
mkdir -p "$(dirname "$LOG")"

# 4TU ResearchData file listing — the dataset hosts per-continent archives.
# Africa file: DeltaDTM_v1.1_Africa.zip (confirmed via 4TU v4 file listing).
# Resolve the real download URL via the 4TU API at runtime.
DOI="10.4121/21997565.v4"

echo "=== NileDeltaBase S03: DeltaDTM v1.1 download ===" | tee "$LOG"
echo "Destination: $DEST" | tee -a "$LOG"
echo "Started: $(date -u +%FT%TZ)" | tee -a "$LOG"

# Resolved from 4TU API v2 file listing (DOI 10.4121/21997565.v4).
# Article UUID:  1da2e70f-6c4d-4b03-86bd-b53e789cc629
# Africa.zip UUID: 22ffa027-184b-4f67-9979-c182f3dfb1ab  (2.77 GB)
# tiles.gpkg UUID: 60a69899-2e67-4f9f-8761-3b57094acd12  (2.9 MB, tile index)
ARTICLE_UUID="1da2e70f-6c4d-4b03-86bd-b53e789cc629"
AFRICA_UUID="22ffa027-184b-4f67-9979-c182f3dfb1ab"
TILES_UUID="60a69899-2e67-4f9f-8761-3b57094acd12"

AFRICA_URL="https://data.4tu.nl/file/$ARTICLE_UUID/$AFRICA_UUID"
TILES_URL="https://data.4tu.nl/file/$ARTICLE_UUID/$TILES_UUID"

# Fetch the tile index (small, useful for knowing which tiles exist)
if [[ ! -f "$DEST/deltadtm_tiles.gpkg" ]]; then
  echo "[get ] tile index: $TILES_URL" | tee -a "$LOG"
  curl -fsSL --retry 3 --retry-delay 5 -o "$DEST/deltadtm_tiles.gpkg" "$TILES_URL" \
    && echo "       -> $DEST/deltadtm_tiles.gpkg ($(du -h "$DEST/deltadtm_tiles.gpkg" | cut -f1))" | tee -a "$LOG"
fi

ARC_NAME="Africa.zip"
OUT="$STAGING/$ARC_NAME"

if [[ -f "$OUT" && -s "$OUT" ]]; then
  echo "[skip] $ARC_NAME — already downloaded ($(du -h "$OUT" | cut -f1))" | tee -a "$LOG"
else
  echo "[get ] $AFRICA_URL" | tee -a "$LOG"
  curl -fsSL --retry 3 --retry-delay 5 -o "$OUT.part" "$AFRICA_URL"
  mv "$OUT.part" "$OUT"
  echo "       -> $(du -h "$OUT" | cut -f1)" | tee -a "$LOG"
fi

# Extract contents
echo "[unzip] extracting all contents of $ARC_NAME" | tee -a "$LOG"
unzip -o "$OUT" -d "$DEST" | tail -20 | tee -a "$LOG"

# Build VRT from all COGs in DEST (any subdirectory)
echo "[vrt ] building deltadtm_nile.vrt from extracted tiles" | tee -a "$LOG"
find "$DEST" -name '*.tif' ! -path "*/_staging/*" > "$ROOT/data/filelists/deltadtm_all.txt"

# Filter to tiles that intersect the bbox (gdaltindex-based or filename-based)
# DeltaDTM tiles use a different naming convention (1° grid, bare-earth only
# for coastal areas). We'll include all matching tiles; gdalbuildvrt handles
# non-overlap gracefully.
gdalbuildvrt -overwrite \
  -input_file_list "$ROOT/data/filelists/deltadtm_all.txt" \
  "$ROOT/data/vrt/deltadtm_global.vrt" 2>&1 | tee -a "$LOG"

# Clip to our bbox for efficiency
gdal_translate \
  -projwin 29.0 31.8 35.1 29.3 \
  -co TILED=YES -co COMPRESS=DEFLATE \
  "$ROOT/data/vrt/deltadtm_global.vrt" \
  "$DEST/deltadtm_nile_clipped.tif" 2>&1 | tee -a "$LOG"

echo "$DEST/deltadtm_nile_clipped.tif" > "$ROOT/data/filelists/deltadtm.txt"
gdalbuildvrt -overwrite \
  -input_file_list "$ROOT/data/filelists/deltadtm.txt" \
  "$ROOT/data/vrt/deltadtm_nile.vrt" 2>&1 | tee -a "$LOG"

echo "Finished: $(date -u +%FT%TZ)" | tee -a "$LOG"
