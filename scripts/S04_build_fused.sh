#!/usr/bin/env bash
# Build the priority-stacked fused VRT from ingested layers.
# Priority (low → high, last wins in gdalbuildvrt):
#   cop_glo30 → fabdem → deltadtm   [fused_30m]
#   + tandemx12                      [fused_12m, future]
#   + earthdem2                      [fused_2m,  future]

set -euo pipefail

ROOT="${NILEDELTABASE_ROOT:-/data4/EGY/NileDeltaBase}"
VRT="$ROOT/data/vrt"
LOG="$ROOT/logs/S04_fused.log"
mkdir -p "$(dirname "$LOG")"

echo "=== NileDeltaBase S04: build fused VRTs ===" | tee "$LOG"
echo "Started: $(date -u +%FT%TZ)" | tee -a "$LOG"

# Verify source VRTs exist
need_ok=1
for V in cop_glo30_nile.vrt fabdem_nile.vrt deltadtm_nile.vrt; do
  if [[ ! -f "$VRT/$V" ]]; then
    echo "[missing] $V — run S01/S02/S03 first" | tee -a "$LOG"
    need_ok=0
  else
    echo "[ok] $V" | tee -a "$LOG"
  fi
done

if [[ $need_ok -eq 0 ]]; then
  echo "Cannot build fused_30m yet — missing sources" | tee -a "$LOG"
  exit 1
fi

# ── fused_30m: priority GLO30 → FABDEM → DeltaDTM ────────────────────
# Build file list in priority order (later entries win)
FLIST="$ROOT/data/filelists/fused_30m_sources.txt"
{
  echo "$VRT/cop_glo30_nile.vrt"
  echo "$VRT/fabdem_nile.vrt"
  echo "$VRT/deltadtm_nile.vrt"
} > "$FLIST"

gdalbuildvrt -overwrite \
  -resolution highest \
  -r bilinear \
  -input_file_list "$FLIST" \
  "$VRT/nile_dem_fused_30m.vrt" 2>&1 | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "fused_30m: $VRT/nile_dem_fused_30m.vrt" | tee -a "$LOG"
gdalinfo "$VRT/nile_dem_fused_30m.vrt" | grep -E "(Size is|Pixel Size)" | tee -a "$LOG"

# ── fused_12m: stub (includes TanDEM-X 12 m when acquired) ───────────
if [[ -f "$VRT/tandemx12_nile.vrt" ]]; then
  FL12="$ROOT/data/filelists/fused_12m_sources.txt"
  {
    echo "$VRT/cop_glo30_nile.vrt"
    echo "$VRT/fabdem_nile.vrt"
    echo "$VRT/deltadtm_nile.vrt"
    echo "$VRT/tandemx12_nile.vrt"
  } > "$FL12"
  gdalbuildvrt -overwrite -resolution highest -r bilinear \
    -input_file_list "$FL12" \
    "$VRT/nile_dem_fused_12m.vrt" 2>&1 | tee -a "$LOG"
  echo "fused_12m: $VRT/nile_dem_fused_12m.vrt" | tee -a "$LOG"
fi

# ── fused_2m: stub (includes EarthDEM 2 m strips when acquired) ──────
if [[ -f "$VRT/earthdem2_nile.vrt" ]]; then
  FL2="$ROOT/data/filelists/fused_2m_sources.txt"
  {
    echo "$VRT/cop_glo30_nile.vrt"
    echo "$VRT/fabdem_nile.vrt"
    echo "$VRT/deltadtm_nile.vrt"
    [[ -f "$VRT/tandemx12_nile.vrt" ]] && echo "$VRT/tandemx12_nile.vrt"
    echo "$VRT/earthdem2_nile.vrt"
  } > "$FL2"
  gdalbuildvrt -overwrite -resolution highest -r bilinear \
    -input_file_list "$FL2" \
    "$VRT/nile_dem_fused_2m.vrt" 2>&1 | tee -a "$LOG"
  echo "fused_2m: $VRT/nile_dem_fused_2m.vrt" | tee -a "$LOG"
fi

echo "Finished: $(date -u +%FT%TZ)" | tee -a "$LOG"
