#!/usr/bin/env julia
# Publication figure of the fused Nile Delta DEM.
# Single RGB image built as hypsometric tint × hillshade over *every* valid
# land pixel (including below-sea-level polders) — the DEM itself defines the
# coastline.  Only true water bodies carry non-DEM colour:
#   - ocean (NaN pixels connected to the bbox boundary) → OCEAN_BLUE
#   - enclosed NaN (lakes / lagoons)                    → LAKE_BLUE
# Everything else is hypsometric(elev) × shade(hillshade, slope).

using Pkg
Pkg.activate(dirname(@__DIR__))

using CairoMakie
using GeoMakie
using NaturalEarth

include("/data/files/pkgs/figstd/src/FigStd.jl")
using .FigStd
FigStd.apply!()

push!(LOAD_PATH, joinpath(dirname(@__DIR__), "src"))
using NileDeltaBase
const NDB = NileDeltaBase
NDB.set_root!(get(ENV, "NILEDELTABASE_ROOT", "/data4/EGY/NileDeltaBase"))

using Statistics

const ROOT = NDB._root()
const DEM_TIF = joinpath(ROOT, "data/derived/display/nile_elevation_display.tif")
const HLS_TIF = joinpath(ROOT, "data/derived/display/nile_hillshade_display.tif")
const SLP_TIF = joinpath(ROOT, "data/derived/display/nile_slope_display.tif")

isfile(DEM_TIF) || error("Missing $DEM_TIF — run S06_derived_display.sh first")
isfile(HLS_TIF) || error("Missing $HLS_TIF — run S06_derived_display.sh first")
isfile(SLP_TIF) || error("Missing $SLP_TIF — run S06_derived_display.sh first")

# ── Raster loading via GDAL CLI (ENVI binary) ───────────────────────────────

function _read_raster(tif::String)
    info = read(`gdalinfo -json $tif`, String)
    m_size = match(r"\"size\":\s*\[\s*(\d+),\s*(\d+)\s*\]", info)
    nx = parse(Int, m_size[1]); ny = parse(Int, m_size[2])
    m_geo = match(r"\"geoTransform\":\s*\[([^\]]+)\]", info)
    gt = parse.(Float64, strip.(split(m_geo[1], ',')))
    xmin = gt[1]; dx = gt[2]
    ymax = gt[4]; dy = gt[6]
    m_nd = match(r"\"noDataValue\":\s*([-\d.eE+]+)", info)
    nd = m_nd !== nothing ? parse(Float32, m_nd[1]) : -9999f0

    bin = tempname() * ".bin"
    run(pipeline(`gdal_translate -of ENVI -ot Float32 $tif $bin`; stderr=devnull))
    data = Array{Float32}(undef, nx, ny)
    read!(bin, data)
    rm(bin; force=true)
    rm(bin * ".aux.xml"; force=true)
    rm(replace(bin, ".bin" => ".hdr"); force=true)

    arr = Matrix{Float32}(undef, ny, nx)
    for r in 1:ny
        arr[r, :] = data[:, ny - r + 1]
    end
    arr[arr .== nd] .= NaN32
    ymin = ymax + dy*ny
    lons = collect(range(xmin + dx/2;    step=dx,  length=nx))
    lats = collect(range(ymin + (-dy)/2; step=-dy, length=ny))
    return arr, lons, lats
end

@info "Loading elevation, hillshade, slope..."
elev, lons, lats = _read_raster(DEM_TIF)
hlsh, _,    _    = _read_raster(HLS_TIF)
slp,  _,    _    = _read_raster(SLP_TIF)
ny, nx = size(elev)

valid    = .!isnan.(elev)
invalid  = .!valid

@info "Grid" size=(ny, nx) lon=(first(lons), last(lons)) lat=(first(lats), last(lats))
@info "Elevation stats" min=minimum(elev[valid]) max=maximum(elev[valid]) median=median(elev[valid])
@info "Hillshade stats" min=minimum(hlsh[valid]) max=maximum(hlsh[valid]) median=median(hlsh[valid])

# Use the raw hillshade — a percentile stretch would compress the desert's
# true variation (the Delta is flat, so p2/p98 covers only a narrow band
# and real escarpments saturate).  With -z 4 the native [0, 255] range
# already separates flat Delta (~180–210) from the jebels (< 120).
hlsh_n = clamp.(hlsh ./ 255f0, 0f0, 1f0)

# Slope-darkening multiplier — steep slopes (≳ 15°) get visibly darker.
# Slopes are in degrees; clamp to [0, 35] and map to [0, 1].
slp_n  = clamp.(slp ./ 35f0, 0f0, 1f0)
slp_n[invalid] .= 0f0

# ── Classify the Mediterranean only ─────────────────────────────────────────
# Flood-fill from the north edge with a strict |elev| ≤ tol criterion so
# the Delta lagoons' negative bathymetry doesn't leak through.  The
# lagoons themselves appear naturally in the hypsometric tint.

function _flood_sea_strict!(sea::AbstractMatrix{Bool},
                            elev::AbstractMatrix{Float32};
                            tol::Float32=0.01f0)
    ny, nx = size(elev)
    q = Tuple{Int,Int}[]
    for j in 1:nx
        e = elev[ny, j]
        if !isnan(e) && abs(e) <= tol
            sea[ny, j] = true
            push!(q, (ny, j))
        end
    end
    head = 1
    while head <= length(q)
        i, j = q[head]; head += 1
        for (di, dj) in ((-1,0),(1,0),(0,-1),(0,1))
            ni, nj = i+di, j+dj
            if 1 <= ni <= ny && 1 <= nj <= nx && !sea[ni, nj]
                e = elev[ni, nj]
                if !isnan(e) && abs(e) <= tol
                    sea[ni, nj] = true
                    push!(q, (ni, nj))
                end
            end
        end
    end
end

ocean = falses(ny, nx)
_flood_sea_strict!(ocean, elev)
ocean .|= invalid
lake = falses(ny, nx)
@info "Water classes" ocean=sum(ocean) lake=sum(lake)

# ── Colour palette ──────────────────────────────────────────────────────────

const OCEAN_BLUE = RGBf(0.55, 0.73, 0.88)
const LAKE_BLUE  = RGBf(0.12, 0.28, 0.52)   # lakes / lagoons (darker, distinct)

# Hypsometric tint — densely sampled across the Delta's *actual* elevation
# range (≈ −25 m to +170 m).  Wider stops make the Delta look like one flat
# green blob, so we squeeze most of the palette into the 0–50 m band where
# almost all the variability lives.  Every stop is a clearly distinct hue
# so the gradient is visible even in the flat plain.
const HYPSO_STOPS = [
    (-25f0,  RGBf(0.30, 0.20, 0.50)),   # deep indigo (subsidence hotspots)
    (-10f0,  RGBf(0.45, 0.40, 0.70)),   # violet
    (-2f0,   RGBf(0.60, 0.68, 0.85)),   # pale periwinkle
    (0f0,    RGBf(0.78, 0.90, 0.80)),   # coastal green
    (3f0,    RGBf(0.60, 0.82, 0.55)),   # bright green (Delta plain)
    (8f0,    RGBf(0.82, 0.90, 0.55)),   # yellow-green
    (15f0,   RGBf(0.96, 0.90, 0.52)),   # pale yellow
    (25f0,   RGBf(0.94, 0.78, 0.45)),   # tan
    (50f0,   RGBf(0.84, 0.62, 0.38)),   # brown
    (100f0,  RGBf(0.68, 0.46, 0.28)),   # dark brown
    (170f0,  RGBf(0.46, 0.32, 0.22)),   # deep brown (desert peaks)
]

@inline function hypsometric(e::Float32)
    e <= HYPSO_STOPS[1][1] && return HYPSO_STOPS[1][2]
    @inbounds for k in 1:length(HYPSO_STOPS)-1
        e0, c0 = HYPSO_STOPS[k]
        e1, c1 = HYPSO_STOPS[k+1]
        if e <= e1
            t = (e - e0) / (e1 - e0)
            return RGBf((1-t)*c0.r + t*c1.r,
                        (1-t)*c0.g + t*c1.g,
                        (1-t)*c0.b + t*c1.b)
        end
    end
    return HYPSO_STOPS[end][2]
end

# Hillshade × slope darkening.  Multiplier in ~[0.1, 1.15]; the low end is
# intentionally dark so desert jebels and escarpments read as real shadows.
@inline function shade(c::RGBf, hn::Float32, sn::Float32)
    hn = isnan(hn) ? 0.5f0 : hn
    s  = 0.10f0 + 1.05f0 * hn - 0.35f0 * sn
    s  = clamp(s, 0.05f0, 1.15f0)
    RGBf(clamp(c.r * s, 0f0, 1f0),
         clamp(c.g * s, 0f0, 1f0),
         clamp(c.b * s, 0f0, 1f0))
end

# ── Build the single RGB image ──────────────────────────────────────────────

rgb = Matrix{RGBf}(undef, ny, nx)
@inbounds for i in 1:ny, j in 1:nx
    if ocean[i, j]
        rgb[i, j] = OCEAN_BLUE
    elseif lake[i, j]
        rgb[i, j] = LAKE_BLUE
    else
        rgb[i, j] = shade(hypsometric(elev[i, j]), hlsh_n[i, j], slp_n[i, j])
    end
end

# ── Build the figure ────────────────────────────────────────────────────────

bbox = NDB.BBOX
map_wh = (bbox.lon_max - bbox.lon_min) / (bbox.lat_max - bbox.lat_min)
# ~22% extra width for the legend column; reflects in figure height.
fig = FigStd.figure(width=:wide, aspect=1.0 / map_wh * 1.22)
ax = Axis(fig[1, 1];
          xlabel="Longitude (°E)", ylabel="Latitude (°N)",
          title="Nile Delta — fused 30 m DEM (COP-GLO30 / FABDEM / DeltaDTM)",
          aspect=DataAspect())

# image! takes (x, y, mat) with mat[i,j] at (x[i], y[j]).  Our rgb is (ny, nx)
# with rows = lats ascending, cols = lons ascending; pass (lons, lats, rgb').
image!(ax, (first(lons), last(lons)), (first(lats), last(lats)),
       permutedims(rgb); interpolate=false)

# DEM-derived shorelines — contour the (ocean ∪ lake) mask.
water_f = Matrix{Float32}(undef, ny, nx)
@inbounds for i in 1:ny, j in 1:nx
    water_f[i, j] = (ocean[i, j] || lake[i, j]) ? 1f0 : 0f0
end
contour!(ax, lons, lats, permutedims(water_f);
         levels=[0.5f0], color=RGBf(0.10, 0.10, 0.10),
         linewidth=0.7)

# Country borders (dashed).  Convert FeatureCollection → basic geometries so
# plotting works on a plain Axis.
function _plot_lines!(ax, fc; color, linewidth=0.5, linestyle=:solid)
    for feature in fc
        geom = GeoMakie.geo2basic(feature.geometry)
        if geom isa AbstractVector
            for g in geom
                lines!(ax, g; color=color, linewidth=linewidth, linestyle=linestyle)
            end
        else
            lines!(ax, geom; color=color, linewidth=linewidth, linestyle=linestyle)
        end
    end
end

try
    _plot_lines!(ax, naturalearth("admin_0_boundary_lines_land", 10);
                 color=RGBAf(0.10, 0.10, 0.10, 0.85),
                 linewidth=0.6, linestyle=:dash)
catch e
    @warn "borders 10m failed" exception=e
end

# ── City labels: black text, no halo ────────────────────────────────────────

label_pts = [
    ("Cairo",        31.24, 30.04),
    ("Alexandria",   29.92, 31.20),
    ("Rosetta",      30.42, 31.40),
    ("Damietta",     31.82, 31.42),
    ("Port Said",    32.28, 31.26),
    ("Ismailia",     32.27, 30.59),
    ("L. Manzala",   31.95, 31.32),
    ("L. Burullus",  30.95, 31.50),
    ("L. Idku",      30.25, 31.28),
    ("Tanta",        30.99, 30.79),
    ("Mansoura",     31.38, 31.04),
    ("Zagazig",      31.52, 30.59),
    ("Banha",        31.18, 30.46),
]
for (name, lon, lat) in label_pts
    scatter!(ax, [lon], [lat]; color=:black, marker=:circle, markersize=5)
    text!(ax, name; position=(lon + 0.05, lat + 0.05),
          fontsize=8, color=:black, font=:regular,
          align=(:left, :bottom))
end

xlims!(ax, bbox.lon_min, bbox.lon_max)
ylims!(ax, bbox.lat_min, bbox.lat_max)

# ── Legend ──────────────────────────────────────────────────────────────────

# Swatches sampled from the actual hypsometric ramp so the legend matches
# what is on the map.  Elevations chosen to hit each stop band.
_sw(e, h) = shade(hypsometric(Float32(e)), Float32(h), 0f0)
elems = [
    PolyElement(color=OCEAN_BLUE),
    PolyElement(color=_sw(-15, 0.75)),
    PolyElement(color=_sw(-3,  0.75)),
    PolyElement(color=_sw(2,   0.75)),
    PolyElement(color=_sw(10,  0.75)),
    PolyElement(color=_sw(30,  0.75)),
    PolyElement(color=_sw(80,  0.75)),
    PolyElement(color=_sw(160, 0.75)),
]
labels = ["Sea (Mediterranean)",
          "Subsidence hotspots (< −10 m)",
          "Reclaimed below sea (−10 to 0 m)",
          "Delta flats (0–5 m)",
          "Delta plain (5–20 m)",
          "Gentle uplands (20–50 m)",
          "Desert margins (50–120 m)",
          "Jebels (> 120 m)"]
Legend(fig[1, 2], elems, labels; framevisible=false, labelsize=9,
       patchsize=(14, 10))
colgap!(fig.layout, 8)

# ── Save ────────────────────────────────────────────────────────────────────

outdir = joinpath(ROOT, "results", "figures")
mkpath(joinpath(outdir, "pdf"))
mkpath(joinpath(outdir, "png"))
name = "nile_delta_dem_shaded"
FigStd.saveboth(joinpath(outdir, name), fig)
@info "Saved" pdf=joinpath(outdir, "pdf", "$name.pdf") png=joinpath(outdir, "png", "$name.png")
