#!/usr/bin/env julia
# Publication figure of the fused Nile Delta + North Sinai DEM.
# Single RGB image built as hypsometric tint × hillshade:
#   - ocean (NaN pixels connected to the bbox boundary) → OCEAN_BLUE
#   - enclosed NaN (lakes / lagoons)                    → LAKE_BLUE
#   - land with elevation < 0                           → SEA_BLUE tinted
#   - land with elevation ≥ 0                           → terrain palette
# Coastline is therefore defined by the DEM itself (no Natural Earth overlay).

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
below    = valid .& (elev .< 0f0)

@info "Grid" size=(ny, nx) lon=(first(lons), last(lons)) lat=(first(lats), last(lats))
@info "Elevation stats" min=minimum(elev[valid]) max=maximum(elev[valid]) median=median(elev[valid])
@info "Hillshade stats" min=minimum(hlsh[valid]) max=maximum(hlsh[valid]) median=median(hlsh[valid])
@info "Pixels below 0 m: $(sum(below))"

# Percentile stretch of the hillshade — most of the bbox is the ~0-slope
# Delta, which compresses the raw gdaldem distribution into a narrow band
# near ~180.  Stretch p2..p98 across [0, 1] so the genuine variation (Sinai
# relief, Delta escarpment, canal levees) shows up on the map.
hlsh_vals = filter(!isnan, hlsh[valid])
h_p2  = Float32(quantile(hlsh_vals, 0.02))
h_p98 = Float32(quantile(hlsh_vals, 0.98))
h_range = max(h_p98 - h_p2, 1f0)
hlsh_n  = clamp.((hlsh .- h_p2) ./ h_range, 0f0, 1f0)
@info "Hillshade stretch" p2=h_p2 p98=h_p98 range=h_range

# Slope-darkening multiplier — steep slopes (≳ 15°) get visibly darker even
# after hillshade stretch.  Slopes are in degrees; clamp to [0, 35] and map
# to an extra attenuation in [0, 0.35].
slp_n  = clamp.(slp ./ 35f0, 0f0, 1f0)
slp_n[invalid] .= 0f0

# ── Classify water: ocean = NaN connected to bbox boundary; lake = enclosed NaN

function _flood_boundary!(ocean::AbstractMatrix{Bool}, invalid::AbstractMatrix{Bool})
    ny, nx = size(invalid)
    q = Tuple{Int,Int}[]
    for j in 1:nx
        if invalid[1,  j]; ocean[1,  j] = true; push!(q, (1,  j)); end
        if invalid[ny, j]; ocean[ny, j] = true; push!(q, (ny, j)); end
    end
    for i in 1:ny
        if invalid[i, 1 ]; ocean[i, 1 ] = true; push!(q, (i, 1 )); end
        if invalid[i, nx]; ocean[i, nx] = true; push!(q, (i, nx)); end
    end
    head = 1
    while head <= length(q)
        i, j = q[head]; head += 1
        for (di, dj) in ((-1,0),(1,0),(0,-1),(0,1))
            ni, nj = i+di, j+dj
            if 1 <= ni <= ny && 1 <= nj <= nx && invalid[ni,nj] && !ocean[ni,nj]
                ocean[ni,nj] = true
                push!(q, (ni, nj))
            end
        end
    end
end

ocean = falses(ny, nx)
_flood_boundary!(ocean, invalid)
lake = invalid .& .!ocean
@info "Water classes" ocean=sum(ocean) lake=sum(lake) below_sea_land=sum(below)

# ── Colour palette ──────────────────────────────────────────────────────────

const OCEAN_BLUE = RGBf(0.55, 0.73, 0.88)
const SEA_BLUE   = RGBf(0.33, 0.57, 0.80)   # land below sea level (SEA_BLUE tint)
const LAKE_BLUE  = RGBf(0.12, 0.28, 0.52)   # lakes / lagoons (darker, distinct)

# Hypsometric tint — pale green coastal → tan → brown → gray peaks
const HYPSO_STOPS = [
    (0f0,    RGBf(0.82, 0.90, 0.72)),   # coastal green
    (20f0,   RGBf(0.92, 0.93, 0.68)),   # pale yellow
    (100f0,  RGBf(0.90, 0.82, 0.58)),   # tan
    (300f0,  RGBf(0.78, 0.64, 0.42)),   # brown
    (700f0,  RGBf(0.58, 0.44, 0.30)),   # dark brown
    (1300f0, RGBf(0.78, 0.72, 0.68)),   # gray peaks
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

# Hillshade-derived intensity in [0.35, 1.15] (values > 1 get clamped after
# channel multiply).  Steep slopes further attenuate by up to 0.35, so
# Sinai escarpments and ridgelines read as proper shadows.
@inline function shade(c::RGBf, hn::Float32, sn::Float32)
    hn = isnan(hn) ? 0.5f0 : hn
    s  = 0.35f0 + 0.80f0 * hn - 0.35f0 * sn
    s  = clamp(s, 0f0, 1.15f0)
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
    elseif below[i, j]
        rgb[i, j] = shade(SEA_BLUE, hlsh_n[i, j], slp_n[i, j])
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

# Country borders only (coastline comes from the DEM mask itself).  Convert
# FeatureCollection → basic geometries so plotting works on a plain Axis.
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
# what's on the map (approximating the low/high land entries).
elev_low  = shade(hypsometric(5f0),   0.75f0, 0f0)
elev_mid  = shade(hypsometric(100f0), 0.75f0, 0f0)
elev_high = shade(hypsometric(700f0), 0.75f0, 0f0)

elems = [
    PolyElement(color=OCEAN_BLUE),
    PolyElement(color=LAKE_BLUE),
    PolyElement(color=shade(SEA_BLUE, 0.75f0, 0f0)),
    PolyElement(color=elev_low),
    PolyElement(color=elev_mid),
    PolyElement(color=elev_high),
]
labels = ["Sea (Mediterranean / Red Sea)",
          "Lakes & lagoons",
          "Land below sea level (< 0 m)",
          "Delta (< 20 m)",
          "Uplands (20–300 m)",
          "Mountains (300–1200 m)"]
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
