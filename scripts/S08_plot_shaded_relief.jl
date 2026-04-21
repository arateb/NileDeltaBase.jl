#!/usr/bin/env julia
# Classic grayscale shaded-relief map of the Nile Delta.
# Hillshade is the hero — no hypsometric tint.  Water carries the only
# colour (ocean, lakes, below-sea-level land), so the eye reads topography
# as shading only.

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

for f in (DEM_TIF, HLS_TIF, SLP_TIF)
    isfile(f) || error("Missing $f — run S06_derived_display.sh first")
end

# ── Raster loader (ENVI binary via GDAL CLI) ────────────────────────────────

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

valid   = .!isnan.(elev)
invalid = .!valid
below   = valid .& (elev .< 0f0)

# Hillshade percentile stretch — same rationale as S07.
hlsh_vals = filter(!isnan, hlsh[valid])
h_p2  = Float32(quantile(hlsh_vals, 0.02))
h_p98 = Float32(quantile(hlsh_vals, 0.98))
h_range = max(h_p98 - h_p2, 1f0)
hlsh_n  = clamp.((hlsh .- h_p2) ./ h_range, 0f0, 1f0)

# Slope-darkening term.
slp_n = clamp.(slp ./ 35f0, 0f0, 1f0)
slp_n[invalid] .= 0f0

@info "Grid" size=(ny, nx) lon=(first(lons), last(lons)) lat=(first(lats), last(lats))
@info "Pixels below 0 m: $(sum(below))"

# ── Flood-fill from boundary to separate ocean from enclosed lakes ──────────

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

# ── Palette ────────────────────────────────────────────────────────────────
# Shaded relief: land is pure grayscale driven by hillshade; only water
# carries chroma.

const OCEAN_BLUE = RGBf(0.58, 0.75, 0.88)   # open sea
const LAKE_BLUE  = RGBf(0.22, 0.42, 0.66)   # enclosed lakes / lagoons
const BELOW_BLUE = RGBf(0.80, 0.87, 0.92)   # below-sea-level land (very pale)

# Grayscale land: lightness ∈ [0.35, 1.0], driven entirely by stretched
# hillshade minus a slope-darkening term.  No elevation colouring.
@inline function land_gray(hn::Float32, sn::Float32)
    hn = isnan(hn) ? 0.5f0 : hn
    # Sun on slopes facing the light are nearly white; shaded slopes + steep
    # terrain darken to ~0.35.  The 0.82 coefficient keeps the brightest
    # sunlit faces from saturating to pure white.
    g = 0.38f0 + 0.82f0 * hn - 0.40f0 * sn
    g = clamp(g, 0.12f0, 1.0f0)
    RGBf(g, g, g)
end

# Blend the below-sea-level blue with hillshade so relief still reads
# inside the polder areas.
@inline function below_shade(hn::Float32, sn::Float32)
    hn = isnan(hn) ? 0.5f0 : hn
    g = 0.55f0 + 0.55f0 * hn - 0.35f0 * sn
    g = clamp(g, 0.35f0, 1.05f0)
    RGBf(clamp(BELOW_BLUE.r * g, 0f0, 1f0),
         clamp(BELOW_BLUE.g * g, 0f0, 1f0),
         clamp(BELOW_BLUE.b * g, 0f0, 1f0))
end

# ── Single RGB image ────────────────────────────────────────────────────────

rgb = Matrix{RGBf}(undef, ny, nx)
@inbounds for i in 1:ny, j in 1:nx
    if ocean[i, j]
        rgb[i, j] = OCEAN_BLUE
    elseif lake[i, j]
        rgb[i, j] = LAKE_BLUE
    elseif below[i, j]
        rgb[i, j] = below_shade(hlsh_n[i, j], slp_n[i, j])
    else
        rgb[i, j] = land_gray(hlsh_n[i, j], slp_n[i, j])
    end
end

# ── Figure ──────────────────────────────────────────────────────────────────

bbox = NDB.BBOX
map_wh = (bbox.lon_max - bbox.lon_min) / (bbox.lat_max - bbox.lat_min)
fig = FigStd.figure(width=:wide, aspect=1.0 / map_wh * 1.22)
ax = Axis(fig[1, 1];
          xlabel="Longitude (°E)", ylabel="Latitude (°N)",
          title="Nile Delta — shaded relief (fused 30 m DEM)",
          aspect=DataAspect())

image!(ax, (first(lons), last(lons)), (first(lats), last(lats)),
       permutedims(rgb); interpolate=false)

# Country borders (no NE coastline — DEM defines the coast).
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
                 color=RGBAf(0.10, 0.10, 0.10, 0.75),
                 linewidth=0.5, linestyle=:dash)
catch e
    @warn "borders 10m failed" exception=e
end

# City & lagoon labels (plain black text, no halo).
label_pts = [
    ("Cairo",       31.24, 30.04),
    ("Alexandria",  29.92, 31.20),
    ("Rosetta",     30.42, 31.40),
    ("Damietta",    31.82, 31.42),
    ("Port Said",   32.28, 31.26),
    ("Ismailia",    32.27, 30.59),
    ("Tanta",       30.99, 30.79),
    ("Mansoura",    31.38, 31.04),
    ("Zagazig",     31.52, 30.59),
    ("Banha",       31.18, 30.46),
    ("L. Idku",     30.25, 31.28),
    ("L. Burullus", 30.95, 31.50),
    ("L. Manzala",  31.95, 31.32),
]
for (name, lon, lat) in label_pts
    scatter!(ax, [lon], [lat]; color=:black, marker=:circle, markersize=5)
    text!(ax, name; position=(lon + 0.04, lat + 0.04),
          fontsize=8, color=:black, font=:regular,
          align=(:left, :bottom))
end

xlims!(ax, bbox.lon_min, bbox.lon_max)
ylims!(ax, bbox.lat_min, bbox.lat_max)

# Legend
elems = [
    PolyElement(color=OCEAN_BLUE),
    PolyElement(color=LAKE_BLUE),
    PolyElement(color=below_shade(0.75f0, 0f0)),
    PolyElement(color=land_gray(0.50f0, 0.20f0)),   # shaded slope
    PolyElement(color=land_gray(0.85f0, 0f0)),      # sunlit slope / flat
]
labels = ["Sea (Mediterranean)",
          "Lakes & lagoons",
          "Land below sea level",
          "Shaded slope",
          "Sunlit / flat land"]
Legend(fig[1, 2], elems, labels;
       framevisible=false, labelsize=9, patchsize=(14, 10))
colgap!(fig.layout, 8)

# ── Save ────────────────────────────────────────────────────────────────────

outdir = joinpath(ROOT, "results", "figures")
mkpath(joinpath(outdir, "pdf"))
mkpath(joinpath(outdir, "png"))
name = "nile_delta_shaded_relief"
FigStd.saveboth(joinpath(outdir, name), fig)
@info "Saved" pdf=joinpath(outdir, "pdf", "$name.pdf") png=joinpath(outdir, "png", "$name.png")
