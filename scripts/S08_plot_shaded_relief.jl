#!/usr/bin/env julia
# Nile Delta shaded-relief basemap — GulfBase aesthetic.
#
# Gray terrain ramp (white at sea level → dark gray at high elevation),
# multidirectional hillshade as an intensity multiplier, pale off-white ocean
# from a rasterized GSHHG L1 landmask, and vector overlays (coastline,
# rivers, borders, cities) drawn on top.
#
# Inputs (from S06):
#   nile_elevation_display.tif   0.003° DEM (Float32, -9999 nodata)
#   nile_hillshade_display.tif   0.003° multidirectional hillshade (Byte)
#   nile_landmask_display.tif    0.003° GSHHG L1 landmask (Byte, 1=land)
# Inputs (from S05b):
#   data/raw/shoreline/gshhg_*_nile.shp

using Pkg
Pkg.activate(dirname(@__DIR__))

using Printf, Statistics
using CairoMakie, Shapefile

include("/data/files/pkgs/figstd/src/FigStd.jl")
using .FigStd
FigStd.apply!()

push!(LOAD_PATH, joinpath(dirname(@__DIR__), "src"))
using NileDeltaBase
const NDB = NileDeltaBase
NDB.set_root!(get(ENV, "NILEDELTABASE_ROOT", "/data4/EGY/NileDeltaBase"))

const ROOT     = NDB._root()
const DISPDIR  = joinpath(ROOT, "data/derived/display")
const SHOREDIR = joinpath(ROOT, "data/raw/shoreline")
const DEM_TIF  = joinpath(DISPDIR, "nile_elevation_display.tif")
const HLS_TIF  = joinpath(DISPDIR, "nile_hillshade_display.tif")
const LMSK_TIF = joinpath(DISPDIR, "nile_landmask_display.tif")

for f in (DEM_TIF, HLS_TIF, LMSK_TIF)
    isfile(f) || error("Missing $f — run S06_derived_display.sh first")
end

const BBOX = NDB.BBOX
const DRES = 0.003

# ── Raster I/O (ENVI binary via GDAL CLI) ───────────────────────────────────

function _read_envi(tif, ::Type{T}) where T
    info = read(`gdalinfo -json $tif`, String)
    m = match(r"\"size\":\s*\[\s*(\d+),\s*(\d+)\s*\]", info)
    nx = parse(Int, m[1]); ny = parse(Int, m[2])
    ot  = T == UInt8 ? "Byte" : "Float32"
    bin = tempname() * ".bin"
    run(pipeline(`gdal_translate -of ENVI -ot $ot $tif $bin`; stderr=devnull))
    data = Array{T}(undef, nx, ny)
    read!(bin, data)
    rm(bin; force=true)
    rm(bin * ".aux.xml"; force=true)
    rm(replace(bin, ".bin" => ".hdr"); force=true)
    arr = Matrix{T}(undef, ny, nx)
    for r in 1:ny
        arr[r, :] = data[:, ny - r + 1]
    end
    return arr, nx, ny
end

function _make_coords(nx, ny)
    lons = collect(range(BBOX.lon_min + DRES/2, BBOX.lon_max - DRES/2; length=nx))
    lats = collect(range(BBOX.lat_min + DRES/2, BBOX.lat_max - DRES/2; length=ny))
    return lons, lats
end

function load_dem()
    arr, nx, ny = _read_envi(DEM_TIF, Float32)
    info = read(`gdalinfo -json $DEM_TIF`, String)
    m_nd = match(r"\"noDataValue\":\s*([-\d.eE+]+)", info)
    nd   = m_nd !== nothing ? parse(Float32, m_nd[1]) : -9999f0
    arr[arr .== nd]          .= NaN32
    arr[abs.(arr) .> 1f4]    .= NaN32
    lons, lats = _make_coords(nx, ny)
    return arr, lons, lats
end

function load_hillshade()
    raw_u8, nx, ny = _read_envi(HLS_TIF, UInt8)
    raw = Float32.(raw_u8) ./ 255f0
    # p1/p99 stretch over the (non-zero) valid pixels — keeps escarpment
    # contrast without blowing out the flat Delta.
    valid = filter(x -> x > 0.01f0, raw[:])
    lo, hi = isempty(valid) ? (0f0, 1f0) : (quantile(valid, 0.01), quantile(valid, 0.99))
    span = max(hi - lo, 0.01f0)
    out = similar(raw)
    for i in eachindex(raw)
        out[i] = raw[i] > 0.01f0 ? clamp((raw[i] - lo) / span, 0f0, 1f0) : 0f0
    end
    return out
end

function load_landmask()
    arr, _, _ = _read_envi(LMSK_TIF, UInt8)
    return arr
end

# ── Shapefile vector loading ────────────────────────────────────────────────

function read_shp_segments(shp_path)
    # Shapefile.Handle reads only the .shp geometry (bypasses DBF, which can
    # fail on some GSHHG-derived files after ogr2ogr clipping).
    segs = Vector{Tuple{Vector{Float64}, Vector{Float64}}}()
    try
        h = open(shp_path) do io
            read(io, Shapefile.Handle)
        end
        for geom in h.shapes
            geom === nothing && continue
            pts    = geom.points
            starts = geom.parts
            for k in eachindex(starts)
                i0 = starts[k] + 1
                i1 = k < length(starts) ? starts[k+1] : length(pts)
                lons = [Float64(pts[i].x) for i in i0:i1]
                lats = [Float64(pts[i].y) for i in i0:i1]
                push!(segs, (lons, lats))
            end
        end
    catch e
        @warn "shapefile read failed" shp_path err=sprint(showerror, e)
    end
    return segs
end

function load_rivers_vector()
    rivers = Vector{Tuple{Vector{Float64}, Vector{Float64}, Int}}()
    for (level, n) in [("L01", 1), ("L02", 2), ("L03", 3), ("L04", 4), ("L05", 5)]
        shp = joinpath(SHOREDIR, "gshhg_WDBII_river_f_$(level)_nile.shp")
        isfile(shp) || continue
        for (lo, la) in read_shp_segments(shp)
            push!(rivers, (lo, la, n))
        end
    end
    return rivers
end

function load_borders_vector()
    segs = Vector{Tuple{Vector{Float64}, Vector{Float64}}}()
    for level in ["L1", "L2"]
        shp = joinpath(SHOREDIR, "gshhg_WDBII_border_i_$(level)_nile.shp")
        isfile(shp) || continue
        append!(segs, read_shp_segments(shp))
    end
    return segs
end

function load_coastline_vector()
    shp = joinpath(SHOREDIR, "gshhg_GSHHS_f_L1_nile.shp")
    isfile(shp) || return Vector{Tuple{Vector{Float64}, Vector{Float64}}}()
    return read_shp_segments(shp)
end

# ── Shaded relief composition (ported from GulfBase S18) ────────────────────

function shaded_relief(elev, hs, landmask; vmin=0f0, vmax=150f0)
    ny, nx = size(elev)
    img = Matrix{RGBAf}(undef, ny, nx)
    for j in 1:nx, i in 1:ny
        v = elev[i, j]
        if isnan(v) || landmask[i, j] == 0x00
            img[i, j] = RGBAf(0.96, 0.96, 0.98, 1.0)           # pale off-white water
        else
            t = clamp((v - vmin) / (vmax - vmin), 0f0, 1f0)
            if v < 0
                # Below-sea-level polders (Manzala / Burullus interiors) —
                # very faint cool tint so the viewer can still see the extent.
                r = 0.90f0; g = 0.92f0; b = 0.96f0
            else
                g_base = Float32(1.0 - 0.82 * t)                # white → dark gray
                r = g_base; g = g_base; b = g_base
            end
            h = 0.35f0 + 0.65f0 * hs[i, j]                      # hillshade multiplier
            img[i, j] = RGBAf(r * h, g * h, b * h, 1.0)
        end
    end
    return img
end

# ── Vector overlays ─────────────────────────────────────────────────────────

const LON_RANGE = (BBOX.lon_min, BBOX.lon_max)
const LAT_RANGE = (BBOX.lat_min, BBOX.lat_max)

function _in_bbox(lo, la)
    any(lo[i] >= LON_RANGE[1] && lo[i] <= LON_RANGE[2] &&
        la[i] >= LAT_RANGE[1] && la[i] <= LAT_RANGE[2] for i in eachindex(lo))
end

function add_coastline!(ax, coastline; color=RGBAf(0.15, 0.15, 0.15, 0.7), linewidth=0.5)
    for (lo, la) in coastline
        _in_bbox(lo, la) || continue
        lines!(ax, lo, la; color=color, linewidth=linewidth)
    end
end

function add_rivers!(ax, rivers; color=RGBAf(0.20, 0.40, 0.70, 0.75))
    lw_map = Dict(1 => 1.2, 2 => 1.0, 3 => 0.7, 4 => 0.4, 5 => 0.25)
    for (lo, la, lev) in rivers
        _in_bbox(lo, la) || continue
        lines!(ax, lo, la; color=color, linewidth=get(lw_map, lev, 0.3))
    end
end

function add_borders!(ax, borders; color=RGBAf(0.25, 0.25, 0.25, 0.55), linewidth=0.6)
    for (lo, la) in borders
        _in_bbox(lo, la) || continue
        lines!(ax, lo, la; color=color, linewidth=linewidth, linestyle=:dash)
    end
end

function add_cities!(ax, cities; fontsize=7, markersize=5)
    for (name, lon, lat) in cities
        (LON_RANGE[1] <= lon <= LON_RANGE[2] && LAT_RANGE[1] <= lat <= LAT_RANGE[2]) || continue
        scatter!(ax, [lon], [lat]; color=:black, markersize=markersize,
                 marker=:circle, strokewidth=0)
        text!(ax, lon + 0.04, lat + 0.04; text=name, fontsize=fontsize,
              color=:black, align=(:left, :bottom))
    end
end

function add_lakes!(ax, lakes; fontsize=7)
    for (name, lon, lat) in lakes
        (LON_RANGE[1] <= lon <= LON_RANGE[2] && LAT_RANGE[1] <= lat <= LAT_RANGE[2]) || continue
        text!(ax, lon, lat; text=name, fontsize=fontsize,
              color=RGBAf(0.10, 0.20, 0.45, 0.85),
              font=:italic, align=(:center, :center))
    end
end

# ── Nile Delta cities & lagoons ─────────────────────────────────────────────

const CITIES = [
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
    ("Kafr el-Sheikh", 30.94, 31.11),
]

const LAKES = [
    ("L. Idku",     30.25, 31.26),
    ("L. Burullus", 30.95, 31.49),
    ("L. Manzala",  31.95, 31.32),
    ("L. Mariout",  29.93, 31.10),
]

# ── Figure ──────────────────────────────────────────────────────────────────

function plot_shaded_basemap()
    @info "Loading elevation, hillshade, landmask..."
    elev, lons, lats = load_dem()
    hs   = load_hillshade()
    lmsk = load_landmask()

    @info "Loading vector overlays..."
    rivers    = load_rivers_vector()
    borders   = load_borders_vector()
    coastline = load_coastline_vector()

    @info "Compositing shaded relief" size=(size(elev,1), size(elev,2)) rivers=length(rivers) borders=length(borders)
    # Nile Delta elevation range: -23 → 700 m (desert escarpments rise
    # sharply at the Delta's margin).  vmax=200 keeps the Delta at near-white
    # (flat ≤10 m) while giving the surrounding desert clear gray gradation.
    img = shaded_relief(elev, hs, lmsk; vmin=0f0, vmax=200f0)

    map_wh = (BBOX.lon_max - BBOX.lon_min) / (BBOX.lat_max - BBOX.lat_min)
    fig = FigStd.figure(width=:wide, aspect=1.0 / map_wh * 1.15)
    ax  = Axis(fig[1, 1];
               title="Nile Delta Shaded Relief — fused 30 m DEM",
               xlabel="Longitude (°E)", ylabel="Latitude (°N)",
               aspect=DataAspect())

    image!(ax, (lons[1], lons[end]), (lats[1], lats[end]),
           permutedims(img); interpolate=false)

    add_coastline!(ax, coastline)
    add_rivers!(ax,    rivers)
    add_borders!(ax,   borders)
    add_cities!(ax,    CITIES; fontsize=7, markersize=4)
    add_lakes!(ax,     LAKES;  fontsize=7)

    xlims!(ax, BBOX.lon_min, BBOX.lon_max)
    ylims!(ax, BBOX.lat_min, BBOX.lat_max)

    cb_cmap = cgrad([RGBf(1,1,1), RGBf(0.18,0.18,0.18)])
    Colorbar(fig[1, 2]; colormap=cb_cmap, limits=(0, 200),
             label="Elevation (m)", width=12)
    colgap!(fig.layout, 8)

    outdir = joinpath(ROOT, "results", "figures")
    mkpath(joinpath(outdir, "pdf"))
    mkpath(joinpath(outdir, "png"))
    name = "nile_delta_shaded_relief"
    FigStd.saveboth(joinpath(outdir, name), fig)
    @info "Saved" pdf=joinpath(outdir, "pdf", "$name.pdf") png=joinpath(outdir, "png", "$name.png")
end

plot_shaded_basemap()
