#!/usr/bin/env julia
# Extended Nile basemap — shaded relief + ESA WorldCover overlays, with a
# zoom-region dispatcher so the same composition can be rendered for the full
# extent (Dabaa → Port Said) and for focused panels (Dabaa, Alexandria, Cairo,
# Rosetta, Damietta, Port Said, Ismailia, Manzala, Burullus, …).
#
# Composition:
#   • Gray elevation ramp (white at sea level → dark gray at high elevation)
#   • Multidirectional hillshade as an intensity multiplier
#   • WorldCover cropland (class 40) → olive-green tint
#   • WorldCover built-up (class 50) → warm brown tint
#   • WorldCover tree/shrub/grass (10/20/30) → faint green tint
#   • WorldCover permanent water (80) & herbaceous wetland (90) → blue
#   • Pale off-white ocean from rasterized GSHHG L1 landmask
#   • GSHHG vectors on top (coast, rivers, borders)
#
# Inputs (from S06b + S05b + S05c):
#   nile_elevation_display_ext.tif   0.003° DEM
#   nile_hillshade_display_ext.tif   hillshade
#   nile_landmask_display_ext.tif    GSHHG L1 land = 1
#   nile_worldcover_display_ext.tif  WorldCover class codes
#   data/raw/shoreline/gshhg_*_ext.shp

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

# Extended bbox (matches S06b output, S05c clip, NDB.BBOX_EXT).
const BBOX_EXT = NDB.BBOX_EXT
const DRES = 0.003

# ── Shapefile vector loading (Handle: skip DBF) ─────────────────────────────

function read_shp_segments(shp_path)
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
        shp = joinpath(SHOREDIR, "gshhg_WDBII_river_f_$(level)_ext.shp")
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
        shp = joinpath(SHOREDIR, "gshhg_WDBII_border_i_$(level)_ext.shp")
        isfile(shp) || continue
        append!(segs, read_shp_segments(shp))
    end
    return segs
end

function load_coastline_vector()
    shp = joinpath(SHOREDIR, "gshhg_GSHHS_f_L1_ext.shp")
    isfile(shp) || return Vector{Tuple{Vector{Float64}, Vector{Float64}}}()
    return read_shp_segments(shp)
end

# ── Basemap composition now lives in NileDeltaBase.Basemap (reusable). ──────
# Pull the WorldCover tints locally for the legend.
using NileDeltaBase.Basemap: WC_CROPLAND, WC_BUILT, WC_TREE, WC_SHRUB,
                             WC_GRASS, WC_WATER, WC_WETLAND, compose_basemap

# ── Overlays ────────────────────────────────────────────────────────────────

function _in_range(lo, la, lon_range, lat_range)
    any(lo[i] >= lon_range[1] && lo[i] <= lon_range[2] &&
        la[i] >= lat_range[1] && la[i] <= lat_range[2] for i in eachindex(lo))
end

function add_coastline!(ax, coastline, lon_range, lat_range;
                        color=RGBAf(0.15, 0.15, 0.15, 0.7), linewidth=0.5)
    for (lo, la) in coastline
        _in_range(lo, la, lon_range, lat_range) || continue
        lines!(ax, lo, la; color=color, linewidth=linewidth)
    end
end

function add_rivers!(ax, rivers, lon_range, lat_range;
                     color=RGBAf(0.15, 0.35, 0.70, 0.85))
    lw_map = Dict(1 => 1.2, 2 => 1.0, 3 => 0.7, 4 => 0.4, 5 => 0.25)
    for (lo, la, lev) in rivers
        _in_range(lo, la, lon_range, lat_range) || continue
        lines!(ax, lo, la; color=color, linewidth=get(lw_map, lev, 0.3))
    end
end

function add_borders!(ax, borders, lon_range, lat_range;
                      color=RGBAf(0.25, 0.25, 0.25, 0.55), linewidth=0.6)
    for (lo, la) in borders
        _in_range(lo, la, lon_range, lat_range) || continue
        lines!(ax, lo, la; color=color, linewidth=linewidth, linestyle=:dash)
    end
end

function add_labels!(ax, pts, lon_range, lat_range;
                     fontsize=7, markersize=4, color=:black,
                     font=:regular, use_markers=true)
    for (name, lon, lat) in pts
        (lon_range[1] <= lon <= lon_range[2] && lat_range[1] <= lat <= lat_range[2]) || continue
        if use_markers
            scatter!(ax, [lon], [lat]; color=:black, markersize=markersize,
                     marker=:circle, strokewidth=0)
            text!(ax, lon + 0.04 * (lon_range[2] - lon_range[1]) / 4.5,
                  lat + 0.04 * (lat_range[2] - lat_range[1]) / 1.7;
                  text=name, fontsize=fontsize, color=color, font=font,
                  align=(:left, :bottom))
        else
            text!(ax, lon, lat; text=name, fontsize=fontsize, color=color, font=font,
                  align=(:center, :center))
        end
    end
end

# ── Cities, lagoons & key sites ─────────────────────────────────────────────

const CITIES = [
    ("Cairo",         31.24, 30.04),
    ("Alexandria",    29.92, 31.20),
    ("Rosetta",       30.42, 31.40),
    ("Damietta",      31.82, 31.42),
    ("Port Said",     32.28, 31.26),
    ("Ismailia",      32.27, 30.59),
    ("Tanta",         30.99, 30.79),
    ("Mansoura",      31.38, 31.04),
    ("Zagazig",       31.52, 30.59),
    ("Banha",         31.18, 30.46),
    ("Kafr el-Sheikh",30.94, 31.11),
    ("Dabaa",         28.48, 31.04),
    ("Marsa Matruh*", 28.05, 31.35),
    ("El Alamein",    28.96, 30.83),
]

const LAKES = [
    ("L. Idku",     30.25, 31.26),
    ("L. Burullus", 30.95, 31.49),
    ("L. Manzala",  31.95, 31.32),
    ("L. Mariout",  29.93, 31.10),
]

# ── Zoom regions ────────────────────────────────────────────────────────────
#
# Each entry: bbox + human-readable title + colorbar vmax (uplands need
# higher vmax to avoid everything clamping to black).

const ZOOMS = Dict(
    :full       => (bbox=BBOX_EXT,                                              title="Nile Delta — Dabaa → Port Said",            vmax=200),
    :dabaa      => (bbox=(lon_min=28.10, lat_min=30.80, lon_max=28.90, lat_max=31.30), title="Dabaa coast",                                vmax=150),
    :alexandria => (bbox=(lon_min=29.70, lat_min=31.10, lon_max=30.20, lat_max=31.35), title="Alexandria & L. Mariout",                    vmax=80),
    :rosetta    => (bbox=(lon_min=30.25, lat_min=31.30, lon_max=30.70, lat_max=31.55), title="Rosetta promontory",                         vmax=60),
    :burullus   => (bbox=(lon_min=30.45, lat_min=31.30, lon_max=31.25, lat_max=31.60), title="Lake Burullus",                              vmax=60),
    :damietta   => (bbox=(lon_min=31.55, lat_min=31.30, lon_max=31.95, lat_max=31.55), title="Damietta promontory",                        vmax=60),
    :manzala    => (bbox=(lon_min=31.55, lat_min=31.00, lon_max=32.35, lat_max=31.45), title="Lake Manzala",                               vmax=60),
    :port_said  => (bbox=(lon_min=32.10, lat_min=31.10, lon_max=32.45, lat_max=31.35), title="Port Said & Suez Canal mouth",               vmax=60),
    :ismailia   => (bbox=(lon_min=32.10, lat_min=30.40, lon_max=32.45, lat_max=30.80), title="Ismailia & Suez Canal",                      vmax=100),
    :cairo      => (bbox=(lon_min=30.90, lat_min=29.95, lon_max=31.60, lat_max=30.35), title="Greater Cairo & Delta apex",                 vmax=200),
    :central    => (bbox=(lon_min=30.70, lat_min=30.50, lon_max=31.70, lat_max=31.20), title="Central Delta — Tanta / Mansoura / Zagazig", vmax=60),
)

# ── Plot dispatcher ─────────────────────────────────────────────────────────

function subset(data::Matrix, lons::Vector, lats::Vector, bbox)
    ic = findall(x -> bbox.lon_min <= x <= bbox.lon_max, lons)
    ir = findall(x -> bbox.lat_min <= x <= bbox.lat_max, lats)
    isempty(ic) && error("No lon pixels in bbox")
    isempty(ir) && error("No lat pixels in bbox")
    return data[ir, ic], lons[ic], lats[ir]
end

function plot_region(region::Symbol;
                     elev, hs, lmsk, wc, bathy, lons, lats,
                     rivers, borders, coastline,
                     outdir="results/figures",
                     show_cropland_legend=true)
    info = get(ZOOMS, region, nothing)
    info === nothing && error("Unknown region $region — options: $(collect(keys(ZOOMS)))")

    bbox  = info.bbox
    title = info.title
    vmax  = Float32(info.vmax)

    e_sub,  lo, la = subset(elev, lons, lats, bbox)
    h_sub,  _,  _  = subset(hs,   lons, lats, bbox)
    l_sub,  _,  _  = subset(lmsk, lons, lats, bbox)
    w_sub,  _,  _  = subset(wc,   lons, lats, bbox)
    b_sub = bathy === nothing ? nothing : subset(bathy, lons, lats, bbox)[1]

    img = compose_basemap(e_sub, h_sub, l_sub, w_sub;
                           bathy=b_sub, vmin=0f0, vmax=vmax)

    map_wh = (bbox.lon_max - bbox.lon_min) / (bbox.lat_max - bbox.lat_min)
    fig = FigStd.figure(width=:wide, aspect=1.0 / map_wh * 1.12)
    ax  = Axis(fig[1, 1];
               title=title,
               xlabel="Longitude (°E)", ylabel="Latitude (°N)",
               aspect=DataAspect())

    image!(ax, (lo[1], lo[end]), (la[1], la[end]),
           permutedims(img); interpolate=false)

    lon_range = (bbox.lon_min, bbox.lon_max)
    lat_range = (bbox.lat_min, bbox.lat_max)

    add_coastline!(ax, coastline, lon_range, lat_range)
    add_rivers!(ax,    rivers,    lon_range, lat_range)
    add_borders!(ax,   borders,   lon_range, lat_range)
    add_labels!(ax, CITIES, lon_range, lat_range; fontsize=7, markersize=4)
    add_labels!(ax, LAKES,  lon_range, lat_range; fontsize=7,
                color=RGBAf(0.10, 0.20, 0.45, 0.85), font=:italic, use_markers=false)

    xlims!(ax, lon_range...)
    ylims!(ax, lat_range...)

    # Elevation colorbar
    cb_cmap = cgrad([RGBf(1,1,1), RGBf(0.18,0.18,0.18)])
    Colorbar(fig[1, 2]; colormap=cb_cmap, limits=(0, vmax),
             label="Elevation (m)", width=12)

    # WorldCover legend
    elems = [
        PolyElement(color=WC_CROPLAND),
        PolyElement(color=WC_BUILT),
        PolyElement(color=WC_TREE),
        PolyElement(color=WC_GRASS),
        PolyElement(color=WC_WATER),
        PolyElement(color=WC_WETLAND),
    ]
    labels = ["Cropland", "Built-up", "Trees", "Grassland", "Water", "Wetland"]
    Legend(fig[2, 1:2], elems, labels; orientation=:horizontal,
           framevisible=false, labelsize=8, patchsize=(14, 10), nbanks=1)
    rowsize!(fig.layout, 2, Fixed(28))
    colgap!(fig.layout, 8)

    name = "nile_basemap_$(region)"
    mkpath(joinpath(outdir, "pdf"))
    mkpath(joinpath(outdir, "png"))
    FigStd.saveboth(joinpath(outdir, name), fig)
    @info "Saved" region name png=joinpath(outdir, "png", "$name.png")
    return fig
end

function main()
    @info "Loading basemap inputs (elev/hillshade/landmask/WC/bathymetry)..."
    bi = NDB.load_basemap(; bbox=BBOX_EXT, dres=DRES, with_bathy=true)
    elev, hs, lmsk, wc, bathy = bi.elev, bi.hs, bi.landmask, bi.wc, bi.bathy
    lons, lats = bi.lons, bi.lats
    bathy === nothing || @info "bathymetry loaded" size=size(bathy)

    @info "Loading vector overlays..."
    rivers    = load_rivers_vector()
    borders   = load_borders_vector()
    coastline = load_coastline_vector()

    outdir = joinpath(ROOT, "results", "figures")
    regions = [:full, :dabaa, :alexandria, :rosetta, :burullus, :damietta,
               :manzala, :port_said, :ismailia, :cairo, :central]
    @info "Rendering" n=length(regions) regions
    for r in regions
        plot_region(r; elev=elev, hs=hs, lmsk=lmsk, wc=wc, bathy=bathy,
                       lons=lons, lats=lats,
                       rivers=rivers, borders=borders, coastline=coastline,
                       outdir=outdir)
    end
    @info "All basemaps saved" outdir=joinpath(outdir, "png")
end

main()
