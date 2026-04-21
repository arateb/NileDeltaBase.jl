#!/usr/bin/env julia
# Publication figure of the fused Nile Delta + North Sinai DEM.
#   - Grayscale multidirectional hillshade base
#   - Elevation gradient overlay (lightens high ground)
#   - Below-sea-level land cells in SEA blue
#   - Lakes and lagoons from Natural Earth 10m, in a DIFFERENT (darker) blue
#   - Coastline + country borders, hotspot labels

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

isfile(DEM_TIF) || error("Missing $DEM_TIF — run S06_derived_display.sh first")
isfile(HLS_TIF) || error("Missing $HLS_TIF — run S06_derived_display.sh first")

# ── Raster loading via GDAL CLI (ENVI binary, same pattern as GulfBase) ─────

function _read_raster(tif::String)
    info = read(`gdalinfo -json $tif`, String)
    m_size = match(r"\"size\":\s*\[\s*(\d+),\s*(\d+)\s*\]", info)
    nx = parse(Int, m_size[1]); ny = parse(Int, m_size[2])
    m_geo = match(r"\"geoTransform\":\s*\[([^\]]+)\]", info)
    gt = parse.(Float64, strip.(split(m_geo[1], ',')))
    xmin = gt[1]; dx = gt[2]
    ymax = gt[4]; dy = gt[6]  # negative
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
    ymin = ymax + dy*ny   # dy is negative
    lons = collect(range(xmin + dx/2;    step=dx,  length=nx))
    lats = collect(range(ymin + (-dy)/2; step=-dy, length=ny))
    return arr, lons, lats
end

@info "Loading elevation and hillshade..."
elev, lons, lats = _read_raster(DEM_TIF)
hlsh, _,    _    = _read_raster(HLS_TIF)

valid = .!isnan.(elev)
@info "Grid" size=size(elev) lon=(first(lons), last(lons)) lat=(first(lats), last(lats))
@info "Elevation stats" min=minimum(elev[valid]) max=maximum(elev[valid]) median=median(elev[valid])
@info "Pixels below 0 m: $(sum(valid .& (elev .< 0)))"

# ── Colours ────────────────────────────────────────────────────────────────
const OCEAN_BLUE = RGBAf(0.55, 0.73, 0.88, 1.0)  # open sea / nodata background
const SEA_BLUE   = RGBAf(0.30, 0.56, 0.78, 1.0)  # land below sea level
const LAKE_BLUE  = RGBAf(0.14, 0.32, 0.55, 1.0)  # lakes / lagoons (darker, distinct)

# ── Build the figure ────────────────────────────────────────────────────────

bbox = NDB.BBOX
# Map width/height ratio — figstd's `aspect` is figure height/width.
map_wh = (bbox.lon_max - bbox.lon_min) / (bbox.lat_max - bbox.lat_min)
# Reserve ~22% width for the legend column; figure height/width matches the
# map aspect with that legend margin included.
fig = FigStd.figure(width=:wide, aspect=1.0 / map_wh * 1.22)
ax = Axis(fig[1, 1];
          xlabel="Longitude (°E)", ylabel="Latitude (°N)",
          title="Nile Delta & North Sinai — fused 30 m DEM (COP-GLO30 / FABDEM / DeltaDTM)",
          aspect=DataAspect())

# 1) Ocean background rectangle (covers everything; land painted over it)
poly!(ax, Rect2f(bbox.lon_min, bbox.lat_min,
                 bbox.lon_max - bbox.lon_min,
                 bbox.lat_max - bbox.lat_min);
      color=OCEAN_BLUE, strokewidth=0)

# 2) Grayscale hillshade as base (mask where DEM is nodata so ocean shows)
hlsh_norm = Float32.(hlsh) ./ 255f0
hlsh_norm[.!valid] .= NaN32
heatmap!(ax, lons, lats, permutedims(hlsh_norm);
         colormap=[RGBf(0.22,0.22,0.22), RGBf(0.95,0.95,0.95)],
         nan_color=:transparent,
         colorrange=(0.12, 0.95),
         rasterize=5)

# 3) Transparent elevation overlay — brightens higher ground
elev_vis = copy(elev)
elev_vis[.!valid] .= NaN32
elev_vis[elev .< 0] .= NaN32   # below-sea-level handled separately
heatmap!(ax, lons, lats, permutedims(clamp.(elev_vis, 0f0, 300f0));
         colormap=[RGBAf(1,1,1,0.0), RGBAf(1,1,1,0.55)],
         nan_color=:transparent,
         colorrange=(0, 300),
         rasterize=5)

# 4) Land below sea level in SEA_BLUE
below_mask = Float32.(valid .& (elev .< 0))
below_mask[below_mask .== 0f0] .= NaN32
heatmap!(ax, lons, lats, permutedims(below_mask);
         colormap=[SEA_BLUE, SEA_BLUE],
         nan_color=:transparent,
         colorrange=(0.5, 1.5),
         rasterize=5)

# ── Helpers to plot Natural Earth features as plain polygons/lines.
# poly!(ax, FeatureCollection) on a non-GeoAxis Axis recurses and blows the
# stack; convert to basic geometries and iterate.
function _plot_polys!(ax, fc; color, strokewidth=0.3, strokecolor=:black)
    for feature in fc
        geom = GeoMakie.geo2basic(feature.geometry)
        if geom isa AbstractVector
            for g in geom
                poly!(ax, g; color=color, strokewidth=strokewidth,
                      strokecolor=strokecolor)
            end
        else
            poly!(ax, geom; color=color, strokewidth=strokewidth,
                  strokecolor=strokecolor)
        end
    end
end
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

# 5) Lakes and lagoons — Natural Earth 10m
@info "Loading Natural Earth lakes 10m..."
try
    _plot_polys!(ax, naturalearth("lakes", 10);
                 color=LAKE_BLUE,
                 strokecolor=RGBAf(0.05,0.15,0.35,0.8))
catch e
    @warn "NE lakes 10m failed, trying 50m fallback" exception=e
    try
        _plot_polys!(ax, naturalearth("lakes", 50);
                     color=LAKE_BLUE,
                     strokecolor=RGBAf(0.05,0.15,0.35,0.8))
    catch; end
end

# 6) Coastlines + country borders on top
try
    _plot_lines!(ax, naturalearth("coastline", 10);
                 color=RGBAf(0.05,0.05,0.05,0.9), linewidth=0.5)
catch e
    @warn "coastline 10m failed" exception=e
end
try
    _plot_lines!(ax, naturalearth("admin_0_boundary_lines_land", 10);
                 color=RGBAf(0.15,0.15,0.15,0.8), linewidth=0.4,
                 linestyle=:dash)
catch e
    @warn "borders 10m failed" exception=e
end

# 7) Hotspot labels
label_pts = [
    ("Cairo",       31.24, 30.04),
    ("Alexandria",  29.92, 31.20),
    ("Rosetta",     30.42, 31.40),
    ("Damietta",    31.82, 31.42),
    ("Port Said",   32.28, 31.26),
    ("Ismailia",    32.27, 30.59),
    ("Suez",        32.55, 29.97),
    ("Lake Manzala", 31.95, 31.26),
    ("L. Burullus", 30.95, 31.50),
    ("L. Bardawil", 33.10, 31.15),
    ("El Arish",    33.80, 31.13),
    ("Rafah",       34.24, 31.29),
    ("Taba",        34.89, 29.50),
]
for (name, lon, lat) in label_pts
    scatter!(ax, [lon], [lat]; color=:black, marker=:circle, markersize=6,
             strokewidth=0.6, strokecolor=:white)
    text!(ax, name; position=(lon+0.05, lat+0.05), fontsize=8,
          color=:black, strokecolor=:white, strokewidth=1.5, font=:regular)
end

xlims!(ax, bbox.lon_min, bbox.lon_max)
ylims!(ax, bbox.lat_min, bbox.lat_max)

# Legend
elems = [
    PolyElement(color=OCEAN_BLUE),
    PolyElement(color=SEA_BLUE),
    PolyElement(color=LAKE_BLUE),
    PolyElement(color=RGBf(0.55,0.55,0.55)),
    PolyElement(color=RGBf(0.92,0.92,0.92)),
]
labels = ["Sea (Mediterranean / Red Sea)",
          "Land below sea level (< 0 m)",
          "Lakes & lagoons",
          "Land — low",
          "Land — high"]
Legend(fig[1, 2], elems, labels; framevisible=false, labelsize=9, patchsize=(14,10))
colgap!(fig.layout, 8)

# ── Save ────────────────────────────────────────────────────────────────────

outdir = joinpath(ROOT, "results", "figures")
mkpath(joinpath(outdir, "pdf"))
mkpath(joinpath(outdir, "png"))
name = "nile_delta_dem_shaded"
FigStd.saveboth(joinpath(outdir, name), fig)
@info "Saved" pdf=joinpath(outdir, "pdf", "$name.pdf") png=joinpath(outdir, "png", "$name.png")
