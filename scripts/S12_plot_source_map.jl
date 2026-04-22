#!/usr/bin/env julia
# Pixel-provenance map — shows which fused layer wins at every location.
# Generated from NileDeltaBase.resolution_source_map at the extended bbox.
#
# Colour key:
#   0 nodata (white)
#   1 GMRT bathymetry / topography (blue)
#   2 COP-GLO30 DSM (tan)
#   3 FABDEM bare-earth (olive)
#   4 DeltaDTM coastal (green)
#   5 TanDEM-X 12 m (purple, future)

using Pkg
Pkg.activate(dirname(@__DIR__))

using CairoMakie

include("/data/files/pkgs/figstd/src/FigStd.jl")
using .FigStd
FigStd.apply!()

push!(LOAD_PATH, joinpath(dirname(@__DIR__), "src"))
using NileDeltaBase
const NDB = NileDeltaBase
NDB.set_root!(get(ENV, "NILEDELTABASE_ROOT", "/data4/EGY/NileDeltaBase"))

const ROOT   = NDB._root()
const OUTDIR = joinpath(ROOT, "results", "figures")

# Build source map on the extended bbox at 0.01° so hand-ops stay fast.
# resolution_source_map uses the narrow BBOX internally; call the generic
# warp pipeline directly for the ext bbox.

const BBOX = (lon_min=28.0, lat_min=30.0, lon_max=32.5, lat_max=31.7)
const DRES = 0.01

function build_source_map()
    nx = round(Int, (BBOX.lon_max - BBOX.lon_min) / DRES)
    ny = round(Int, (BBOX.lat_max - BBOX.lat_min) / DRES)
    src = zeros(Int8, ny, nx)
    layers = [
        (:gmrt,      Int8(1)),
        (:cop_glo30, Int8(2)),
        (:fabdem,    Int8(3)),
        (:deltadtm,  Int8(4)),
        (:tandemx12, Int8(5)),
    ]
    for (layer, code) in layers
        vrt = layer === :deltadtm ?
              joinpath(ROOT, "data/vrt/deltadtm_nile_masked.vrt") :
              NDB.vrt_path(layer)
        isfile(vrt) || continue
        bin = tempname() * ".bin"
        try
            run(pipeline(`gdalwarp -overwrite
                -t_srs EPSG:4326
                -te $(BBOX.lon_min) $(BBOX.lat_min) $(BBOX.lon_max) $(BBOX.lat_max)
                -tr $DRES $DRES
                -r near -wo NUM_THREADS=ALL_CPUS
                -of ENVI -ot Float32 -dstnodata -9999
                $vrt $bin`; stderr=devnull))
            info = read(`gdalinfo -json $bin`, String)
            m = match(r"\"size\":\s*\[\s*(\d+),\s*(\d+)\s*\]", info)
            nnx = parse(Int, m[1]); nny = parse(Int, m[2])
            raw = Array{Float32}(undef, nnx, nny)
            read!(bin, raw)
            for c in 1:min(nnx, nx), r in 1:min(nny, ny)
                v = raw[c, nny - r + 1]
                if v != -9999f0 && abs(v) < 1f4
                    src[r, c] = code
                end
            end
        catch; end
        rm(bin; force=true); rm(bin * ".aux.xml"; force=true)
        rm(replace(bin, ".bin" => ".hdr"); force=true)
    end
    lons = collect(range(BBOX.lon_min + DRES/2, BBOX.lon_max - DRES/2; length=nx))
    lats = collect(range(BBOX.lat_min + DRES/2, BBOX.lat_max - DRES/2; length=ny))
    return src, lons, lats
end

function main()
    @info "Computing source provenance..."
    src, lons, lats = build_source_map()

    # Report coverage fractions
    total = length(src)
    for code in 0:5
        n = sum(src .== code)
        label = ("nodata", "GMRT bathy/topo", "COP-GLO30", "FABDEM",
                 "DeltaDTM", "TanDEM-X 12m")[code + 1]
        @info "coverage" label code n pct=round(100*n/total, digits=2)
    end

    colors = [
        RGBAf(1.00, 1.00, 1.00, 0.0),   # 0 nodata
        RGBAf(0.25, 0.45, 0.75, 0.85),  # 1 GMRT
        RGBAf(0.85, 0.70, 0.45, 0.85),  # 2 GLO30 (tan)
        RGBAf(0.55, 0.65, 0.30, 0.85),  # 3 FABDEM (olive)
        RGBAf(0.30, 0.65, 0.35, 0.85),  # 4 DeltaDTM (green)
        RGBAf(0.55, 0.30, 0.70, 0.85),  # 5 TanDEM-X (purple)
    ]
    labels = ["nodata", "GMRT bathy/topo", "COP-GLO30 DSM",
              "FABDEM bare-earth", "DeltaDTM coastal", "TanDEM-X 12 m"]

    rgb = Matrix{RGBAf}(undef, size(src)...)
    for j in axes(src, 2), i in axes(src, 1)
        rgb[i, j] = colors[src[i, j] + 1]
    end

    map_wh = (BBOX.lon_max - BBOX.lon_min) / (BBOX.lat_max - BBOX.lat_min)
    fig = FigStd.figure(width=:wide, aspect=1.0 / map_wh * 1.1)
    ax = Axis(fig[1, 1];
              title="Fused-DEM source provenance (last wins)",
              xlabel="Longitude (°E)", ylabel="Latitude (°N)",
              aspect=DataAspect())
    image!(ax, (lons[1], lons[end]), (lats[1], lats[end]),
           permutedims(rgb); interpolate=false)
    xlims!(ax, BBOX.lon_min, BBOX.lon_max)
    ylims!(ax, BBOX.lat_min, BBOX.lat_max)

    # Legend
    elems = [PolyElement(color=colors[i]) for i in eachindex(labels)]
    Legend(fig[2, 1], elems, labels; orientation=:horizontal,
           framevisible=false, labelsize=9, patchsize=(16, 12), nbanks=2)
    rowsize!(fig.layout, 2, Fixed(60))

    mkpath(joinpath(OUTDIR, "pdf"))
    mkpath(joinpath(OUTDIR, "png"))
    FigStd.saveboth(joinpath(OUTDIR, "nile_source_provenance"), fig)
    @info "Saved" path=joinpath(OUTDIR, "png", "nile_source_provenance.png")
end

main()
