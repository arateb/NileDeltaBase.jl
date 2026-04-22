#!/usr/bin/env julia
# Inter-layer residual maps — validates datum consistency and surfaces
# systematic biases between COP-GLO30 (DSM, includes buildings/trees),
# FABDEM (DSM with veg/buildings ML-removed), and DeltaDTM (bare-earth).
#
# Residuals computed at the extended display grid (0.003°).  Plots:
#   diff_fabdem_minus_glo30  — expect mostly ≤0 (veg/buildings subtracted)
#   diff_deltadtm_minus_fabdem — Delta tidal-correction bias + ICESat-2 tie
#   diff_deltadtm_minus_glo30 — total Delta-DTM offset from raw DSM
#
# Each figure shows a map + residual histogram side-by-side.

using Pkg
Pkg.activate(dirname(@__DIR__))

using Statistics, Printf
using CairoMakie

include("/data/files/pkgs/figstd/src/FigStd.jl")
using .FigStd
FigStd.apply!()

push!(LOAD_PATH, joinpath(dirname(@__DIR__), "src"))
using NileDeltaBase
const NDB = NileDeltaBase
NDB.set_root!(get(ENV, "NILEDELTABASE_ROOT", "/data4/EGY/NileDeltaBase"))

const ROOT     = NDB._root()
const OUTDIR   = joinpath(ROOT, "results", "figures")
const TMPDIR   = joinpath(ROOT, "data", "derived", "residuals")
mkpath(TMPDIR)

# Extended bbox & resolution — share with S06b / S09
const BBOX = (lon_min=28.0, lat_min=30.0, lon_max=32.5, lat_max=31.7)
const DRES = 0.003

function warp_to_grid(vrt, out_tif)
    (isfile(out_tif) && stat(out_tif).size > 0) && return out_tif
    run(`gdalwarp -overwrite
         -t_srs EPSG:4326
         -te $(BBOX.lon_min) $(BBOX.lat_min) $(BBOX.lon_max) $(BBOX.lat_max)
         -tr $DRES $DRES
         -r bilinear -wo NUM_THREADS=ALL_CPUS
         -of GTiff -ot Float32 -dstnodata -9999
         -co TILED=YES -co COMPRESS=DEFLATE
         $vrt $out_tif`)
    return out_tif
end

function read_float(tif)
    info = read(`gdalinfo -json $tif`, String)
    m = match(r"\"size\":\s*\[\s*(\d+),\s*(\d+)\s*\]", info)
    nx, ny = parse(Int, m[1]), parse(Int, m[2])
    bin = tempname() * ".bin"
    run(pipeline(`gdal_translate -of ENVI -ot Float32 $tif $bin`; stderr=devnull))
    raw = Array{Float32}(undef, nx, ny)
    read!(bin, raw)
    rm(bin; force=true); rm(bin * ".aux.xml"; force=true)
    rm(replace(bin, ".bin" => ".hdr"); force=true)
    arr = Matrix{Float32}(undef, ny, nx)
    for r in 1:ny
        arr[r, :] = raw[:, ny - r + 1]
    end
    arr[arr .== -9999f0]   .= NaN32
    arr[abs.(arr) .> 1f4]  .= NaN32
    return arr
end

function load_aligned(layer::Symbol)
    vrt = layer === :deltadtm ? joinpath(ROOT, "data/vrt/deltadtm_nile_masked.vrt") :
                                 NDB.vrt_path(layer)
    out = joinpath(TMPDIR, "$(layer)_aligned_ext.tif")
    warp_to_grid(vrt, out)
    return read_float(out)
end

function plot_residual(diff, title_str, fname; vlim=30f0)
    lons = collect(range(BBOX.lon_min + DRES/2, BBOX.lon_max - DRES/2; step=DRES))
    lats = collect(range(BBOX.lat_min + DRES/2, BBOX.lat_max - DRES/2; step=DRES))
    nx, ny = length(lons), length(lats)
    # Handle size mismatch from gdalwarp rounding
    dr = diff[1:min(ny, size(diff,1)), 1:min(nx, size(diff,2))]
    lons = lons[1:size(dr,2)]; lats = lats[1:size(dr,1)]

    vals = filter(!isnan, dr[:])
    μ = mean(vals); σ = std(vals); med = median(vals)
    p05, p95 = quantile(vals, 0.05), quantile(vals, 0.95)

    map_wh = (BBOX.lon_max - BBOX.lon_min) / (BBOX.lat_max - BBOX.lat_min)
    fig = FigStd.figure(width=:wide, aspect=1.0 / map_wh * 0.7)

    ax = Axis(fig[1, 1]; title=title_str,
              xlabel="Longitude (°E)", ylabel="Latitude (°N)",
              aspect=DataAspect())
    hm = heatmap!(ax, lons, lats, permutedims(dr);
                  colormap=FigStd.berlin_gray(), colorrange=(-vlim, vlim),
                  lowclip=:black, highclip=:red, nan_color=RGBAf(1, 1, 1, 0))
    xlims!(ax, BBOX.lon_min, BBOX.lon_max)
    ylims!(ax, BBOX.lat_min, BBOX.lat_max)
    Colorbar(fig[1, 2], hm; label="residual (m)", width=12)

    ax2 = Axis(fig[1, 3]; title="Residual distribution",
               xlabel="Δ (m)", ylabel="density",
               limits=(-vlim, vlim, nothing, nothing))
    bins = range(-vlim, vlim; length=101)
    hist!(ax2, vals; bins=bins, color=RGBAf(0.3, 0.45, 0.7, 0.75),
          normalization=:pdf, strokecolor=:black, strokewidth=0.3)
    vlines!(ax2, [μ], color=:red, linewidth=1.5, label="μ = $(round(μ, digits=2))")
    vlines!(ax2, [med], color=:orange, linewidth=1.5, linestyle=:dash,
            label="med = $(round(med, digits=2))")
    axislegend(ax2; position=:rt, labelsize=7, framevisible=false)

    colsize!(fig.layout, 1, Relative(0.55))
    colsize!(fig.layout, 3, Relative(0.30))

    @info "residual stats" title=title_str μ med σ p05 p95 n=length(vals)
    mkpath(joinpath(OUTDIR, "pdf")); mkpath(joinpath(OUTDIR, "png"))
    FigStd.saveboth(joinpath(OUTDIR, fname), fig)
    return (; μ, med, σ, p05, p95, n=length(vals))
end

function main()
    @info "Warping and loading aligned rasters..."
    glo = load_aligned(:cop_glo30)
    fab = load_aligned(:fabdem)
    del = load_aligned(:deltadtm)

    stats = Dict{Symbol, NamedTuple}()

    @info "FABDEM - GLO30"
    stats[:fab_minus_glo] = plot_residual(
        fab .- glo,
        "FABDEM − COP-GLO30 (bare-earth vs DSM)",
        "nile_residual_fabdem_minus_glo30"; vlim=20f0)

    @info "DeltaDTM - FABDEM"
    stats[:del_minus_fab] = plot_residual(
        del .- fab,
        "DeltaDTM − FABDEM (delta-tuned vs generic bare-earth)",
        "nile_residual_deltadtm_minus_fabdem"; vlim=10f0)

    @info "DeltaDTM - GLO30"
    stats[:del_minus_glo] = plot_residual(
        del .- glo,
        "DeltaDTM − COP-GLO30 (end-to-end offset)",
        "nile_residual_deltadtm_minus_glo30"; vlim=25f0)

    # Summary table → stdout
    println("\n=== Residual summary (metres) ===")
    @printf("%-22s %8s %8s %8s %8s %8s %10s\n", "diff", "μ", "median", "σ", "p05", "p95", "n")
    for (k, s) in stats
        @printf("%-22s %8.2f %8.2f %8.2f %8.2f %8.2f %10d\n",
                string(k), s.μ, s.med, s.σ, s.p05, s.p95, s.n)
    end
end

main()
