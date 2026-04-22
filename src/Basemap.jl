"""
    NileDeltaBase.Basemap

Reusable shaded-relief + WorldCover + bathymetry compositor. Produces an
`RGBAf` matrix suitable for `image!` in CairoMakie. Used by S09 (landcover
basemap) and exposable to external projects (SARMotion figures, texas_vlm,
etc.) that want the same Delta basemap as a backdrop.
"""
module Basemap

using Statistics
using ColorTypes
const RGBAf = RGBA{Float32}
const RGBf  = RGB{Float32}

export shaded_relief_base, compose_basemap, wc_tint, alpha_blend,
       BasemapInputs, load_basemap_inputs, read_envi_raster

# ── WorldCover tint palette (class codes → low-sat overlays) ────────────
const WC_CROPLAND = RGBAf(0.55, 0.70, 0.30, 0.55)
const WC_BUILT    = RGBAf(0.70, 0.35, 0.25, 0.70)
const WC_TREE     = RGBAf(0.25, 0.55, 0.25, 0.50)
const WC_SHRUB    = RGBAf(0.50, 0.60, 0.30, 0.40)
const WC_GRASS    = RGBAf(0.65, 0.75, 0.40, 0.40)
const WC_WATER    = RGBAf(0.30, 0.50, 0.75, 0.80)
const WC_WETLAND  = RGBAf(0.40, 0.60, 0.65, 0.60)

@inline function wc_tint(code::UInt8)
    code == 0x28 && return WC_CROPLAND   # 40
    code == 0x32 && return WC_BUILT      # 50
    code == 0x0A && return WC_TREE       # 10
    code == 0x14 && return WC_SHRUB      # 20
    code == 0x1E && return WC_GRASS      # 30
    code == 0x50 && return WC_WATER      # 80
    code == 0x5A && return WC_WETLAND    # 90
    return RGBAf(0, 0, 0, 0)
end

@inline function alpha_blend(base::RGBAf, top::RGBAf)
    a = top.alpha
    a == 0 && return base
    r = base.r * (1 - a) + top.r * a
    g = base.g * (1 - a) + top.g * a
    b = base.b * (1 - a) + top.b * a
    return RGBAf(r, g, b, 1.0)
end

# ── Bathymetry colour ramp (shallow → deep, pale cyan → deep indigo) ───
@inline function bathy_color(depth::Real; dmax::Real=2000)
    # depth positive downward; clamp to [0, dmax]
    d = clamp(Float32(depth), 0f0, Float32(dmax))
    t = sqrt(d / Float32(dmax))             # compress dynamic range
    # shallow: RGB(0.78, 0.88, 0.95), deep: RGB(0.06, 0.12, 0.35)
    r = Float32(0.78 - 0.72 * t)
    g = Float32(0.88 - 0.76 * t)
    b = Float32(0.95 - 0.60 * t)
    return RGBAf(r, g, b, 1.0)
end

# ── Shaded relief base (hypsometric + hillshade) ────────────────────────
"""
    shaded_relief_base(elev, hs, landmask; bathy=nothing, vmin=0, vmax=200,
                       dmax=2000, ocean=RGBAf(0.78, 0.88, 0.95, 1))

Return an `RGBAf` matrix sized like `elev`.

  * Land pixels (`landmask==1`): gray hypsometric ramp scaled to `[vmin, vmax]`,
    multiplied by `0.35 + 0.65*hs` for shading.
  * Ocean pixels (`landmask==0`):
      – if `bathy` is provided and < 0 at the pixel, render with `bathy_color`,
        multiplied by the same hillshade modulation for texture.
      – otherwise render as `ocean` (pale off-white).
"""
function shaded_relief_base(elev::AbstractMatrix, hs::AbstractMatrix,
                            landmask::AbstractMatrix;
                            bathy::Union{Nothing, AbstractMatrix}=nothing,
                            vmin::Real=0f0, vmax::Real=200f0,
                            dmax::Real=2000f0,
                            ocean::RGBAf=RGBAf(0.78, 0.88, 0.95, 1.0))
    ny, nx = size(elev)
    img = Matrix{RGBAf}(undef, ny, nx)
    @inbounds for j in 1:nx, i in 1:ny
        v = elev[i, j]
        h = 0.35f0 + 0.65f0 * Float32(hs[i, j])
        if landmask[i, j] == 0x01 && !isnan(v)
            t = clamp((Float32(v) - Float32(vmin)) /
                      (Float32(vmax) - Float32(vmin)), 0f0, 1f0)
            g_base = Float32(1.0 - 0.82 * t)
            if v < 0
                g_base = 0.96f0
            end
            img[i, j] = RGBAf(g_base * h, g_base * h, g_base * h, 1.0)
        else
            # Ocean / outside landmask
            if bathy !== nothing
                b = bathy[i, j]
                if !isnan(b) && b < 0
                    c = bathy_color(-b; dmax=dmax)
                    img[i, j] = RGBAf(c.r * h, c.g * h, c.b * h, 1.0)
                    continue
                end
            end
            img[i, j] = ocean
        end
    end
    return img
end

"""
    compose_basemap(elev, hs, landmask, wc; bathy=nothing, kwargs...)

Shaded relief base, then overlay WorldCover tints on land. Water stays as
bathymetry/ocean. Returns `Matrix{RGBAf}`.
"""
function compose_basemap(elev::AbstractMatrix, hs::AbstractMatrix,
                         landmask::AbstractMatrix, wc::AbstractMatrix;
                         bathy::Union{Nothing, AbstractMatrix}=nothing,
                         kwargs...)
    base = shaded_relief_base(elev, hs, landmask; bathy=bathy, kwargs...)
    ny, nx = size(base)
    img = similar(base)
    @inbounds for j in 1:nx, i in 1:ny
        if landmask[i, j] == 0x01
            img[i, j] = alpha_blend(base[i, j], wc_tint(wc[i, j]))
        else
            img[i, j] = base[i, j]
        end
    end
    return img
end

# ── Raster I/O via GDAL CLI (no GDAL.jl dep) ────────────────────────────

"""
    read_envi_raster(tif, T) → (Matrix{T}, nx, ny)

GDAL-translate `tif` to raw ENVI, then `read!` into memory. Returns a
south-up matrix of element type `T` (`Float32` or `UInt8`).
"""
function read_envi_raster(tif::AbstractString, ::Type{T}) where T
    info = read(`gdalinfo -json $tif`, String)
    m = match(r"\"size\":\s*\[\s*(\d+),\s*(\d+)\s*\]", info)
    nx, ny = parse(Int, m[1]), parse(Int, m[2])
    ot = T === UInt8 ? "Byte" : "Float32"
    bin = tempname() * ".bin"
    run(pipeline(`gdal_translate -of ENVI -ot $ot $tif $bin`; stderr=devnull))
    raw = Array{T}(undef, nx, ny)
    read!(bin, raw)
    rm(bin; force=true)
    rm(bin * ".aux.xml"; force=true)
    rm(replace(bin, ".bin" => ".hdr"); force=true)
    arr = Matrix{T}(undef, ny, nx)
    @inbounds for r in 1:ny
        arr[r, :] = raw[:, ny - r + 1]
    end
    return arr, nx, ny
end

"""Container for a loaded basemap at a chosen bbox/grid."""
struct BasemapInputs
    elev::Matrix{Float32}
    hs::Matrix{Float32}
    landmask::Matrix{UInt8}
    wc::Matrix{UInt8}
    bathy::Union{Nothing, Matrix{Float32}}
    lons::Vector{Float64}
    lats::Vector{Float64}
end

"""
    load_basemap_inputs(dispdir, bbox, dres; with_bathy=true) → BasemapInputs

Load elevation / hillshade / landmask / WorldCover / bathymetry from display
TIFs at `dispdir`, naming convention `nile_<kind>_display_ext.tif`.
"""
function load_basemap_inputs(dispdir::AbstractString, bbox::NamedTuple,
                             dres::Real; with_bathy::Bool=true)
    dem_tif  = joinpath(dispdir, "nile_elevation_display_ext.tif")
    hls_tif  = joinpath(dispdir, "nile_hillshade_display_ext.tif")
    lmsk_tif = joinpath(dispdir, "nile_landmask_display_ext.tif")
    wc_tif   = joinpath(dispdir, "nile_worldcover_display_ext.tif")
    bth_tif  = joinpath(dispdir, "nile_bathymetry_display_ext.tif")

    elev_raw, nx, ny = read_envi_raster(dem_tif, Float32)
    info = read(`gdalinfo -json $dem_tif`, String)
    m_nd = match(r"\"noDataValue\":\s*([-\d.eE+]+)", info)
    nd   = m_nd !== nothing ? parse(Float32, m_nd[1]) : -9999f0
    elev_raw[elev_raw .== nd]        .= NaN32
    elev_raw[abs.(elev_raw) .> 1f4]  .= NaN32

    hs_u8, _, _ = read_envi_raster(hls_tif, UInt8)
    hs_raw = Float32.(hs_u8) ./ 255f0
    valid = filter(x -> x > 0.01f0, hs_raw[:])
    lo, hi = isempty(valid) ? (0f0, 1f0) : (Float32(quantile(valid, 0.01)),
                                             Float32(quantile(valid, 0.99)))
    span = max(hi - lo, 0.01f0)
    hs = similar(hs_raw)
    @inbounds for i in eachindex(hs_raw)
        hs[i] = hs_raw[i] > 0.01f0 ? clamp((hs_raw[i] - lo) / span, 0f0, 1f0) : 0f0
    end

    lmsk, _, _ = read_envi_raster(lmsk_tif, UInt8)
    wc,   _, _ = read_envi_raster(wc_tif,   UInt8)

    bathy = nothing
    if with_bathy && isfile(bth_tif)
        b, _, _ = read_envi_raster(bth_tif, Float32)
        b[b .== -9999f0]        .= NaN32
        b[abs.(b) .> 1f5]       .= NaN32
        bathy = b
    end

    lons = collect(range(bbox.lon_min + dres/2, bbox.lon_max - dres/2; length=nx))
    lats = collect(range(bbox.lat_min + dres/2, bbox.lat_max - dres/2; length=ny))
    return BasemapInputs(elev_raw, hs, lmsk, wc, bathy, lons, lats)
end

end # module
