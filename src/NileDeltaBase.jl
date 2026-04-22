"""
    NileDeltaBase

Access layer for a priority-stacked DEM covering the Nile Delta, Egypt.

Current stack (low → high priority, last wins):
  GMRT bathymetry/topography (multibeam + SRTM, ~100 m offshore)   [ocean only]
  Copernicus GLO-30 (30 m DSM, TanDEM-X-derived)
  FABDEM v1.2 (30 m bare-earth, ML veg/building removal from GLO-30)
  DeltaDTM v1.1 (30 m coastal bare-earth, Pronk et al. 2024, ICESat-2/GEDI-corrected)

Future layers:
  TanDEM-X 12 m (DLR science proposal, pending access)

Vertical datum: orthometric heights relative to EGM2008 (source-native for all
30 m layers). GMRT bathymetry is MSL-referenced (~equivalent to EGM2008 in the
Mediterranean to <1 m). No datum conversion applied — the Delta is <20 m,
vertical differences between EGM2008 and WGS84 ellipsoid across the region are
~10 m absolute but near-constant relative to InSAR baselines.

# Setup
    using NileDeltaBase
    NileDeltaBase.set_root!("/data4/EGY/NileDeltaBase")   # or ENV["NILEDELTABASE_ROOT"]

# Quick start
    h = sample_elevation(:fused_30m, 31.25, 30.05)       # Cairo, ~22 m
    data, lons, lats = extract_region(:fused_30m;
                                      lon_range=(30.5, 32.0),
                                      lat_range=(30.5, 31.5),
                                      res=0.00027778)     # ~30 m
"""
module NileDeltaBase

export BBOX, BBOX_EXT, HOTSPOTS, FRAMES, LAYERS,
       vrt_path, display_path, results_path,
       sample_elevation, sample_elevation_batch,
       sample_elevation_cached, warm_cache!, clear_cache!,
       load_display, load_basemap,
       classify_lecz, in_hotspot, in_frame, best_layer, in_bbox,
       subset, subset_hotspot, extract_region, extract_hotspot, extract_frame,
       resolution_source_map,
       dem_for_frame, sample_dem_at_points

using Statistics

include("Basemap.jl")
using .Basemap
export Basemap, BasemapInputs, compose_basemap, shaded_relief_base

# ── Configurable project root ────────────────────────────────────────

const _ROOT = Ref{String}("")

function _root()
    isempty(_ROOT[]) && _init_root!()
    return _ROOT[]
end

function _init_root!()
    r = get(ENV, "NILEDELTABASE_ROOT", "")
    if !isempty(r) && isdir(r)
        _ROOT[] = r
        return
    end
    for candidate in [
        "/data4/EGY/NileDeltaBase",
        joinpath(homedir(), "NileDeltaBase"),
    ]
        if isdir(joinpath(candidate, "data/vrt"))
            _ROOT[] = candidate
            return
        end
    end
    error("NileDeltaBase project root not found. Set ENV[\"NILEDELTABASE_ROOT\"] or call NileDeltaBase.set_root!(path)")
end

"""
    NileDeltaBase.set_root!(path)

Set the project root directory (must contain `data/vrt/`).
"""
function set_root!(path::AbstractString)
    isdir(path) || error("Not a directory: $path")
    isdir(joinpath(path, "data/vrt")) ||
        error("Missing data/vrt/ in $path")
    _ROOT[] = String(path)
    return path
end

_vrtdir()   = joinpath(_root(), "data/vrt")
_rawdir()   = joinpath(_root(), "data/raw")
_dispdir()  = joinpath(_root(), "data/derived/display")
_hydrodir() = joinpath(_root(), "data/derived/hydro")
_resdir()   = joinpath(_root(), "results")
_logdir()   = joinpath(_root(), "logs")

# ── Constants ────────────────────────────────────────────────────────

# Nile Delta only — Cairo (south) → Mediterranean coast (north),
# Alexandria / West Nubaria (west) → Port Said / Suez Canal mouth (east).
# ~310 km × 180 km ≈ 56,000 km².
const BBOX = (lon_min=29.6, lat_min=30.0, lon_max=32.4, lat_max=31.65)
# Extended display bbox — Dabaa → Port Said, shared with S06b / S09.
const BBOX_EXT = (lon_min=28.0, lat_min=30.0, lon_max=32.5, lat_max=32.5)
const DRES_DISPLAY = 0.01   # ~1 km display grid

# Sentinel-1 InSAR frames currently processed
const FRAMES = (
    nubaria_asc_131A = (lon=(29.0, 31.2), lat=(29.7, 31.5), track="131A"),
    delta_asc_058A   = (lon=(30.5, 32.7), lat=(29.7, 31.9), track="058A"),
    nubaria_dsc_065D = (lon=(29.0, 31.2), lat=(29.7, 31.5), track="065D"),
    delta_dsc_167D   = (lon=(30.5, 32.7), lat=(29.7, 31.9), track="167D"),
)

# Named regions of interest (subsidence hotspots, major cities, lagoons)
const HOTSPOTS = (
    # Nile Delta
    cairo      = (lon=(30.9, 31.5), lat=(29.9, 30.2), desc="Greater Cairo"),
    alexandria = (lon=(29.7, 30.1), lat=(31.1, 31.3), desc="Alexandria"),
    rosetta    = (lon=(30.3, 30.5), lat=(31.4, 31.5), desc="Rosetta promontory"),
    damietta   = (lon=(31.7, 31.9), lat=(31.4, 31.5), desc="Damietta promontory"),
    port_said  = (lon=(32.2, 32.4), lat=(31.2, 31.4), desc="Port Said"),
    manzala    = (lon=(31.6, 32.3), lat=(31.1, 31.4), desc="Lake Manzala"),
    burullus   = (lon=(30.5, 31.3), lat=(31.3, 31.6), desc="Lake Burullus"),
    nubaria    = (lon=(29.5, 30.5), lat=(30.3, 30.9), desc="West Nubaria reclamation"),
    # Eastern edge (Suez Canal mouth only — Sinai is out of the bbox)
    ismailia   = (lon=(32.2, 32.4), lat=(30.5, 30.7), desc="Ismailia"),
)

const LAYERS = (
    # Priority-fused VRTs
    fused_12m  = "nile_dem_fused_12m.vrt",   # + TanDEM-X 12 m (future)
    fused_30m  = "nile_dem_fused_30m.vrt",   # GMRT ← GLO30 ← FABDEM ← DeltaDTM (last wins)
    # Per-source VRTs
    gmrt       = "gmrt_nile.vrt",
    cop_glo30  = "cop_glo30_nile.vrt",
    fabdem     = "fabdem_nile.vrt",
    deltadtm   = "deltadtm_nile.vrt",
    tandemx12  = "tandemx12_nile.vrt",
)

# ── Path accessors ───────────────────────────────────────────────────

"""
    vrt_path(layer::Symbol) → String

Full path to a named VRT layer.
"""
vrt_path(layer::Symbol) = joinpath(_vrtdir(), getfield(LAYERS, layer))
vrt_path(name::String)  = joinpath(_vrtdir(), name)

"""
    display_path(product::Symbol) → String

Path to a display-resolution GeoTIFF.
"""
function display_path(product::Symbol)
    d = Dict(
        :elevation  => "nile_elevation_0p01deg.tif",
        :slope      => "nile_slope_0p01deg.tif",
        :dist_coast => "nile_dist_to_coast_0p01deg.tif",
    )
    haskey(d, product) || error("Unknown display product: $product. Options: $(keys(d))")
    return joinpath(_dispdir(), d[product])
end

results_path(name::String) = joinpath(_resdir(), name)

# ── Spatial queries ──────────────────────────────────────────────────

"""
    in_bbox(lon, lat) → Bool
"""
in_bbox(lon, lat) = BBOX.lon_min <= lon <= BBOX.lon_max && BBOX.lat_min <= lat <= BBOX.lat_max

"""
    in_hotspot(lon, lat) → Union{Symbol, Nothing}
"""
function in_hotspot(lon::Real, lat::Real)
    for (name, hs) in pairs(HOTSPOTS)
        if hs.lon[1] <= lon <= hs.lon[2] && hs.lat[1] <= lat <= hs.lat[2]
            return name
        end
    end
    return nothing
end

"""
    in_frame(lon, lat) → Vector{Symbol}

Returns all InSAR frame names that contain the given point.
"""
function in_frame(lon::Real, lat::Real)
    matches = Symbol[]
    for (name, fr) in pairs(FRAMES)
        if fr.lon[1] <= lon <= fr.lon[2] && fr.lat[1] <= lat <= fr.lat[2]
            push!(matches, name)
        end
    end
    return matches
end

"""
    best_layer(lon, lat) → Symbol

Highest-resolution fused layer available. Currently `:fused_30m` everywhere;
future `:fused_12m` once TanDEM-X is ingested.
"""
function best_layer(lon::Real, lat::Real)
    # Future: check TanDEM-X footprint → :fused_12m
    return :fused_30m
end

# ── Elevation sampling via gdallocationinfo ──────────────────────────

"""
    sample_elevation_batch(layer::Symbol, coords) → Vector{Union{Float64,Nothing}}

Sample a VRT at multiple `(lon, lat)` points via a single `gdallocationinfo`
subprocess. Returns `nothing` for nodata/out-of-coverage.
"""
function sample_elevation_batch(layer::Symbol, coords::Vector{Tuple{Float64,Float64}})
    return _sample_vrt_batch(vrt_path(layer), coords)
end

function sample_elevation_batch(vrt_file::String, coords::Vector{Tuple{Float64,Float64}})
    return _sample_vrt_batch(vrt_file, coords)
end

function _sample_vrt_batch(vrt, coords)
    isfile(vrt) || return fill(nothing, length(coords))
    tmp_in  = tempname() * ".txt"
    tmp_out = tempname() * ".out"
    open(tmp_in, "w") do io
        for (lon, lat) in coords
            println(io, "$lon $lat")
        end
    end
    try
        run(pipeline(`gdallocationinfo -wgs84 -valonly $vrt`;
                     stdin=tmp_in, stdout=tmp_out, stderr=devnull))
    catch; end
    out = isfile(tmp_out) ? read(tmp_out, String) : ""
    rm(tmp_in; force=true); rm(tmp_out; force=true)
    lines = split(out, '\n', keepempty=true)
    vals = Vector{Union{Nothing,Float64}}(undef, length(coords))
    for i in eachindex(coords)
        s = i <= length(lines) ? strip(lines[i]) : ""
        if isempty(s)
            vals[i] = nothing
        else
            v = tryparse(Float64, s)
            if v === nothing || abs(v) > 1e30 || v == -99999 || v == -999999
                vals[i] = nothing
            else
                vals[i] = v
            end
        end
    end
    return vals
end

"""
    sample_elevation(layer::Symbol, lon, lat) → Union{Float64, Nothing}
"""
function sample_elevation(layer::Symbol, lon::Real, lat::Real)
    v = sample_elevation_batch(layer, [(Float64(lon), Float64(lat))])
    return v[1]
end

# ── Display raster loading ───────────────────────────────────────────

"""
    load_display(product::Symbol) → (data, lons, lats)

Load a display-resolution raster. Row 1 = lat_min (south-up). NaN for nodata.
"""
function load_display(product::Symbol)
    tif = display_path(product)
    isfile(tif) || error("Display raster not found: $tif\nRun derived-product script first.")
    return _read_tif(tif, DRES_DISPLAY)
end

function _read_tif(tif, dres)
    info = read(`gdalinfo -json $tif`, String)
    m_size = match(r"\"size\":\s*\[\s*(\d+),\s*(\d+)\s*\]", info)
    nx = parse(Int, m_size[1]); ny = parse(Int, m_size[2])
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
    arr[abs.(arr) .> 1e4] .= NaN32
    lons = collect(range(BBOX.lon_min + dres/2, BBOX.lon_max - dres/2; length=nx))
    lats = collect(range(BBOX.lat_min + dres/2, BBOX.lat_max - dres/2; length=ny))
    return arr, lons, lats
end

# ── LECZ classification ──────────────────────────────────────────────

"""
    classify_lecz(h) → Symbol

LECZ zones relative to EGM2008 orthometric height. The Delta is almost
entirely <10 m; most of it is <2 m and at risk from sea-level rise.
"""
function classify_lecz(h::Real)
    isnan(h) && return :nodata
    h < 2  && return :below_2m
    h < 5  && return :below_5m
    h < 10 && return :below_10m
    h < 50 && return :below_50m
    return :above_50m
end

# ── Subsetting (from already-loaded arrays) ─────────────────────────

"""
    subset(data, lons, lats; lon_range, lat_range) → (sub, lo, la)
"""
function subset(data::Matrix, lons::Vector, lats::Vector;
                lon_range::Tuple{Real,Real}, lat_range::Tuple{Real,Real})
    ic = findall(x -> lon_range[1] <= x <= lon_range[2], lons)
    ir = findall(x -> lat_range[1] <= x <= lat_range[2], lats)
    isempty(ic) && error("No longitude pixels in range $lon_range")
    isempty(ir) && error("No latitude pixels in range $lat_range")
    return data[ir, ic], lons[ic], lats[ir]
end

"""
    subset_hotspot(name::Symbol, data, lons, lats) → (sub, lo, la)
"""
function subset_hotspot(name::Symbol, data::Matrix, lons::Vector, lats::Vector)
    hs = getfield(HOTSPOTS, name)
    return subset(data, lons, lats; lon_range=hs.lon, lat_range=hs.lat)
end

# ── High-resolution extraction (from VRT via gdalwarp) ──────────────

"""
    extract_region(layer; lon_range, lat_range, res=0.00027778) → (data, lons, lats)

Extract a subset directly from a VRT using `gdalwarp`.
Default `res=1/3600 ≈ 30 m` matches the source resolution of the 30 m layers.
Returns `(Matrix{Float32}, Vector{Float64}, Vector{Float64})`, south-up, NaN for nodata.
"""
function extract_region(layer::Symbol;
                        lon_range::Tuple{Real,Real},
                        lat_range::Tuple{Real,Real},
                        res::Real=1/3600)
    return _extract_vrt(vrt_path(layer), lon_range, lat_range, Float64(res))
end

"""
    extract_hotspot(name, layer; res=1/3600) → (data, lons, lats)
"""
function extract_hotspot(name::Symbol, layer::Symbol; res::Real=1/3600)
    hs = getfield(HOTSPOTS, name)
    return extract_region(layer; lon_range=hs.lon, lat_range=hs.lat, res=res)
end

"""
    extract_frame(name, layer; res=1/3600) → (data, lons, lats)

Extract the bbox of an InSAR frame (see `FRAMES`) at the given resolution.
"""
function extract_frame(name::Symbol, layer::Symbol; res::Real=1/3600)
    fr = getfield(FRAMES, name)
    return extract_region(layer; lon_range=fr.lon, lat_range=fr.lat, res=res)
end

function _extract_vrt(vrt, lon_range, lat_range, res)
    isfile(vrt) || error("VRT not found: $vrt")
    bin = tempname() * ".bin"
    hdr = replace(bin, ".bin" => ".hdr")
    xmin, xmax = Float64.(lon_range)
    ymin, ymax = Float64.(lat_range)
    try
        run(pipeline(`gdalwarp -overwrite
            -t_srs EPSG:4326
            -te $xmin $ymin $xmax $ymax
            -tr $res $res
            -r bilinear -wo NUM_THREADS=ALL_CPUS
            -of ENVI -ot Float32 -dstnodata -9999
            $vrt $bin`; stderr=devnull))
    catch e
        rm(bin; force=true); rm(hdr; force=true)
        rm(bin * ".aux.xml"; force=true)
        rethrow(e)
    end
    info = read(`gdalinfo -json $bin`, String)
    m_size = match(r"\"size\":\s*\[\s*(\d+),\s*(\d+)\s*\]", info)
    nx = parse(Int, m_size[1]); ny = parse(Int, m_size[2])
    data = Array{Float32}(undef, nx, ny)
    read!(bin, data)
    rm(bin; force=true); rm(hdr; force=true)
    rm(bin * ".aux.xml"; force=true)
    arr = Matrix{Float32}(undef, ny, nx)
    for r in 1:ny
        arr[r, :] = data[:, ny - r + 1]
    end
    arr[arr .== -9999f0] .= NaN32
    arr[abs.(arr) .> 1e4] .= NaN32
    lons = collect(range(xmin + res/2, xmax - res/2; length=nx))
    lats = collect(range(ymin + res/2, ymax - res/2; length=ny))
    return arr, lons, lats
end

# ── Resolution source map ───────────────────────────────────────────

"""
    resolution_source_map(; res=0.01) → (src_map, lons, lats)

Pixel-wise source-provenance map for the fused product.
Codes: 0=nodata, 1=GMRT (bathy), 2=COP-GLO30 (30 m DSM),
       3=FABDEM (30 m bare-earth), 4=DeltaDTM (30 m coastal DTM),
       5=TanDEM-X (12 m, future).

Coverage is inferred from the presence of each layer's VRT/filelist.
"""
function resolution_source_map(; res::Real=0.01)
    nx = round(Int, (BBOX.lon_max - BBOX.lon_min) / res)
    ny = round(Int, (BBOX.lat_max - BBOX.lat_min) / res)
    lons = collect(range(BBOX.lon_min + res/2, BBOX.lon_max - res/2; length=nx))
    lats = collect(range(BBOX.lat_min + res/2, BBOX.lat_max - res/2; length=ny))

    src_map = zeros(Int8, ny, nx)

    # Priority: later layers overwrite earlier ones (matches fusion order)
    layer_codes = [(:gmrt, Int8(1)), (:cop_glo30, Int8(2)), (:fabdem, Int8(3)),
                   (:deltadtm, Int8(4)), (:tandemx12, Int8(5))]

    for (layer, code) in layer_codes
        vrt = vrt_path(layer)
        isfile(vrt) || continue
        bin = tempname() * ".bin"
        hdr = replace(bin, ".bin" => ".hdr")
        xmin, xmax = Float64(BBOX.lon_min), Float64(BBOX.lon_max)
        ymin, ymax = Float64(BBOX.lat_min), Float64(BBOX.lat_max)
        try
            run(pipeline(`gdalwarp -overwrite
                -t_srs EPSG:4326
                -te $xmin $ymin $xmax $ymax
                -tr $res $res
                -r near -wo NUM_THREADS=ALL_CPUS
                -of ENVI -ot Float32 -dstnodata -9999
                $vrt $bin`; stderr=devnull))
            info = read(`gdalinfo -json $bin`, String)
            m_s = match(r"\"size\":\s*\[\s*(\d+),\s*(\d+)\s*\]", info)
            nnx = parse(Int, m_s[1]); nny = parse(Int, m_s[2])
            raw = Array{Float32}(undef, nnx, nny)
            read!(bin, raw)
            for c in 1:min(nnx, nx), r in 1:min(nny, ny)
                v = raw[c, nny - r + 1]
                if v != -9999f0 && v != -999999f0 && abs(v) < 1e4
                    src_map[r, c] = code
                end
            end
        catch; end
        rm(bin; force=true); rm(hdr; force=true)
        rm(bin * ".aux.xml"; force=true)
    end

    return src_map, lons, lats
end

# ── In-memory cache for fast repeated sampling ──────────────────────
#
# `sample_elevation_batch` spawns one `gdallocationinfo` subprocess per
# call (~100 ms per invocation on a warm filesystem).  For bulk InSAR PS
# sampling (10k+ points) that adds up.  The cache loads a full display
# raster into RAM once and does in-memory bilinear interpolation.

mutable struct _RasterCache
    layer::Symbol
    res::Float64
    bbox::NamedTuple
    data::Matrix{Float32}
    lons::Vector{Float64}
    lats::Vector{Float64}
end

const _CACHE = Ref{Union{Nothing, _RasterCache}}(nothing)

"""
    warm_cache!(; layer=:fused_30m, res=0.0025, bbox=BBOX_EXT) → Matrix

Load a layer into RAM for fast repeated point sampling. Default grid is
~275 m at the extended bbox (~1800×680 Float32, ~5 MB).  Returns the
loaded matrix.
"""
function warm_cache!(; layer::Symbol=:fused_30m, res::Real=0.0025,
                      bbox::NamedTuple=BBOX_EXT)
    arr, lons, lats = _extract_vrt(vrt_path(layer),
                                    (bbox.lon_min, bbox.lon_max),
                                    (bbox.lat_min, bbox.lat_max),
                                    Float64(res))
    _CACHE[] = _RasterCache(layer, Float64(res), bbox, arr, lons, lats)
    return arr
end

"""
    clear_cache!()
"""
clear_cache!() = (_CACHE[] = nothing; nothing)

@inline function _bilinear(c::_RasterCache, lon::Float64, lat::Float64)
    (lon < c.bbox.lon_min || lon > c.bbox.lon_max ||
     lat < c.bbox.lat_min || lat > c.bbox.lat_max) && return NaN32
    fx = (lon - c.lons[1]) / c.res
    fy = (lat - c.lats[1]) / c.res
    j0 = floor(Int, fx) + 1
    i0 = floor(Int, fy) + 1
    j1 = min(j0 + 1, size(c.data, 2))
    i1 = min(i0 + 1, size(c.data, 1))
    j0 = max(j0, 1); i0 = max(i0, 1)
    ax = Float32(fx - (j0 - 1)); ay = Float32(fy - (i0 - 1))
    v00 = c.data[i0, j0]; v01 = c.data[i0, j1]
    v10 = c.data[i1, j0]; v11 = c.data[i1, j1]
    (isnan(v00) || isnan(v01) || isnan(v10) || isnan(v11)) && return NaN32
    return (1-ax)*(1-ay)*v00 + ax*(1-ay)*v01 + (1-ax)*ay*v10 + ax*ay*v11
end

"""
    sample_elevation_cached(lons, lats) → Vector{Float32}

Bilinear-interpolate the warmed cache at N points. ~3000× faster than
`sample_elevation_batch` for 10k+ points (no subprocess, no disk I/O).
Call `warm_cache!()` first.
"""
function sample_elevation_cached(lons::AbstractVector, lats::AbstractVector)
    length(lons) == length(lats) || throw(DimensionMismatch("lons/lats length"))
    c = _CACHE[]
    c === nothing && error("Cache not warmed — call warm_cache!() first")
    out = Vector{Float32}(undef, length(lons))
    @inbounds for i in eachindex(lons)
        out[i] = _bilinear(c, Float64(lons[i]), Float64(lats[i]))
    end
    return out
end

sample_elevation_cached(lon::Real, lat::Real) =
    sample_elevation_cached([Float64(lon)], [Float64(lat)])[1]

# ── SARMotion / InSAR integration helpers ──────────────────────────

"""
    dem_for_frame(frame::Symbol, layer::Symbol=:fused_30m;
                  res=1/3600, out_path=nothing) → String

Clip the fused DEM to an InSAR frame's footprint and write a GeoTIFF
ready to be used as `EQA.dem_par` input (SARMotion / LiCSBAS).  Defaults
to `results/<frame>_<layer>.tif` if `out_path` is nothing.
"""
function dem_for_frame(frame::Symbol, layer::Symbol=:fused_30m;
                       res::Real=1/3600,
                       out_path::Union{String, Nothing}=nothing)
    fr = getfield(FRAMES, frame)
    out = out_path !== nothing ? out_path :
          joinpath(_resdir(), "dem", "$(frame)_$(layer).tif")
    mkpath(dirname(out))
    run(`gdalwarp -overwrite
         -t_srs EPSG:4326
         -te $(Float64(fr.lon[1])) $(Float64(fr.lat[1]))
             $(Float64(fr.lon[2])) $(Float64(fr.lat[2]))
         -tr $(Float64(res)) $(Float64(res))
         -r bilinear -wo NUM_THREADS=ALL_CPUS
         -of GTiff -ot Float32 -dstnodata -9999
         -co TILED=YES -co COMPRESS=DEFLATE
         $(vrt_path(layer)) $out`)
    return out
end

"""
    sample_dem_at_points(lons, lats; layer=:fused_30m, use_cache=true) → Vector{Float32}

Sample the DEM at N point locations (e.g., InSAR PS or SBAS pixels).
If `use_cache=true` and the cache is warmed, uses in-memory bilinear;
otherwise falls back to `gdallocationinfo` batch sampling.

Returns a `Vector{Float32}` with `NaN32` for out-of-coverage points.
"""
function sample_dem_at_points(lons::AbstractVector, lats::AbstractVector;
                              layer::Symbol=:fused_30m,
                              use_cache::Bool=true)
    if use_cache
        c = _CACHE[]
        (c !== nothing && c.layer === layer) ||
            warm_cache!(; layer=layer)
        return sample_elevation_cached(lons, lats)
    else
        coords = [(Float64(lons[i]), Float64(lats[i])) for i in eachindex(lons)]
        vs = sample_elevation_batch(layer, coords)
        return Float32[v === nothing ? NaN32 : Float32(v) for v in vs]
    end
end

# ── Basemap convenience wrapper ─────────────────────────────────────

"""
    load_basemap(; bbox=BBOX_EXT, dres=0.003, with_bathy=true) → BasemapInputs

Load all display rasters (elevation/hillshade/landmask/WorldCover/bathymetry)
from `data/derived/display/` and return a `BasemapInputs` struct. Used by
S09 and external projects that want the Delta basemap as a backdrop.
"""
function load_basemap(; bbox::NamedTuple=BBOX_EXT, dres::Real=0.003,
                       with_bathy::Bool=true)
    return Basemap.load_basemap_inputs(_dispdir(), bbox, dres;
                                        with_bathy=with_bathy)
end

end # module
