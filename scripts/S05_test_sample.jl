#!/usr/bin/env julia
# Smoke test: sample elevation at known landmarks across the bbox.

using Pkg
Pkg.activate(dirname(@__DIR__))

push!(LOAD_PATH, joinpath(dirname(@__DIR__), "src"))
using NileDeltaBase
using Statistics

NileDeltaBase.set_root!(get(ENV, "NILEDELTABASE_ROOT", "/data4/EGY/NileDeltaBase"))

# Known ~elevations (EGM2008 orthometric, m)
points = [
    ("Cairo (Tahrir)",     31.236, 30.046, 23.0),
    ("Alexandria centre",  29.920, 31.200,  4.0),
    ("Rosetta",            30.417, 31.405,  2.0),
    ("Damietta",           31.812, 31.418,  1.0),
    ("Port Said",          32.284, 31.265,  3.0),
    ("Ismailia",           32.272, 30.588,  7.0),
    ("Lake Manzala (ctr)", 31.950, 31.250,  0.0),
    ("Lake Burullus (ctr)",31.000, 31.450,  0.0),
    ("Nubaria town",       30.176, 30.676, 25.0),
    ("Lake Bardawil (ctr)",33.100, 31.150,  0.0),
    ("El Arish",           33.802, 31.131, 30.0),
    ("Rafah",              34.240, 31.290, 70.0),
    ("Taba (village)",     34.895, 29.495,  5.0),
    ("Suez (N end canal)", 32.550, 29.967,  5.0),
]

println("\n=== NileDeltaBase sample_elevation test ===\n")
println(rpad("Location", 22), " ", rpad("lon", 9), " ", rpad("lat", 9), " ",
        rpad("expected", 10), " ", rpad("fused_30m", 11), " layer")
println("-"^95)

for layer in (:fused_30m, :cop_glo30, :fabdem, :deltadtm)
    vrt = NileDeltaBase.vrt_path(layer)
    if !isfile(vrt)
        @info "skipping $layer — VRT not built yet" vrt
        continue
    end
    println("\n--- layer: $layer ---")
    for (name, lon, lat, expected) in points
        h = sample_elevation(layer, lon, lat)
        hstr = h === nothing ? "nodata" : string(round(h, digits=2), " m")
        println(rpad(name, 22), " ", rpad(string(lon), 9), " ",
                rpad(string(lat), 9), " ", rpad("~$(expected) m", 10), " ",
                rpad(hstr, 11))
    end
end

println("\n=== resolution_source_map ===")
src, lons, lats = resolution_source_map(; res=0.05)
pixels = sum(src .> 0)
total  = length(src)
println("Total pixels (0.05°): $total, covered: $pixels ($(round(100*pixels/total, digits=1))%)")
for code in 1:5
    n = sum(src .== code)
    label = ("COP-GLO30", "FABDEM", "DeltaDTM", "TanDEM-X 12m", "EarthDEM 2m")[code]
    n > 0 && println("  code $code ($label): $n pixels")
end
