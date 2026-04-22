using Test, NileDeltaBase
const NDB = NileDeltaBase

NDB.set_root!(get(ENV, "NILEDELTABASE_ROOT", "/data4/EGY/NileDeltaBase"))

@testset "NileDeltaBase — geometry" begin
    @test NDB.in_bbox(31.25, 30.05)            # Cairo
    @test !NDB.in_bbox(50.0, 10.0)             # Arabia
    @test NDB.in_hotspot(31.25, 30.05) === :cairo
    @test NDB.in_hotspot(29.92, 31.20) === :alexandria
    @test NDB.in_hotspot(40.0, 10.0) === nothing

    frames = NDB.in_frame(31.25, 30.05)
    @test !isempty(frames)
end

@testset "NileDeltaBase — layers & paths" begin
    @test haskey(NDB.LAYERS, :fused_30m)
    @test haskey(NDB.LAYERS, :gmrt)
    @test haskey(NDB.LAYERS, :cop_glo30)
    @test !haskey(NDB.LAYERS, :earthdem2)        # ensure cleanup
    @test !haskey(NDB.LAYERS, :fused_2m)

    @test endswith(NDB.vrt_path(:fused_30m), ".vrt")
    @test isfile(NDB.vrt_path(:fused_30m))
    @test isfile(NDB.vrt_path(:cop_glo30))

    @test NDB.best_layer(31.25, 30.05) === :fused_30m
end

@testset "NileDeltaBase — classification" begin
    @test NDB.classify_lecz(1.0)   === :below_2m
    @test NDB.classify_lecz(3.0)   === :below_5m
    @test NDB.classify_lecz(7.0)   === :below_10m
    @test NDB.classify_lecz(30.0)  === :below_50m
    @test NDB.classify_lecz(100.0) === :above_50m
    @test NDB.classify_lecz(NaN)   === :nodata
end

@testset "NileDeltaBase — sample_elevation (VRT)" begin
    # Each row: (lon, lat, expected (m), tolerance (m))
    # Delta bare-earth values are 2–10 m typically; uplands up to 200 m.
    # Tolerances are generous — we just want to catch outright breakage
    # (datum flip, nodata leakage, or DeltaDTM 30 m saturation regression).
    points = [
        ("Cairo",      31.236, 30.046,  20.0, 10.0),
        ("Alexandria", 29.920, 31.200,   6.0,  8.0),
        ("Damietta",   31.812, 31.418,   3.0,  8.0),
        ("Manzala",    31.950, 31.250,   0.0,  3.0),
        ("Rafah",      34.240, 31.290,  50.0, 40.0),
        ("Suez",       32.550, 29.967,   5.0, 10.0),
    ]
    for (name, lon, lat, expected, tol) in points
        h = sample_elevation(:fused_30m, lon, lat)
        @test h !== nothing
        @test isfinite(h)
        @test -20.0 < h < 500.0
        @test abs(h - expected) ≤ tol
    end

    # Hotspot-scale sampling
    for hs_name in (:cairo, :alexandria, :damietta, :manzala)
        hs = getfield(NDB.HOTSPOTS, hs_name)
        h = sample_elevation(:fused_30m, (hs.lon[1] + hs.lon[2])/2,
                                         (hs.lat[1] + hs.lat[2])/2)
        @test h !== nothing
        @test isfinite(h)
    end
end

@testset "NileDeltaBase — DeltaDTM 30 m saturation regression" begin
    # Upland points that previously returned exactly 30.0 due to DeltaDTM
    # saturation. After the S04 pre-mask fix the fused VRT should fall
    # through to FABDEM/GLO30 and return the true upland elevation.
    h1 = sample_elevation(:fused_30m, 29.7, 30.2)   # NW desert
    h2 = sample_elevation(:fused_30m, 32.0, 30.2)   # E desert near Ismailia
    @test h1 !== nothing && h1 > 30.1              # not stuck at 30 m
    @test h2 !== nothing && h2 > 30.1
end

@testset "NileDeltaBase — cache" begin
    NDB.warm_cache!(; layer=:fused_30m, res=0.005, bbox=NDB.BBOX_EXT)
    lons = [31.236, 29.920, 31.812, 31.950]
    lats = [30.046, 31.200, 31.418, 31.250]
    vs = sample_elevation_cached(lons, lats)
    @test length(vs) == 4
    @test all(isfinite, vs)
    @test -10 < vs[1] < 200       # Cairo
    @test -5  < vs[2] < 30        # Alexandria
    @test -5  < vs[4] < 5         # Manzala

    # Round-trip helper
    vs2 = sample_dem_at_points(lons, lats; use_cache=true)
    @test all(isfinite, vs2)
    @test vs ≈ vs2

    NDB.clear_cache!()
end

@testset "NileDeltaBase — resolution_source_map" begin
    src, lons, lats = resolution_source_map(; res=0.05)
    @test size(src) == (length(lats), length(lons))
    @test sum(src .> 0) > 0
    # In the narrow Delta bbox DeltaDTM (code 4) dominates.
    @test sum(src .== 4) > sum(src .== 2)
end
