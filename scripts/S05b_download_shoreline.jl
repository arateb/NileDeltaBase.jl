#!/usr/bin/env julia
# Shoreline vectors for NileDeltaBase — GSHHG 2.3.7 (coastline + rivers + borders),
# clipped to the Nile Delta bbox.  Used by S08 as vector overlays and by S06 to
# rasterize a landmask (so ocean/lagoon pixels render as pale off-white instead
# of leaking through the elevation-based flood fill).
#
# Run: julia --project=/data/files/pkgs/NileDeltaBase scripts/S05b_download_shoreline.jl

using Pkg
Pkg.activate(dirname(@__DIR__))

using Downloads, Dates

const ROOT   = get(ENV, "NILEDELTABASE_ROOT", "/data4/EGY/NileDeltaBase")
const OUTDIR = joinpath(ROOT, "data/raw/shoreline")

# Narrow bbox (Delta only — S08 plot) and extended bbox (Dabaa → Suez — S09 plot).
# Two passes so both S08 and S09 can run without re-downloading GSHHG.
const BBOXES = [
    (suffix="nile", bbox=(lon_min=29.6, lat_min=30.0, lon_max=32.4, lat_max=31.65)),
    (suffix="ext",  bbox=(lon_min=28.0, lat_min=30.0, lon_max=32.5, lat_max=32.50)),
]

const SOURCES = [
    (label    = "gshhg",
     url      = "https://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip",
     zip_name = "gshhg-shp-2.3.7.zip"),
]

function download_zip(url, dst; timeout=1800)
    isfile(dst) && filesize(dst) > 100_000 && (@info "cached" dst; return)
    tmp = dst * ".part"
    @info "downloading" url
    Downloads.download(url, tmp; timeout=timeout)
    mv(tmp, dst; force=true)
end

function extract_zip(zip, extract_dir)
    mkpath(extract_dir)
    run(pipeline(`unzip -o -q $zip -d $extract_dir`, stderr=devnull))
end

function clip_shp_dir(src_dir, out_prefix, bbox, suffix)
    shps = String[]
    for (root, _, files) in walkdir(src_dir)
        for f in files
            endswith(lowercase(f), ".shp") && push!(shps, joinpath(root, f))
        end
    end
    out_files = String[]
    for shp in shps
        bn  = replace(basename(shp), r"\.shp$"i => "")
        out = joinpath(OUTDIR, "$(out_prefix)_$(bn)_$(suffix).shp")
        isfile(out) && filesize(out) > 0 && (push!(out_files, out); continue)
        try
            run(`ogr2ogr -overwrite
                 -clipsrc $(bbox.lon_min) $(bbox.lat_min) $(bbox.lon_max) $(bbox.lat_max)
                 -skipfailures
                 $out $shp`)
            isfile(out) && filesize(out) > 0 && push!(out_files, out)
        catch e
            @warn "clip failed" shp err=sprint(showerror, e)
        end
    end
    return out_files
end

function process_source(s)
    ext_dir = joinpath(OUTDIR, "_ext_$(s.label)")
    mkpath(ext_dir)
    zip = joinpath(OUTDIR, "_zips", s.zip_name); mkpath(dirname(zip))
    try
        download_zip(s.url, zip)
        extract_zip(zip, ext_dir)
        all_out = String[]
        for b in BBOXES
            out_files = clip_shp_dir(ext_dir, s.label, b.bbox, b.suffix)
            append!(all_out, out_files)
        end
        total_mb  = isempty(all_out) ? 0 : sum(filesize.(all_out)) ÷ (1024^2)
        rm(zip; force=true)
        rm(ext_dir; recursive=true, force=true)
        return (label=s.label, status="ok", n_shp=length(all_out),
                size_mb=total_mb, url=s.url,
                files=join(basename.(all_out), ";"))
    catch e
        @warn "failed" s.label err=sprint(showerror, e)
        return (label=s.label, status="fail: " * sprint(showerror, e),
                n_shp=0, size_mb=0, url=s.url, files="")
    end
end

function write_manifest(results)
    path = joinpath(OUTDIR, "shoreline_manifest.csv")
    ts = string(now())
    open(path, "w") do io
        println(io, "label,status,n_shp,size_mb,url,files,downloaded_at")
        for r in results
            status = replace(r.status, ',' => ';')
            println(io, "$(r.label),$status,$(r.n_shp),$(r.size_mb),$(r.url),$(r.files),$ts")
        end
    end
    @info "manifest" path
end

function main()
    mkpath(OUTDIR)
    results = [process_source(s) for s in SOURCES]
    write_manifest(results)
    zdir = joinpath(OUTDIR, "_zips")
    isdir(zdir) && isempty(readdir(zdir)) && rm(zdir; recursive=true)
    n_ok = count(r -> r.status == "ok", results)
    @info "shoreline done" n_ok n_fail=length(results)-n_ok
end

main()
