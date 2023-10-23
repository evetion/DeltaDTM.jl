# using Revise
using Distributed
using GeoDataFrames
using ProgressMeter
using DataFrames
using ProgressMeter

# Coastal lowland tiles
fn = joinpath(@__DIR__, "data/coastal_lowland_tiles_v4.gpkg")
df = GeoDataFrames.read(fn)

output_folder = "DeltaDTM/data/deltadtm/v1"

# Only process tiles that are not yet processed
subset!(df, :tile => tile -> .!isfile.(joinpath.(output_folder, basename.(tile))))


# For multiprocessing from commandline
nproc = 5
p = Progress(nrow(df), dt=1.0)
for subdf in collect(Iterators.partition(df, nproc * 10))
    t = addprocs(nproc, exeflags="--heap-size-hint=8G -JDeltaDTM_all.so")
    @everywhere begin
        using DeltaDTM
        using ProgressMeter
        output_folder = "DeltaDTM/data/deltadtm/v1"
    end

    pmap(procs()) do _
        folder = joinpath(@__DIR__, "../../data")
        global tin = DeltaDTM.setup_bias(joinpath(folder, "biasv5_2.gpkg"))
        nothing
    end

    pmap(1:nrow(subdf)) do i
        folder = joinpath("data")

        tilename = basename(subdf.dem[i])
        tilef = splitext(basename(subdf.tile[i]))[1]
        tilefilename = basename(subdf.tile[i])
        tile = joinpath(folder, "copernicus/copernicus2", tilename)
        otile = joinpath(folder, "copernicus/copernicus_orig/elevation", tilename)
        etile = joinpath(folder, "copernicus/copernicus_orig/error", tilename)
        mtile = joinpath(folder, "copernicus/copernicus_orig/masks", tilename)
        mtiles = joinpath(folder, "copernicus/copernicus_orig/masks", tilefilename)
        ctile = joinpath(folder, "esa-worldcover-2021/copernicus2", tilename)
        icefn = joinpath(folder, "bias3/$(tilef)/bias.txt")
        gedifn = joinpath(folder, "biasg2/$(tilef)/bias.txt")

        ofn = joinpath(output_folder, tilefilename)
        isfile(ofn) && return nothing

        try
            DeltaDTM.deltadtm(tile, subdf.x[i], subdf.y[i], icefn, gedifn, mtile, etile, otile, ctile, tin, output_folder; cropsize=0.1, lowlimit=subdf.lowlimit[i], crossval=false, debug=false)
        catch e
            e isa InterruptException && throw(e)
            Base.showerror(stdout, e)
            @error "Something went wrong at $i: $e"
        finally
            GC.gc()
        end
    end
    next!(p, step=nrow(subdf))
    rmprocs(t)
    nothing
end
