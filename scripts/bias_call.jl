using Distributed

addprocs(6)#, exeflags="--project=$(Base.active_project())")
@everywhere using DeltaDTM
@everywhere using DelimitedFiles
@everywhere using CSV
@everywhere using DataFrames
@everywhere using ProgressMeter
@everywhere using GeoDataFrames

folder = "data/bias3"
tiles = readdir(folder, join=true)
filter!(isdir, tiles)
@info length(tiles)

ofn(tilename) = !isfile(joinpath(@__DIR__, "bias", basename(tilename) * ".txt"))
filter!(ofn, tiles)
@info length(tiles)

outputfolder = "data/biasv1/"

@showprogress 1 "Computing..." pmap(tiles) do tilename
    tile = joinpath(folder, tilename)
    fn = joinpath(outputfolder, basename(tilename) * ".gpkg")
    if !isfile(fn)
        try
            vnt = DeltaDTM.process_tile(tile; divideby=2)
            df = DataFrame(vnt)
            GeoDataFrames.write(fn, df; geom_columns=(:geom,))
        catch e
            if e isa InterruptException
                throw(e)
            else
                @error tilename, e
            end
        end
    end
end
