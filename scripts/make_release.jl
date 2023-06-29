using GeoDataFrames
using ProgressMeter
using DataFrames
using GeoArrays
using DataAPI

# Coastal lowland tiles
fn = joinpath(@__DIR__, "data/coastal_lowland_tiles_v4.gpkg")
df = GeoDataFrames.read(fn)

input_folder = "DeltaDEM/data/deltadem/v7"
output_folder = "DeltaDEM/data/deltadem/v1.0"

md = [
    ("TIFFTAG_DOCUMENTNAME", "DeltaDEM v1.0"),
    ("TIFFTAG_IMAGEDESCRIPTION", "A global coastal digital terrain model, based on CopernicusDEM, ESA WorldCover, ICESat-2 and GEDI data. For more information, see X et al. (2023) DeltaDEM: A global coastal digital terrain model."),
    ("TIFFTAG_SOFTWARE", "DeltaDEM.jl"),
    ("TIFFTAG_ARTIST", "X"),
    ("TIFFTAG_COPYRIGHT", "CC-BY-SA 4.0"),
]
limit = 10

"""Copernicus_DSM_10_N00_00_E106_00_DEM.tif becomes DeltaDEM_v1_0_N00E099.tif"""
function todd(copdemtile)
    corner = splitext(basename(copdemtile))[1]
    lats, lons = split(corner, "_")[[4, 6]]
    "DeltaDEM_v1_0_" * lats * lons * ".tif"
end

subset!(df, :tile => tile -> isfile.(joinpath.(input_folder, tile)))

@showprogress for i in 1:nrow(df)
    folder = joinpath("data")
    tilename = basename(df.dem[i])
    tilef = splitext(basename(df.tile[i]))[1]
    tilefilename = basename(df.tile[i])
    mtile = joinpath(folder, "copernicus/copernicus_orig/masks", tilefilename)
    itile = joinpath(input_folder, tilefilename)

    otile = joinpath(output_folder, todd(tilefilename))
    omtile = joinpath(output_folder, "masks", todd(tilefilename))

    # continue on existing output
    isfile(otile) && isfile(omtile) && continue
    # continue on non-existing input
    isfile(itile) || continue

    ga = GeoArrays.coalesce(GeoArrays.read(itile), -9999)
    gam = GeoArrays.read(mtile)

    m = ga .> limit
    ga[m] .= limit
    gam[m] .= 255

    for (k, v) in md
        DataAPI.metadata!(ga, k, v)
    end

    epsg!(ga, 9518)
    epsg!(gam, 9518)

    GeoArrays.write(otile, ga; nodata=-9999, shortname="COG", options=Dict("COMPRESS" => "ZSTD", "PREDICTOR" => "3"))
    GeoArrays.write(omtile, gam; nodata=-9999, shortname="COG", options=Dict("COMPRESS" => "ZSTD", "PREDICTOR" => "2"))
end

# Make deltadem_tiles.gpkg
df.tile = todd.(df.tile)
subset!(df, :tile => tile -> isfile.(joinpath.(output_folder, tile)))
rename!(df, Dict(:loc => :copernicus, :x => :longitude, :y => :latitude))
select!(df, Not(Cols(:dem, :project, :artefacts, :path, :ice, :id, :val)))
df.zipfile = df.copernicus .|> x -> split(x, "/")[2] .* ".zip"
GeoDataFrames.write("deltadem_tiles.gpkg", df)

# Zip tiles per continent
for sdf in groupby(df, :zipfile)
    @info "zipping $(sdf.zipfile[1])"
    files = joinpath.(output_folder, sdf.tile)[1:2]
    cmd = `zip -0 -j $(joinpath(output_folder, sdf.zipfile[1])) $files`
    @info cmd
end

# Zip for all masks files (within folder)
cd(output_folder)
mask_files = filter!(x -> endswith(x, ".tif"), readdir("masks", join=true))
run(`zip -0 $(joinpath(output_folder, "mask_tiles.zip")) $(mask_files)`)
