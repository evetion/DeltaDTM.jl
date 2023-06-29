using DelimitedFiles
using GeoDataFrames

tiles = readdir(joinpath(@__DIR__, "bias"), join=true)
filter!(isfile, tiles)
filter!(x -> endswith(x, ".txt"), tiles)

@info length(tiles)

vals = Float64[]
vals2 = Float64[]
geoms = GeoDataFrames.AG.IGeometry[]

for tile in tiles
    corner = splitext(basename(tile))[1]
    lats, lons = split(corner, "_")[[4, 6]]
    lat =
        lats[1] == 'N' ? parse(Float64, lats[2:end]) :
        -parse(Float64, lats[2:end])
    lon =
        lons[1] == 'E' ? parse(Float64, lons[2:end]) :
        -parse(Float64, lons[2:end])
    val, val2 = readdlm(tile)
    push!(vals, val)
    push!(vals2, val2)
    point = GeoDataFrames.AG.createpoint([lon + 0.5, lat + 0.5])
    push!(geoms, point)
end

df = DataFrame(id=basename.(tiles), bias=vals, n=vals2, geom=geoms)
GeoDataFrames.write(joinpath(@__DIR__, "biasv1.gpkg"), df)
