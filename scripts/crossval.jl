using GeoParquet
using GeoDataFrames
using Statistics
using ProgressMeter
using CSV
using DataFrames
using Printf

withinrange(x, y) = sum(-y .<= x .<= y) / length(x)
mae(x) = sum(abs.(x)) / length(x)
rmse(x) = sqrt(sum(x .^ 2) / length(x))
mad(x) = median(abs.(x .- median(x)))


output_folder = "DeltaDTM/data/deltadtm/v3"
files = readdir(output_folder, join=true)
filter!(endswith(".pq"), files)


@enum landcover begin
    None = 0
    Trees = 10
    Shrubland = 20
    Grassland = 30
    Cropland = 40
    Urban = 50
    Bare = 60
    Snow = 70
    Water = 80
    Wetland = 90
    Mangroves = 95
    Moss = 100
end

allvalues = Float64[]
allcover = UInt8[]
rows = NamedTuple[]
@showprogress for tile in files
    corner = splitext(basename(tile))[1]
    lats, lons = split(corner, "_")[[4, 6]]
    lat =
        lats[1] == 'N' ? parse(Float64, lats[2:end]) :
        -parse(Float64, lats[2:end])
    lon =
        lons[1] == 'E' ? parse(Float64, lons[2:end]) :
        -parse(Float64, lons[2:end])
    point = GeoDataFrames.AG.createpoint([lon + 0.5, lat + 0.5])
    df = GeoParquet.read(tile)
    df = subset(df, :height => h -> isfinite.(h))
    df = subset(df, :deltadtm => h -> isfinite.(h))
    df = subset(df, :height => d -> -50 .<= d .<= 10)
    df = subset(df, :deltadtm => d -> -1000 .<= d .<= 1000)
    df = subset(df, :cover => c -> c .!== 0)
    df = subset(df, :cover => c -> c .!== 80)
    df.diff = df.deltadtm - df.height
    df = subset(df, :diff => d -> isfinite.(d))
    nrow(df) == 0 && continue
    values = df.diff
    filter!(isfinite, values)
    append!(allvalues, values)
    append!(allcover, df.cover)
    row = (;
        id=basename(tile),
        geom=point,
        n=length(values),
        mean=mean(values),
        mae=mae(values),
        mad=mad(values),
        rmse=rmse(values),
        b1m=withinrange(values, 1),
        b2m=withinrange(values, 2),
        b5m=withinrange(values, 5))
    push!(rows, row)
end
idf = DataFrame(rows)
GeoDataFrames.write("crossval_tiles.gpkg", idf, geom_columns=(:geom,))

# filter!(isfinite, allvalues)
rows = NamedTuple[]
@showprogress for cover in Integer.(instances(landcover))
    cover == 0 && continue
    cover == 80 && continue
    m = allcover .== cover
    v = allvalues[m]
    isempty(v) && continue
    row = (;
        cover=string(landcover(cover)),
        n=length(v),
        mean=mean(v),
        mae=mae(v),
        mad=mad(v),
        rmse=rmse(v),
        b1m=withinrange(v, 1),
        b2m=withinrange(v, 2),
        b5m=withinrange(v, 5)
    )
    push!(rows, row)
end
allvalues_f = allvalues[(allcover.!=0 .& allcover.!=80)]  # ignore evil stuff
row = (;
    cover="Overall",
    n=length(allvalues_f),
    mean=mean(allvalues_f),
    mae=mae(allvalues_f),
    mad=mad(allvalues_f),
    rmse=rmse(allvalues_f),
    b1m=withinrange(allvalues_f, 1),
    b2m=withinrange(allvalues_f, 2),
    b5m=withinrange(allvalues_f, 5)
)
push!(rows, row)

allvalues_fs = allvalues[(allcover.!=0 .& allcover.!=70 .& allcover.!=80)]  # ignore evil stuff
row = (;
    cover="Overall*",
    n=length(allvalues_fs),
    mean=mean(allvalues_fs),
    mae=mae(allvalues_fs),
    mad=mad(allvalues_fs),
    rmse=rmse(allvalues_fs),
    b1m=withinrange(allvalues_fs, 1),
    b2m=withinrange(allvalues_fs, 2),
    b5m=withinrange(allvalues_fs, 5)
)

push!(rows, row)
ndf = DataFrame(rows)


function quantize(x, v)
    v isa AbstractFloat || return v
    x >= 7 && return round(Int, v)
    @sprintf "%.2f" v
end

ndf.b1m .*= 100
ndf.b2m .*= 100
ndf.b5m .*= 100

CSV.write("crossval_overall2.csv", ndf, transform=quantize)
