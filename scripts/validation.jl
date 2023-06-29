using GeoDataFrames
using DataFrames
using GeoArrays
using Statistics
using CSV
using ProgressMeter
using Printf

function compare(a::AbstractString, b::AbstractString, c, threshold=10)
    ga = GeoArrays.read(a)
    ga2 = GeoArrays.read(b)
    ga3 = GeoArrays.read(c)
    A = Float32.(coalesce(ga, -9999))  # can be Int
    A[A.==-9999] .= Inf
    B = coalesce(ga2, Inf)
    C = coalesce(ga3, 0)

    m = isfinite.(A) .& isfinite.(B)
    lm = (B .<= threshold) .& m
    @info sum(lm), sum(m), length(A)
    diffl = A[lm] - B[lm]
    return diffl, C[lm]
end

withinrange(x, y) = sum(-y .<= x .<= y) / length(x)
mae(x) = sum(abs.(x)) / length(x)
rmse(x) = sqrt(sum(x .^ 2) / length(x))
mad(x) = length(x) == 0 ? NaN : median(abs.(x .- median(x)))

fn = joinpath(@__DIR__, "data/coastal_lowland_tiles_v4.gpkg")
df = GeoDataFrames.read(fn)
subset!(df, :val => val -> val .!= "");  # only validation tiles

dfolder = joinpath("data")
esa = joinpath(dfolder, "esa-worldcover-2021/copernicus2")
valfolder = joinpath(dfolder, "validation/")

hlimit = 10
llimit = 0
scalebar_size = 0.1
dems = Dict(
    "NASADEM" => "data/nasadem",
    "CopernicusDEM" => "data/copernicus/copernicus2",
    "MERIT" => "data/merit",
    "CoastalDEM" => "data/coastaldem",
    "FABDEM" => "data/FABDEM",
    "DeltaDEM" => "DeltaDEM/data/demo/flevodem",
    "DeltaDEM" => "DeltaDEM/data/deltadem/v1",
)

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

data = DataFrame[]
for (dem, folder) in dems
    values, landcovers = Float32[], UInt8[]
    @showprogress map(1:nrow(df)) do i
        tilename = first(splitext(basename(df.dem[i])))
        demfn = joinpath(folder, tilename)
        coverfn = joinpath(esa, tilename)
        d, c = compare(demfn, df.val[i], coverfn)
        append!(values, d)
        append!(landcovers, c)
    end
    rows = NamedTuple[]
    for cover in Integer.(instances(landcover))
        m = landcovers .== cover
        v = values[m]
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
    length(values) == 0 && continue
    values = values[(landcovers.!=0).&(landcovers.!=80)]
    row = (;
        cover="Overall",
        n=length(values),
        mean=mean(values),
        mae=mae(values),
        mad=mad(values),
        rmse=rmse(values),
        b1m=withinrange(values, 1),
        b2m=withinrange(values, 2),
        b5m=withinrange(values, 5)
    )
    push!(rows, row)
    ndf = DataFrame(rows)
    ndf.dem .= dem
    push!(data, ndf)
end

vdf = reduce(vcat, data)


data_areas = DataFrame[]
for project in unique(df.project)
    for (dem, folder) in dems
        values, landcovers = Float32[], UInt8[]
        sdf = df[df.project.==project, :]
        @showprogress map(1:nrow(sdf)) do i
            tilename = first(splitext(basename(sdf.dem[i])))
            demfn = joinpath(folder, tilename)
            coverfn = joinpath(esa, tilename)
            d, c = compare(demfn, sdf.val[i], coverfn)
            append!(values, d)
            append!(landcovers, c)
        end
        rows = NamedTuple[]
        for cover in Integer.(instances(landcover))
            m = landcovers .== cover
            v = values[m]
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
        length(values) == 0 && continue
        values = values[(landcovers.!=0).&(landcovers.!=80)]
        row = (;
            cover="Overall",
            n=length(values),
            mean=mean(values),
            mae=mae(values),
            mad=mad(values),
            rmse=rmse(values),
            b1m=withinrange(values, 1),
            b2m=withinrange(values, 2),
            b5m=withinrange(values, 5)
        )
        push!(rows, row)
        ndf = DataFrame(rows)
        ndf.project .= project
        push!(data_areas, ndf)
    end
end
vadf = reduce(vcat, data_areas)


# order for table
coverf = Dict(
    "Trees" => 90,
    "Shrubland" => 50,
    "Grassland" => 30,
    "Cropland" => 40,
    "Urban" => 51,
    "Bare" => 1,
    "Snow" => 70,
    "Water" => 80,
    "Wetland" => 2,
    "Mangroves" => 95,
    "Moss" => 99,
    "Overall" => 100,
)
demf = Dict(
    "NASADEM" => 1,
    "CopernicusDEM" => 2,
    "MERIT" => 3,
    "CoastalDEM" => 4,
    "FABDEM" => 5,
    "DeltaDEM" => 6,
)

function quantize(x, v)
    v isa AbstractFloat || return v
    x >= 7 && return round(Int, v)
    @sprintf "%.2f" v
end

subset!(vdf, :cover => cover -> cover .!= "None")
subset!(vdf, :cover => cover -> cover .!= "Water")
sort!(vdf, [order(:cover, by=x -> coverf[x]), order(:dem, by=x -> demf[x])])

vdf.b1m .*= 100
vdf.b2m .*= 100
vdf.b5m .*= 100

CSV.write("validation_overall.csv", vdf, transform=quantize)

subset!(vadf, :cover => cover -> cover .!= "None")
subset!(vadf, :cover => cover -> cover .!= "Water")
sort!(vadf, [order(:project), order(:cover, by=x -> coverf[x])])

vadf.b1m .*= 100
vadf.b2m .*= 100
vadf.b5m .*= 100

CSV.write("validation_per_area.csv", vadf, transform=quantize)
