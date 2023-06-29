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

colors = LittleDict(
    None => "#ffffff",
    Trees => "#006400",
    Shrubland => "#ffbb22",
    Grassland => "#ffff4c",
    Cropland => "#f096ff",
    Urban => "#fa0000",
    Bare => "#b4b4b4",
    Snow => "#f0f0f0",
    Water => "#0064c8",
    Wetland => "#0096a0",
    Mangroves => "#00cf75",
    Moss => "#fae6a0",
)

landcovercolors(landcovers) = getindex.(Ref(colors), sort(unique(landcovers)))
center(nt::NamedTuple{(:min_x, :min_y, :max_x, :max_y),NTuple{4,T}}) where {T<:Real} = ((nt.min_x + nt.max_x) / 2, (nt.min_y + nt.max_y) / 2)

function biasmask(sample, cover)
    covermask =
        in.(
            landcover.(cover),
            Ref((Shrubland, Grassland, Cropland, Bare, Moss, Snow)),
        )
    nodatamask = sample .!= -9999 .& .!isnan.(sample) .& isfinite.(sample)
    mask = covermask .& nodatamask
    return mask
end

function plot_bias(dff, b, tile, sub=1)
    groupedhist(
        dff.diff,
        group=dff.covername,
        bar_position=:stack,
        xlim=(-2, 2),
        markerstrokewidth=0,
        linewidth=0,
        xlabel="CopernicusDEM - ICESat-2 [m]",
        ylabel="Count",
        xticks=-2:0.5:2,
        # yaxis=false,
        dpi=300,
        color_palette=landcovercolors(dff.covername),
    )
    vline!([b], label="correction", color="black")
    StatsPlots.annotate!(
        b + 0.05,
        0,
        text(@sprintf("%.2f", b), :left, pointsize=8),
    )
    # StatsPlots.savefig(joinpath(@__DIR__, "../scripts/bias/plots_2", first(splitext(basename(tile))) * "_$sub.png"))
    closeall()
    return nothing
end

function stack(tile)
    ices = Vector{DataFrame}()
    for f in readlines(joinpath(tile, "bias_single.txt"))
        df = GeoParquet.read(f)
        push!(ices, df)
    end
    ice = reduce(vcat, ices)
end

function tiles(lon, lat, divideby=3)
    out = []
    step = 1 // divideby
    ((min_x=x, min_y=y, max_x=x + step, max_y=y + step) for x in lon:step:lon+1-step, y in lat:step:lat+1-step)
end

function process_tile(tile; divideby=3)
    lon, lat = parse_name(tile)
    df = stack(tile)

    mask = biasmask(df.sample, df.cover)
    dff = df[mask, :]

    dff.covername = (landcover.(dff.cover))
    dff.diff = dff.sample - dff.height  # copernicusdem - icesat2

    subset!(dff, :diff => x -> -30 .< x .< 30)

    out = NamedTuple[]
    for (i, subtile) in enumerate(tiles(Float64(lon), Float64(lat), divideby))
        sdf = in_bbox(dff, subtile)
        bias = if nrow(sdf) > 0
            den = kde(sdf.diff)
            mval = findmax(den.density)
            bias = den.x[mval[2]]::Float64
            plot_bias(sdf, bias, tile, i)
            bias
        else
            NaN
        end
        push!(out, (geom=center(subtile), bias=bias, n=nrow(sdf)))
    end
    out
end
