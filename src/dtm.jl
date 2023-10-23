function setup_bias(bias_fn=joinpath(@__DIR__, "../data/biasv4.gpkg"))
    bias = GeoDataFrames.read(bias_fn)
    dropmissing!(bias)
    subset!(bias, :n => x -> x .> 250)
    subset!(bias, :bias => x -> -2.5 .< x .< 2.5)
    locs = GeoInterface.coordinates.(bias.geom)
    bias.x = Float64.(first.(locs))
    bias.y = Float64.(last.(locs))
    interpolator(bias)
end

function parse_name(tile)
    corner = splitext(basename(tile))[1]
    lats, lons = split(corner, "_")[[4, 6]]
    l = parse(Int64, lats[2:end])
    lat = lats[1] == 'N' ? l : -l
    l = parse(Int64, lons[2:end])
    lon = lons[1] == 'E' ? l : -l
    return lon, lat
end

function interpolator(t)
    tt = hcat(t.x, t.y, t.bias)
    ttt = permutedims(tt)
    dt = DT()
    insert!(dt, ttt)
    return dt
end

"""
Create a SpaceLiDAR data `df` only DTM, for a given GeoArray `ga` size and extents, 
but downsampled by a factor `f` (default 30, so 30 m * 30 = 900 m)
"""
function multiresdtm(ga, df; maxdepth=5, n=3, reducer=median, r2=std)
    local g
    local pg
    for div in reverse(0:maxdepth)
        w, h = size(ga)[1:2] .÷ 2^div
        if (w == 1) || (h == 1)
            maxdepth -= 1
            continue
        end
        A = fill(0.0, w, h, 2)
        g = GeoArray(A, ga.f, ga.crs, Dict{String,Any}())
        GeoArrays.bbox!(g, GeoArrays.bbox(ga))
        GeoArrays.flipud!(g)

        d = Dictionary{CartesianIndex{2},Vector{Float32}}()

        for p in Tables.rows(df)
            I = indices(g, SVector{2}(p.longitude, p.latitude), GeoArrays.Center())
            if I in keys(d)
                push!(d[I], p.height)
            else
                insert!(d, I, [p.height])
            end
        end
        if div != maxdepth
            A = imresize(pg[:, :, 1], size(g)[1:2], method=Interpolations.BSpline(Interpolations.Linear()))
            B = imresize(pg[:, :, 2], size(g)[1:2], method=Interpolations.BSpline(Interpolations.Linear()))

            m = trues(size(A))
            for (k, v) in pairs(d)
                g[k, 1] = reducer(v)
                g[k, 2] = r2(v)
                m[k] = length(v) <= n
            end
            g[m, 1] .= A[m]
            g[m, 2] .= B[m]
        end
        pg = g
    end
    g
end


function downsample(A, f=24, reduce=any)
    TI = TileIterator(axes(A), (f, f))
    new = similar(A, size(TI))
    for I in eachindex(TI)
        @inbounds new[I] = reduce(@view A[TI[I]...])
    end
    return new
end

function apply_bias!(ga, dt)
    ui, uj = size(ga)[1:2]
    m = ismissing.(ga)
    c = collect(GeoArrays.coords(ga))
    for i = 1:ui, j = 1:uj
        m[i, j] && continue
        @inbounds ga[i, j, 1] -= interpolate_linear(dt, c[i, j][1], c[i, j][2])
    end
end


function burn_icesat(ga, df, watermask)
    lx, ly, _ = size(ga)
    mask = falses((lx, ly))
    ind = zeros(Int, (lx, ly))
    values = fill(Inf, size(ga)[1:2])
    differences = fill(Inf, size(ga)[1:2])
    for n = 1:nrow(df)
        isnan(df.height[n]) && continue
        i, j = indices(ga, (df.longitude[n], df.latitude[n])).I
        ((0 < i <= lx) && (0 < j <= ly)) && !ismissing(ga.A[i, j, 1]) ||
            continue
        watermask[i, j] && continue
        diff = df.height[n] - ga[i, j, 1]
        if diff < 0
            values[i, j] = df.height[n]
            ind[i, j] = n
            mask[i, j] = true
            differences[i, j] = diff
        end
    end
    return mask, values, differences, ind
end

function burn_gedi(ga, df, watermask)
    lx, ly, _ = size(ga)
    mask = falses((lx, ly))
    ind = zeros(Int, (lx, ly))
    values = fill(Inf, size(ga)[1:2])
    differences = fill(Inf, size(ga)[1:2])
    for n = 1:nrow(df)
        isnan(df.height[n]) && continue
        i, j = indices(ga, (df.longitude[n], df.latitude[n])).I
        ((0 < i <= lx) && (0 < j <= ly)) && !ismissing(ga.A[i, j, 1]) ||
            continue
        watermask[i, j] && continue
        diff = df.height[n] - ga[i, j, 1]
        if diff < 0
            values[i, j] = df.height[n]
            ind[i, j] = n
            mask[i, j] = true
            differences[i, j] = diff
        end
    end
    return mask, values, differences, ind
end

function sample_mask(df, ga)
    lx, ly, _ = size(ga)
    values = trues(nrow(df))
    for n = 1:nrow(df)
        i, j = indices(ga, (df.longitude[n], df.latitude[n])).I
        ((0 < i <= lx) && (0 < j <= ly)) && !ismissing(ga.A[i, j, 1]) ||
            continue
        values[n] = ga.A[i, j, 1]
    end
    return values
end


function sample!(df, ga; name=:csample)
    lx, ly, _ = size(ga)
    df[!, name] = Vector{eltype(ga)}(undef, nrow(df))
    for n = 1:nrow(df)
        i, j = indices(ga, (df.longitude[n], df.latitude[n])).I
        if ((0 < i <= lx) && (0 < j <= ly)) && !ismissing(ga.A[i, j, 1])
            df[!, name][n] = ga.A[i, j, 1]
        else
            df[!, name][n] = NaN
        end
    end
    nothing
end

function GeoArrays.crop(ga, x, y, buffer=0.1)
    lx, my = indices(ga, (x - buffer, y - buffer), GeoArrays.Center()).I
    mx, ly = indices(ga, (x + 1 + buffer, y + 1 + buffer), GeoArrays.Center()).I
    mmx, mmy = size(ga)[1:2]
    ix, iy = max(1, lx):min(mmx, mx), max(1, ly):min(mmy, my)
    length(ix) % 2 == 0 || (ix = ix[1:end-1])
    length(iy) % 2 == 0 || (iy = iy[1:end-1])
    return ga[ix, iy]
end


function tileddeviation(A; w=25, wo=w ÷ 2)
    myview = TiledView(A, (w, w), (wo, wo))
    function filter_pits(tile)
        t = filter(isfinite, skipmissing(tile))
        length(t) > 0 || return tile
        μ = median(t)
        σ = length(t) > 1 ? std(t) : 1
        @.(tile - μ) / σ
    end
    tiled_processing(myview, filter_pits, verbose=false).parent
end

function tileddiff(A; w=25, wo=w ÷ 2)
    myview = TiledView(A, (w, w), (wo, wo))
    function filter_pits(tile)
        t = filter(isfinite, skipmissing(tile))
        length(t) > 0 || return tile
        μ = median(t)
        σ = length(t) > 1 ? std(t) : 1
        @.(tile - μ)
    end
    tiled_processing(myview, filter_pits, verbose=false).parent
end



function deltadtm(tile, x, y, icefn, gedifn, mtile, etile, otile, ctile, tin, output_folder; cropsize::Union{Nothing,Float64}=0.05, lowlimit=-Inf, crossval=false, debug=false)
    output_fn = splitext(basename(tile))[1]
    ofn = joinpath(output_folder, output_fn)
    ofn = joinpath(output_folder, output_fn)

    ga = GeoArrays.read(tile)
    gao = GeoArrays.read(otile)
    gam = GeoArrays.read(mtile)
    ea = DeltaDTM.coalesce2d(gam, 255)
    any(ea .== 0) || return output_fn

    # Open error file
    isfile(etile) || return output_fn
    ega = GeoArrays.read(etile)

    # Open land cover
    lga = GeoArrays.read(ctile)

    # Crop to smaller buffer, complete degrees are a bit much...
    if !isnothing(cropsize)
        ga = crop(ga, x, y, cropsize)
        gam = crop(gam, x, y, cropsize)
        ega = crop(ega, x, y, cropsize)
        lga = crop(lga, x, y, cropsize)
        gao = crop(gao, x, y, cropsize)
        @assert size(ga) == size(gam) == size(ega) == size(lga) == size(gao)
    end

    ea = DeltaDTM.coalesce2d(ega, Inf)
    cover = DeltaDTM.coalesce2d(lga, 0)

    # Apply bias correction
    DeltaDTM.apply_bias!(ga, tin)

    # Clamp original roughness (TPI) for pattern burning later on
    O = DeltaDTM.coalesce2d(gao, 0)
    pattern = GeoArrayOps.TPI(O)
    A = DeltaDTM.coalesce2d(ga, Inf)

    # Keep all data higher than 100m
    high_mask = A .> 100
    nm = isfinite.(A)
    high_mask .&= nm

    # Std deviation from median to find low outliers in CopernicusDEM
    dev = tileddeviation(A; w=25)
    low_mask = (dev .< -2.0) .| (A .< lowlimit)
    A[low_mask] .= Inf  # so it doesn't influence filters
    orig = copy(A)

    ocean_watermask = coalesce2d(gam, 0) .== 1  # only ocean
    watermask = coalesce2d(gam, 0) .> 0
    dwatermask = dilate(watermask)
    gawm = GeoArray(watermask, ga.f, ga.crs, Dict{String,Any}())
    dgawm = GeoArray(dwatermask, ga.f, ga.crs, Dict{String,Any}())

    # Find all ICESat-2 data (stack of geobuffers)
    ices = Vector{DataFrame}()
    for f in readlines(icefn)
        df = GeoParquet.read(f)
        push!(ices, df)
    end
    ice = reduce(vcat, ices)
    if nrow(ice) > 0
        in_bbox!(ice, GeoArrays.bbox(ga))
        subset!(ice, :height_error => u -> u .< 5.0)
        subset!(ice, :height => z -> z .> lowlimit)
        subset!(ice, :cover => c -> c .!= 50)  # remove urban
    end
    if crossval && (nrow(ice) > 0)
        crossd = subset(ice, :track => t -> Base.contains.(t, "gt2"))
        subset!(ice, :track => t -> .!Base.contains.(t, "gt2"))
    end
    icemask, icevalues, icediff, iceind = burn_icesat(ga, ice, dwatermask)

    # Find all GEDI data (stack of geobuffers)
    gedis = Vector{DataFrame}()
    if isfile(gedifn)
        for f in readlines(gedifn)
            df = GeoParquet.read(f)
            push!(gedis, df)
        end
    end
    gedi = reduce(vcat, gedis)
    if nrow(gedi) > 0
        in_gbbox!(gedi, GeoArrays.bbox(ga))
        subset!(gedi, :sensitivity => s -> 0.95 .<= s .<= 1.0)
        subset!(gedi, :height => z -> z .> lowlimit)
    end
    gedimask, gedivalues, gedidiff, gediind = burn_gedi(ga, gedi, dwatermask)

    # Combine ICESat-2 and GEDI
    combined_mask = icemask .| gedimask
    gedivalues[icemask] .= icevalues[icemask]
    # dev = tileddeviation(gedivalues, w=50)
    dev = tileddiff(gedivalues, w=50)
    nnice = (dev .< -2.0) #.| (dev .> 2.0)
    combined_mask[nnice] .= false
    icemask[nnice] .= false

    # Filtered SpaceLiDAR points
    oind = filter!(x -> x > 0, gediind[nnice])
    @assert allunique(oind)
    gedi = deleteat!(gedi, sort(oind))
    oind = filter!(x -> x > 0, iceind[nnice])
    @assert allunique(oind)
    ice = deleteat!(ice, sort(oind))

    # Burn spaceborne lidar data
    covermask =
        in.(
            landcover.(cover),
            Ref((Bare, Cropland, Grassland, Moss, Snow, Water)),
        )
    lower_mask = gedivalues .<= A
    burnmask = combined_mask .& lower_mask .& .!covermask
    A[burnmask] .= gedivalues[burnmask]
    if debug
        gac = GeoArray(A, ga.f, ga.crs)
        gac = crop(gac, x, y, 0.0)
        GeoArrays.write(splitext(ofn)[1] * "_burn.tif", gac)
    end

    gedivalues[nnice] .= Inf
    gnnice = GeoArray(nnice, ga.f, ga.crs)
    glv = GeoArray(gedivalues, ga.f, ga.crs)
    if debug
        glv = crop(glv, x, y, 0.0)
        GeoArrays.write(splitext(ofn)[1] * "_lidar.tif", glv)
    end

    if nrow(ice) > 0
        ip = ice[:, [:longitude, :latitude, :height]]
        if nrow(gedi) > 0
            gp = gedi[:, [:longitude, :latitude, :height]]
            c = vcat(ip, gp)
        else
            c = ip
        end

        # Remove points that are removed by grid filter
        n = sample_mask(c, gnnice)
        c.nnice = n
        subset!(c, :nnice => x -> .!x)

        # Remove spacelidar values next to and on water
        v = sample_mask(c, gawm)
        c.mask = v
        subset!(c, :mask => x -> .!x)
        c.geometry = collect(zip(c.longitude, c.latitude))

        # Generate low-resolution DTM based on spaceborne lidar
        f = 15  # downsample factor
        exaggeration = 1  # exaggerate slope
        lwatermask = downsample(ocean_watermask, f, x -> all(x))
        lwaterval = downsample(gao[:, :, 1], f, minimum)
        w, h = round.(Int, size(ga)[1:2] ./ f, RoundUp)
        nga = GeoArray(Array{Union{Missing,Float32}}(missing, (w, h, 1)))
        bbox!(nga, GeoArrays.bbox(ga))
        GeoArrays.flipud!(nga)
        nga.crs = ga.crs
        lowbound = multiresdtm(nga, c, n=5; reducer=median, r2=mad)  # median of 5 points

        og = GeoArray(GeoArrayOps.mapwindow(median, lowbound[:, :, 1], (3, 3)), lowbound.f, lowbound.crs)
        ogdh0 = GeoArray(GeoArrayOps.mapwindow(median, lowbound[:, :, 2], (3, 3)), lowbound.f, lowbound.crs)
        ogf = copy(og)
        ogf[lwatermask, 1] .= lwaterval[lwatermask, 1]
        ogf[lwatermask, 1] .+= 2.0

        # Derive slope from low-resolution DTM
        sl = Float32.(tand.(slope(ogf[:, :, 1], cellsize=30 * f / exaggeration, method=MDG())))
        if debug
            ogf = crop(ogf, x, y, 0.0)
            GeoArrays.write(splitext(ofn)[1] * "_lowbound.tif", ogf)
        end

        # Resample to high-resolution
        lowbound = GeoArrays.warp(og, ga, Dict("r" => "bilinear"))
        lowbound[watermask, 1] .= -Inf
        if debug
            ogf = crop(lowbound, x, y, 0.0)
            GeoArrays.write(splitext(ofn)[1] * "_lowboundh.tif", ogf)
        end
        ss = imresize(sl, size(A), method=Interpolations.BSpline(Interpolations.Linear()))
        if debug
            gac = GeoArray(ss, ga.f, ga.crs, Dict{String,Any}())
            gac = crop(gac, x, y, 0.0)
            GeoArrays.write(splitext(ofn)[1] * "_slope.tif", gac)
        end
        dh0 = GeoArrays.warp(ogdh0, ga)[:, :, 1]
        if debug
            gac = GeoArray(dh0, ga.f, ga.crs, Dict{String,Any}())
            gac = crop(gac, x, y, 0.0)
            GeoArrays.write(splitext(ofn)[1] * "_dh0.tif", gac)
        end
    else
        @warn "No ICESat-2 or GEDI points found!"
        ss = fill(Float32(0.002), size(A))
        dh0 = fill(Float32(1.0), size(A))
        lowbound = GeoArray(zeros(size(A)), ga.f, ga.crs)
    end

    # Set inland waters to Inf, so not to influence filters
    A[watermask.&.!ocean_watermask, 1] .= Inf

    # Normalize DSM with low-resolution DTM
    X = copy(A)
    X .-= lowbound[:, :, 1]
    XX = copy(X)

    mindh0 = 1.0
    mins = 0.001

    # pattern = GeoArrayOps.TPI(normalized)
    cdh0 = clamp.(mindh0, dh0, Inf)
    pattern = clamp.(pattern, -cdh0 / 4, cdh0 / 4)

    # Run first filter for closed landcover
    B, flags, _ = GeoArrayOps.apsf(
        XX,
        ωₘ=2_000.0,
        slope=clamp.(mins, ss * 1, Inf),
        dhₘ=100.0,
        dh₀=clamp.(mindh0, dh0 * 1, Inf),
        cellsize=30.0,
        circular=true,
    )
    fm = X .> B

    if debug
        gax = GeoArray(fm, ga.f, ga.crs, Dict{String,Any}())
        gax = crop(gax, x, y, 0.0)
        GeoArrays.write(splitext(ofn)[1] * "_fm.tif", gax)
        gax = GeoArray(B, ga.f, ga.crs, Dict{String,Any}())
        GeoArrays.write(splitext(ofn)[1] * "_fmB.tif", gax)
        # gax = GeoArray(dmask, ga.f, ga.crs, Dict{String,Any}())
        # GeoArrays.write(splitext(ofn)[1] * "_division.tif", gax)
    end

    # Set closed landcover to zero in normalized DSM
    closed_landcover_mask = fm .& .!covermask
    X[closed_landcover_mask] .= 0
    XX = clamp.(X, 0, Inf)

    if debug
        gax = GeoArray(XX, ga.f, ga.crs, Dict{String,Any}())
        gax = crop(gax, x, y, 0.0)
        GeoArrays.write(splitext(ofn)[1] * "_norm.tif", gax)
    end

    # Run second filter for open landcover
    B, flags = GeoArrayOps.apsf(
        XX,
        ωₘ=1_000.0,
        slope=clamp.(mins, ss * 20, Inf),
        dhₘ=100.0,
        dh₀=clamp.(mindh0, dh0 * 2, Inf),
        cellsize=30.0,
        circular=true,
    )
    sfm = X .> B

    if debug
        gax = GeoArray(sfm, ga.f, ga.crs, Dict{String,Any}())
        gax = crop(gax, x, y, 0.0)
        GeoArrays.write(splitext(ofn)[1] * "_sfm.tif", gax)
        gax = GeoArray(B, ga.f, ga.crs, Dict{String,Any}())
        GeoArrays.write(splitext(ofn)[1] * "_sfmB.tif", gax)
    end

    urbanmask = landcover.(cover) .== Urban
    wc_watermask = landcover.(cover) .== Water

    # Quality is ICESat-2 burned points and low error CopernicusDEM
    quality_mask = icemask .& .!covermask .& .!urbanmask  # ICESat-2 can't penetrate buildings

    # Error high error and nodata (filled) CopernicusDEM
    # but only for the cells that are not covered by spaceborne lidar
    error_mask = (ea .>= 0.75) .| isinf.(ea) .& .!combined_mask

    # Create mask of non-terrain points by combining all masks
    m = nm .& ((fm .& .!covermask) .| error_mask .| (low_mask .& .!combined_mask) .| sfm)
    m .&= .!quality_mask  # values to interpolate

    # Keep burned ICESat-2 values in removed open terrain
    pbm = m .& icemask .& covermask .& lower_mask
    A[pbm] .= gedivalues[pbm]
    m .&= .!pbm

    # Dont't interpolate values higher than 100m
    m .&= .!high_mask  # values to interpolate

    A[m] .= Inf
    A[watermask, 1] .= Inf

    gac = GeoArray(A, ga.f, ga.crs, Dict{String,Any}())
    if debug
        gacs = crop(gac, x, y, 0.0)
        GeoArrays.write(splitext(ofn)[1] * "_f.tif", gacs)
    end

    # Fill filtered out cells with AIDW interpolation
    mm = reshape3d(m)
    C = interpolate(AIDW(knn=25, power=1), gac, mm)

    # Add back TPI based noise
    C[mm] .+= reshape3d(pattern)[mm]

    # Remove values higher than original
    too_high = C .> reshape3d(orig)
    C[too_high] .= reshape3d(orig)[too_high]

    # Put back watervalues
    C[watermask, 1] .= Inf

    # Crop to actual 1x1 degree center
    gaf = GeoArray(C.A, ga.f, ga.crs, Dict{String,Any}())
    if !isnothing(cropsize)
        gafc = crop(gaf, x, y, 0.0)
        cx, cy = size(gaf)[1:2] .% 2 .!= 0
        gafc = gafc[1:end-cx, 1:end-cy]  # cutoff extra row/column if exists
    else
        gafc = gaf
    end

    if crossval && (nrow(ice) > 0)
        sample!(crossd, gafc)
        dropmissing!(crossd, :csample)
        if nrow(crossd) > 0
            GeoParquet.write(splitext(ofn)[1] * "_crossval.pq", crossd, (:geometry,))
        end
    end

    nodata = -9999
    gafc[isinf.(gafc)] .= nodata
    return GeoArrays.write(ofn, gafc; nodata=nodata, shortname="COG", options=Dict("COMPRESS" => "ZSTD", "PREDICTOR" => "3"))
end

function coalesce2d(ga, v)
    dropdims(coalesce(ga, v), dims=3)
end

function reshape3d(A::AbstractMatrix)
    reshape(A, (size(A)..., 1))
end

function reshape2d(A::AbstractArray)
    dropdims(A, dims=3)
end

function in_bbox!(df::DataFrame, bbox)
    subset!(df, :longitude => x -> (bbox.min_x .<= x .<= bbox.max_x), :latitude => y -> (bbox.min_y .<= y .<= bbox.max_y))
end

function in_bbox(df::DataFrame, bbox)
    subset(df, :longitude => x -> (bbox.min_x .<= x .<= bbox.max_x), :latitude => y -> (bbox.min_y .<= y .<= bbox.max_y))
end

function in_gbbox!(df::DataFrame, bbox)
    subset!(df, :longitude => x -> (bbox.min_x .<= x .<= bbox.max_x), :latitude => y -> (bbox.min_y .<= y .<= bbox.max_y))
end
