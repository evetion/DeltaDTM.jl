using NearestNeighbors
using Proj

const width_of_degree_in_km = 111.32

function calc_angles(coords, coordinates, work)
    for (i, c) in enumerate(eachcol(coordinates))
        work[i] = atand(c[2] - coords[2], c[1] - coords[1]) % 180
    end
    return work
end

function calc_k(
    coords,
    distances,
    coordinates,
    work,
    w,
    metric=Euclidean(),
)::Vector{Float64}
    # assumes these coordinates and angles are already sorted
    # with respect to their distance to the point to interpolate
    angles = calc_angles(coords, coordinates, work)
    shield_angle = 360.0 / length(w)
    w .= 1.0

    # For each point, check if closerby points shield it
    # and if so calculate a lower weight for it
    for i in 2:length(w)
        for j in 1:i-1
            α = abs(angles[i] - angles[j])
            if α == 0
                w[i] = 0
                break
            elseif α < shield_angle
                ci = @view coordinates[:, i]
                cj = @view coordinates[:, j]
                dij = Distances.evaluate(metric, ci, cj)
                if dij == 0  # i is directly behind or on j
                    w[i] = 0
                    break
                end
                θ = distances[j] * sind(α / 2) / (dij / 2)
                w[i] *= abs(θ)
            end
        end
    end
    return w
end

function aidw(img, I, patch)
    pc = Tuple.(patch)
    distances = euclidean.(pc, Ref(I.I))
    distancesv = reshape(distances, :)
    s = sortperm(distancesv)
    coordinates = hcat(collect.(reshape(pc, :)[s])...)::Matrix{Int64}
    k = calc_k(distancesv[s], I.I, coordinates)
    w = k ./ distancesv[s]

    A = reshape(view(img, patch), :)[s]
    m = isfinite.(A)
    if sum(m) == 0
        return NaN
    else
        Σw = sum(w[m])
        w[m] = w[m] / Σw
        return sum(A[m] .* w[m])
    end
end

function _aidw(coord, dists, nb, values, work, work2, p=2)
    k = calc_k(coord, dists, nb, work, work2)
    Σw = 0.0
    out = 0.0
    for I in eachindex(k)
        k[I] = k[I] / dists[I]^p
        Σw += k[I]
    end
    for I in eachindex(k)
        out += values[I] * k[I] / Σw
    end
    return out
end

function _idw(dists, values, p=2)
    w = dists .^ p
    Σw = sum(w)
    w /= Σw
    return sum(values .* w)
end


trans = Proj.Transformation("EPSG:4326", "EPSG:3857", always_xy=true)

function straight2!(X)
    X .= trans(X)
end

function project!(data)
    size(data)[1] == 2 || error("data must be 2D")
    for col in eachcol(data)
        @inbounds straight2!(col)
    end
end

abstract type AbstractInterpolation end


Base.@kwdef struct AIDW{N<:Integer,P<:AbstractFloat} <: AbstractInterpolation
    knn::N = 25
    power::P = 2.0
end


function interpolate(i::AbstractInterpolation, A, m)
    return interpolate!(i, copy(A), m)
end

"""
    interpolate!(::AIDW, A, mask)

Interpolate using the AIDW interpolation method. This uses a KDTree based on the finite
data points in `A` and interpolates the points in `mask`.

Adjusted IDW interpolation limits the weight of points that are obscured by other points.
"""
function interpolate!(a::AIDW, A, mask)
    coordinates = collect(GeoArrays.coords(A))

    # Setup valid input locations
    datamask3 = isfinite.(A)
    any(datamask3) || return A
    datamask = dropdims(datamask3; dims=3)
    data = reinterpret(reshape, Float64, coordinates[datamask])
    project!(data)

    k = min(a.knn, size(data)[2])

    # Setup output data locations
    outmask3 = mask .== 1.0
    any(outmask3) || return A
    outmask = dropdims(outmask3; dims=3)
    points = reinterpret(reshape, Float64, coordinates[outmask])
    project!(points)

    tree = KDTree(data, Euclidean())
    (idxs, dists) = knn(tree, points, k, true)

    zp = zeros(size(points)[2])
    z = A[datamask3]
    work = zeros(k)  # for angles
    work2 = zeros(k)  # for weights
    @inbounds for i in eachindex(idxs)
        coord = @view points[:, i]
        nb = @view data[:, idxs[i]]
        distances = dists[i]
        values = @view z[idxs[i]]
        zp[i] = _aidw(coord, distances, nb, values, work, work2, a.power)
    end
    A[outmask3] .= zp
    return A
end
