abstract QuasiMetric <: PreMetric   # positivity, positive definiteness, and triangle inequality (no symmetry)

### Quasimetric brute force
type QuasiMetricNN_BruteForce{S<:State,T<:FloatingPoint,U<:ControlInfo} <: QuasiMetricNN
    V::Vector{S}
    DS::Matrix{T}
    cacheF::Vector{Neighborhood{T}}
    cacheB::Vector{Neighborhood{T}}
    kNNrF::Vector{T}
    kNNrB::Vector{T}
    dist::QuasiMetric
    US::Matrix{U}

    function QuasiMetricNN_BruteForce(V::Vector{S}, dist::QuasiMetric)
        N = length(V)
        new(V, pairwise(dist, hcat(V...)), 
            Array(Neighborhood{T}, N),
            Array(Neighborhood{T}, N),
            zeros(T, N),
            zeros(T, N),
            dist,
            Array(u, N, N))
    end
end
QuasiMetricNN_BruteForce{S<:State}(V::Vector{S}, dist::QuasiMetric) =
    QuasiMetricNN_BruteForce{S,eltype(S),controltype(dist)}(V,dist) # so constructor works without {}

## Forwards

function inballF(NN::QuasiMetricNN_BruteForce, v::Int, r)
    @devec nn_bool = NN.DS[v,:] .<= r
    nn_bool[v] = false
    inds = find(nn_bool)
    Neighborhood(inds, NN.DS[v,inds])
end

function knnF(NN::QuasiMetricNN_BruteForce, v::Int, k = 1)    # ds sorted increasing
    r = select!(vec(NN.DS[v,:]), k+1)     # TODO: could be faster with an ArrayView and leveraging quickselect algo
    nbhd = inballF(NN, v, r)
    p = sortperm(nbhd.ds)
    permute!(nbhd.inds, p)
    permute!(nbhd.ds, p)
    nbhd
end

## Backwards

function inballB(NN::QuasiMetricNN_BruteForce, v::Int, r)
    @devec nn_bool = NN.DS[:,v] .<= r
    nn_bool[v] = false
    inds = find(nn_bool)
    Neighborhood(inds, NN.DS[inds,v])
end

function knnB(NN::QuasiMetricNN_BruteForce, v::Int, k = 1)    # ds sorted increasing
    r = select!(NN.DS[:,v], k+1)     # TODO: could be faster with an ArrayView and leveraging quickselect algo
    nbhd = inball(NN, v, r)
    p = sortperm(nbhd.ds)
    permute!(nbhd.inds, p)
    permute!(nbhd.ds, p)
    nbhd
end

### General Methods (HUGE TODO: learn metaprogramming)

## Forwards

function inballF!(NN::QuasiMetricNN, v::Int, r)
    if !isdefined(NN.cacheF, v)
        NN.cacheF[v] = inballF(NN, v, r)
        NN.kNNrF[v] = r
    end
    return NN.cacheF[v]
end

function knnF!(NN::QuasiMetricNN, v::Int, k = 1)
    if !isdefined(NN.cacheF, v)
        nbhd = knnF(NN, v, k)
        NN.cacheF[v] = nbhd
        NN.kNNrF[v] = nbhd.ds[end]
    end
    return NN.cacheF[v]
end

inballF!(NN::QuasiMetricNN, v::Int, r, f::BitVector) = filter_neighborhood(inballF!(NN, v, r), f)
knnF!(NN::QuasiMetricNN, v::Int, k, f::BitVector) = filter_neighborhood(knnF!(NN, v, k), f)
function mutualknnF!(NN::QuasiMetricNN, v::Int, k, f::BitVector)
    n = knnF!(NN, v, k)
    for w in n.inds
        knnB!(NN, w, k)       # to ensure neighbors' respective knn sets have been computed
    end
    @devec mutual_inds = f[n.inds] & (n.ds .<= NN.kNNrB[n.inds])
    return Neighborhood(n.inds[mutual_inds], n.ds[mutual_inds])
end

## Backwards

function inballB!(NN::QuasiMetricNN, v::Int, r)
    if !isdefined(NN.cacheB, v)
        NN.cacheB[v] = inballB(NN, v, r)
        NN.kNNrB[v] = r
    end
    return NN.cacheB[v]
end

function knnB!(NN::QuasiMetricNN, v::Int, k = 1)
    if !isdefined(NN.cacheB, v)
        nbhd = knnB(NN, v, k)
        NN.cacheB[v] = nbhd
        NN.kNNrB[v] = nbhd.ds[end]
    end
    return NN.cacheB[v]
end

inballB!(NN::QuasiMetricNN, v::Int, r, f::BitVector) = filter_neighborhood(inballB!(NN, v, r), f)
knnB!(NN::QuasiMetricNN, v::Int, k, f::BitVector) = filter_neighborhood(knnB!(NN, v, k), f)
function mutualknnB!(NN::QuasiMetricNN, v::Int, k, f::BitVector)
    n = knnB!(NN, v, k)
    for w in n.inds
        knnF!(NN, w, k)       # to ensure neighbors' respective knn sets have been computed
    end
    @devec mutual_inds = f[n.inds] & (n.ds .<= NN.kNNrF[n.inds])
    return Neighborhood(n.inds[mutual_inds], n.ds[mutual_inds])
end