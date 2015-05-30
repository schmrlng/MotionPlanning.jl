export QuasiMetric, QuasiMetricNN_BruteForce

abstract QuasiMetric <: PreMetric   # positivity, positive definiteness, and triangle inequality (no symmetry)

### Quasimetric brute force
type QuasiMetricNN_BruteForce{S<:State,T<:FloatingPoint,U<:ControlInfo} <: QuasiMetricNN
    V::Vector{S}
    dist::QuasiMetric
    DS::SparseMatrixCSC{T,Int}
    US::SparseMatrixCSC{U,Int}
    cacheF::Vector{Neighborhood{T}}
    cacheB::Vector{Neighborhood{T}}
    kNNrF::Vector{T}
    kNNrB::Vector{T}

    function QuasiMetricNN_BruteForce(V::Vector{S}, dist::QuasiMetric)
        NN = new(V, dist)
        initcache(NN)
    end
end
function initcache{S,T,U}(NN::QuasiMetricNN_BruteForce{S,T,U})
    N = length(NN.V)
    NN.DS, NN.US = pairwise_distances(NN.dist, NN.V)
    NN.cacheF, NN.kNNrF = Array(Neighborhood{T}, N), zeros(T, N)
    NN.cacheB, NN.kNNrB = Array(Neighborhood{T}, N), zeros(T, N)
    NN
end
QuasiMetricNN_BruteForce{S<:State}(V::Vector{S}, dist::QuasiMetric) =
    QuasiMetricNN_BruteForce{S,eltype(S),controltype(dist)}(V,dist) # so constructor works without {}
function loadNN!{S,T,U}(NN::QuasiMetricNN_BruteForce{S,T,U}, NNDC::NNDatastructureCache, init::State, goal::Goal, CC::CollisionChecker)
    filter = Bool[is_free_state(v, CC) for v in NNDC.V]
    N = sum(filter)
    NN.V = NNDC.V[filter]
    NN.DS = NNDC.DS[filter,filter]
    NN.US = NNDC.US[filter,filter]
    if NN.V[1] != init
        # TODO
    end
    NN.cacheF, NN.kNNrF = Array(Neighborhood{T}, N), zeros(T, N)
    NN.cacheB, NN.kNNrB = Array(Neighborhood{T}, N), zeros(T, N)
    # TODO: this info should be cached too but it's late and this hack will work
    r = maximum(NN.DS)
    for v in 1:N
        inballF!(NN,v,r)
        inballB!(NN,v,r)
    end
    NN
end

## Forwards

function inballF(NN::QuasiMetricNN_BruteForce, v::Int, r)
    @devec nn_bool = (0 .< NN.DS[v,:] .<= r)
    nn_bool[v] = false
    inds = find(nn_bool)
    Neighborhood(inds, vec(full(NN.DS[v,inds])))
end

function knnF(NN::QuasiMetricNN_BruteForce, v::Int, k = 1)    # ds sorted increasing
    error("TODO: this code needs a fix to work with sparse DS")
    r = select!(vec(NN.DS[v,:]), k+1)     # TODO: could be faster with an ArrayView and leveraging quickselect algo
    nbhd = inballF(NN, v, r)
    p = sortperm(nbhd.ds)
    permute!(nbhd.inds, p)
    permute!(nbhd.ds, p)
    nbhd
end

## Backwards

function inballB(NN::QuasiMetricNN_BruteForce, v::Int, r)
    @devec nn_bool = (0 .< NN.DS[:,v] .<= r)
    nn_bool[v] = false
    inds = find(nn_bool)
    Neighborhood(inds, vec(full(NN.DS[inds,v])))
end

function knnB(NN::QuasiMetricNN_BruteForce, v::Int, k = 1)    # ds sorted increasing
    error("TODO: this code needs a fix to work with sparse DS")
    r = select!(NN.DS[:,v], k+1)     # TODO: could be faster with an ArrayView and leveraging quickselect algo
    nbhd = inballB(NN, v, r)
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