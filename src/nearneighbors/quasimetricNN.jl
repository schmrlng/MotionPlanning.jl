export QuasiMetric, QuasiMetricNN_BruteForce, generate_QM_cache

abstract QuasiMetric <: PreMetric   # positivity, positive definiteness, and triangle inequality (no symmetry)

### Quasimetric brute force
type QuasiMetricNN_BruteForce{S<:State,T<:FloatingPoint,U<:ControlInfo} <: QuasiMetricNN
    V::Vector{S}
    dist::QuasiMetric
    init::S
    DS::SparseMatrixCSC{T,Int}
    US::SparseMatrixCSC{U,Int}
    cacheF::Vector{Neighborhood{T}}
    cacheB::Vector{Neighborhood{T}}
    kNNrF::Vector{T}
    kNNrB::Vector{T}

    function QuasiMetricNN_BruteForce(V::Vector{S}, dist::QuasiMetric, init::S)
        NN = new(V, dist, init)
        initcache(NN)
    end
end
function initcache{S,T,U}(NN::QuasiMetricNN_BruteForce{S,T,U})
    N = length(NN.V)
    if N > 0
        NN.DS, NN.US = pairwise_distances(NN.dist, NN.V)
    end
    NN.cacheF, NN.kNNrF = Array(Neighborhood{T}, N), zeros(T, N)
    NN.cacheB, NN.kNNrB = Array(Neighborhood{T}, N), zeros(T, N)
    NN
end
QuasiMetricNN_BruteForce{S<:State}(V::Vector{S}, dist::QuasiMetric, init::S) =
    QuasiMetricNN_BruteForce{S,eltype(S),controltype(dist)}(V,dist,init) # so constructor works without {}
function loadNN!{S,T,U}(NN::QuasiMetricNN_BruteForce{S,T,U}, NNDC::NNDatastructureCache, init::State, goal::Goal, CC::CollisionChecker)
    filter = Bool[is_free_state(v, CC) for v in NNDC.V]
    N = sum(filter)
    NN.V = NNDC.V[filter]
    NN.DS = NNDC.DS[filter,filter]
    NN.US = NNDC.US[filter,filter]
    if NN.V[1] != init
        # TODO, or make this part of fmtstar!, or even sample_free!
        # also goal check
        # also Ross's ML stuff is maybe worth a shot if dist eval is super duper slow
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
function loadNN!{S,T,U}(NN::QuasiMetricNN_BruteForce{S,T,U}, NNDC::NNDatastructureCache)
    N = length(NNDC.V)
    NN.V = NNDC.V
    NN.DS = NNDC.DS
    NN.US = NNDC.US
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
function generate_QM_cache(SS, N, init, goal, r, fname = Pkg.dir("MotionPlanning")*"/data/DI_$(N)_$(r)")  # TODO: init, goal only here as a stopgap until quasimetricNN.jl:37 addressed
    V = [sample_space(SS) for i in 1:N]
    V[1] = init
    V[end] = vector_to_state(sample_goal(goal), SS)
    setup_steering(SS, r)
    NNDC = NNDatastructureCache(QuasiMetricNN_BruteForce(V, SS.dist))
    saveNN(NNDC, fname)
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

### Quasimetric LB by Euclidean (TODO: general bounded NN, and really a redesign of everything NN)
type QMArcLength_Pruned{S<:State,T<:FloatingPoint,U<:ControlInfo} <: QuasiMetricNN
    V::Vector{S}
    dist::QuasiMetric
    init::S
    DS::KDTree
    US::SparseMatrixCSC{U,Int}      # TODO: switch to a sparse matrix structure better suited for insertion, at least while populating
    cacheF::Vector{Neighborhood{T}}
    cacheB::Vector{Neighborhood{T}}
    kNNrF::Vector{T}
    kNNrB::Vector{T}

    function QMArcLength_Pruned(V::Vector{S}, dist::QuasiMetric, init::S)
        NN = new(V, dist, init)
        initcache(NN)
    end
end
function initcache{S,T,U}(NN::QMArcLength_Pruned{S,T,U})
    N = length(NN.V)
    N > 0 && (NN.DS = KDTree(hcat(NN.V...)))  # TODO: leafsize, reorder?
    NN.US = sparse([],[],U[],N,N)
    NN.cacheF, NN.kNNrF = Array(Neighborhood{T}, N), zeros(T, N)
    NN.cacheB, NN.kNNrB = Array(Neighborhood{T}, N), zeros(T, N)
    NN
end
QMArcLength_Pruned{S<:State}(V::Vector{S}, dist::QuasiMetric, init::S) =
    QMArcLength_Pruned{S,eltype(S),controltype(dist)}(V,dist,init) # so constructor works without {}

## Forwards

function inballF{S,T,U}(NN::QMArcLength_Pruned{S,T,U}, v::S, r)
    inds = KDTrees.inball(NN.DS, convert(Vector{T}, v), r, true)
    ds = T[evaluate(NN.dist, v, NN[i]) for i in inds]
    @devec pruned = ds .<= r
    Neighborhood(inds[pruned], ds[pruned])
end

function inballF{S,T,U}(NN::QMArcLength_Pruned{S,T,U}, v::Int, r)
    inds = KDTrees.inball(NN.DS, convert(Vector{T}, NN[v]), r, true)
    inds = deleteat!(inds, searchsortedfirst(inds, v))
    ds = T[evaluate(NN.dist, NN[v], NN[i]) for i in inds]
    @devec pruned = ds .<= r
    Neighborhood(inds[pruned], ds[pruned])
end

function knnF{S,T,U}(NN::QMArcLength_Pruned{S,T,U}, v::Int, k = 1)    # ds sorted increasing
    error("TODO")
end

## Backwards

function inballB{S,T,U}(NN::QMArcLength_Pruned{S,T,U}, v::Int, r)
    inds = KDTrees.inball(NN.DS, convert(Vector{T}, NN[v]), r, true)
    inds = deleteat!(inds, searchsortedfirst(inds, v))
    ds = T[evaluate(NN.dist, NN[i], NN[v]) for i in inds]
    @devec pruned = ds .<= r
    Neighborhood(inds[pruned], ds[pruned])
end

function knnB{S,T,U}(NN::QMArcLength_Pruned{S,T,U}, v::Int, k = 1)    # ds sorted increasing
    error("TODO")
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