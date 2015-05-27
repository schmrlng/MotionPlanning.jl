export EuclideanNN_KDTree, MetricNN_BruteForce, ArcLength_Pruned

### Metric brute force
type MetricNN_BruteForce{S<:State,T<:FloatingPoint,U<:ControlInfo} <: MetricNN
    V::Vector{S}
    dist::Metric
    DS::Matrix{T}
    US::Matrix{U}
    cache::Vector{Neighborhood{T}}
    kNNr::Vector{T}

    function MetricNN_BruteForce(V::Vector{S}, dist::Metric = Euclidean())
        NN = new(V, dist)
        initcache(NN)
    end
end
function initcache{S,T,U}(NN::MetricNN_BruteForce{S,T,U})
    N = length(NN.V)
    NN.DS, NN.US = pairwise_distances(NN.dist, NN.V)
    NN.cache, NN.kNNr = Array(Neighborhood{T}, N), zeros(T, N)
    NN
end
MetricNN_BruteForce{S<:State}(V::Vector{S}, dist::Metric = Euclidean()) =
    MetricNN_BruteForce{S,eltype(S),controltype(dist)}(V,dist) # so constructor works without {}

function inball(NN::MetricNN_BruteForce, v::Int, r)
    @devec nn_bool = NN.DS[:,v] .<= r
    nn_bool[v] = false
    inds = find(nn_bool)
    Neighborhood(inds, NN.DS[inds,v])
end

function knn(NN::MetricNN_BruteForce, v::Int, k = 1)    # ds sorted increasing
    r = select!(NN.DS[:,v], k+1)     # TODO: could be faster with an ArrayView and leveraging quickselect algo
    nbhd = inball(NN, v, r)
    p = sortperm(nbhd.ds)
    permute!(nbhd.inds, p)
    permute!(nbhd.ds, p)
    nbhd
end

### Euclidean k-d tree
type EuclideanNN_KDTree{S<:AbstractVector,T<:FloatingPoint,U<:ControlInfo} <: MetricNN
    V::Vector{S}
    dist::Metric
    DS::KDTree
    cache::Vector{Neighborhood{T}}
    kNNr::Vector{T}

    function EuclideanNN_KDTree(V::Vector{S}, dist::Metric = Euclidean())
        dist != Euclidean() && error("Distance metric must be Euclidean for EuclideanNN_KDTree")
        NN = new(V, Euclidean())
        initcache(NN)
    end
end
function initcache{S,T,U}(NN::EuclideanNN_KDTree{S,T,U})
    N = length(NN.V)
    NN.DS = KDTree(hcat(NN.V...))  # TODO: leafsize, reorder?
    NN.cache, NN.kNNr = Array(Neighborhood{T}, N), zeros(T, N)
    NN
end
EuclideanNN_KDTree{S<:AbstractVector}(V::Vector{S}, dist::Metric = Euclidean()) =
    EuclideanNN_KDTree{S,eltype(S),controltype(dist)}(V,dist) # so constructor works without {}

function inball{S,T,U}(NN::EuclideanNN_KDTree{S,T,U}, v::Int, r)
    inds = KDTrees.inball(NN.DS, convert(Vector{T}, NN[v]), r, true)
    inds = deleteat!(inds, searchsortedfirst(inds, v))
    Neighborhood(inds, colwise(Euclidean(), NN[v], hcat(NN[inds]...)))
end

function knn{S,T,U}(NN::EuclideanNN_KDTree{S,T,U}, v::Int, k = 1)    # ds sorted increasing
    inds, ds = KDTrees.knn(NN.DS, convert(Vector{T}, NN[v]), k+1)
    shift!(inds)
    shift!(ds)
    Neighborhood(inds, ds)
end

### Arc length pruned KDTree
type ArcLength_Pruned{S<:State,T<:FloatingPoint,U<:ControlInfo} <: MetricNN
    V::Vector{S}
    dist::Metric
    DS::KDTree
    US::SparseMatrixCSC{U,Int}
    cache::Vector{Neighborhood{T}}
    kNNr::Vector{T}

    function ArcLength_Pruned(V::Vector{S}, dist::Metric = Euclidean())
        NN = new(V, dist)
        initcache(NN)
    end
end
function initcache{S,T,U}(NN::ArcLength_Pruned{S,T,U})
    N = length(NN.V)
    NN.DS = KDTree(hcat(NN.V...))  # TODO: leafsize, reorder?
    NN.US = sparse([],[],U[],N,N)
    NN.cache, NN.kNNr = Array(Neighborhood{T}, N), zeros(T, N)
    NN
end
ArcLength_Pruned{S<:State}(V::Vector{S}, dist::Metric = Euclidean()) =
    ArcLength_Pruned{S,eltype(S),controltype(dist)}(V,dist) # so constructor works without {}

function inball{S,T,U}(NN::ArcLength_Pruned{S,T,U}, v::Int, r)
    inds = KDTrees.inball(NN.DS, convert(Vector{T}, NN[v]), r, true)
    inds = deleteat!(inds, searchsortedfirst(inds, v))
    ds = T[v < i ? evaluate(NN.dist, NN[v], NN[i]) : evaluate(NN.dist, NN[i], NN[v])   # in case metric evaluation isn't quite symmetric (e.g. interpolated)
           for i in inds]
    @devec pruned = ds .<= r
    Neighborhood(inds[pruned], ds[pruned])
end

function knn{S,T,U}(NN::ArcLength_Pruned{S,T,U}, v::Int, k = 1)    # ds sorted increasing
    error("NN datastructure ArcLength_Pruned does not support k-nearest neighbors. Try brute force?")
end

# type EuclideanNN_FLANN{S<:AbstractVector,T<:FloatingPoint} <: MetricNN
#     V::Vector{S}
#     DS::FLANNIndex
#     cache::Vector{Neighborhood{T}}
#     kNNr::Vector{T}

#     function EuclideanNN_FLANN(V::Vector{S})
#         NN = new(V, flann(hcat(V...), FLANNParameters()), Array(Neighborhood{T}, length(V)), zeros(T, length(V)))
#         finalizer(NN, x -> close(x.DS))
#         NN
#     end
# end

# type CachedMetricNN{S<:State,T<:FloatingPoint,U<:ControlInfo} <: MetricNN  # TODO: refactor cache -> ball vs. knn
#     V::Vector{S}
#     dist::Metric
#     cache::Vector{Neighborhood{T}}
#     kNNr::Vector{T}
# end
# function CachedMetricNN(NN::MetricNN)
#     # with above TODO: make sure cache is populated
#     CachedMetricNN(NN.V, NN.dist, NN.cache, NN.kNNr)
# end
# function CachedMetricNN(fname::String)
#     open(fname, "r") do f
#         deserialize(f)
#     end
# end
# inball(NN::CachedMetricNN{S,T,U}, v::Int, r) = error("Online inball eval in CachedMetricNN: check intialization?")
# knn(NN::CachedMetricNN{S,T,U}, v::Int, k) = error("Online knn eval in CachedMetricNN: check intialization?")
# function saveNN(NN::MetricNN, fname::String)

### General Methods
function inball!(NN::MetricNN, v::Int, r)
    if !isdefined(NN.cache, v)
        NN.cache[v] = inball(NN, v, r)
        NN.kNNr[v] = r
    end
    return NN.cache[v]
end

function knn!(NN::MetricNN, v::Int, k = 1)
    if !isdefined(NN.cache, v)
        nbhd = knn(NN, v, k)
        NN.cache[v] = nbhd
        NN.kNNr[v] = nbhd.ds[end]
    end
    return NN.cache[v]
end

inball!(NN::MetricNN, v::Int, r, f::BitVector) = filter_neighborhood(inball!(NN, v, r), f)
knn!(NN::MetricNN, v::Int, k, f::BitVector) = filter_neighborhood(knn!(NN, v, k), f)
function mutualknn!(NN::MetricNN, v::Int, k, f::BitVector)
    n = knn!(NN, v, k)
    for w in n.inds
        knn!(NN, w, k)       # to ensure neighbors' respective knn sets have been computed
    end
    @devec mutual_inds = f[n.inds] & (n.ds .<= NN.kNNr[n.inds])
    return Neighborhood(n.inds[mutual_inds], n.ds[mutual_inds])
end

### Forwards- and backwards-NN are the same for metric spaces (TODO: figure out how metaprogramming works)
inballF(NN::MetricNN, args...) = inball(NN, args...)
knnF(NN::MetricNN, args...) = knn(NN, args...)
inballF!(NN::MetricNN, args...) = inball!(NN, args...)
knnF!(NN::MetricNN, args...) = knn!(NN, args...)
mutualknnF!(NN::MetricNN, args...) = mutualknn!(NN, args...)
inballB(NN::MetricNN, args...) = inball(NN, args...)
knnB(NN::MetricNN, args...) = knn(NN, args...)
inballB!(NN::MetricNN, args...) = inball!(NN, args...)
knnB!(NN::MetricNN, args...) = knn!(NN, args...)
mutualknnB!(NN::MetricNN, args...) = mutualknn!(NN, args...)