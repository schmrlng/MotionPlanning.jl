export EuclideanNN_KDTree, MetricNN_BruteForce

### Metric brute force
type MetricNN_BruteForce{S<:State,T<:FloatingPoint,U<:ControlInfo} <: MetricNN
    V::Vector{S}
    DS::Matrix{T}
    cache::Vector{Neighborhood{T}}
    kNNr::Vector{T}
    dist::Metric
    US::Matrix{U}

    function MetricNN_BruteForce(V::Vector{S}, dist::Metric = Euclidean())
        N = length(V)
        DS, US = pairwise_distances(dist, V)
        new(V, DS, Array(Neighborhood{T}, N), zeros(T, N), dist, US)
    end
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
    DS::KDTree
    cache::Vector{Neighborhood{T}}
    kNNr::Vector{T}
    dist::Metric
    US::Matrix{U}

    function EuclideanNN_KDTree(V::Vector{S}, dist::Metric = Euclidean())
        N = length(V)
        dist != Euclidean() && error("Distance metric must be Euclidean for EuclideanNN_KDTree")
        new(V, KDTree(hcat(V...)), Array(Neighborhood{T}, N), zeros(T, N), Euclidean(), Array(NullControl, 0, 0)) # TODO: leafsize, reorder?
    end
end
EuclideanNN_KDTree{S<:State}(V::Vector{S}, dist::Metric = Euclidean()) =
    EuclideanNN_KDTree{S,eltype(S),controltype(dist)}(V,dist) # so constructor works without {}

function inball{S,T}(NN::EuclideanNN_KDTree{S,T}, v::Int, r)
    inds = KDTrees.inball(NN.DS, convert(Vector{T}, NN.V[v]), r, true)
    inds = deleteat!(inds, searchsortedfirst(inds, v))
    Neighborhood(inds, colwise(Euclidean(), NN.V[v], NN.v[inds]))
end

function knn{S,T}(NN::EuclideanNN_KDTree{S,T}, v::Int, k = 1)    # ds sorted increasing
    inds, ds = KDTrees.knn(NN.DS, convert(Vector{T}, NN.V[v]), k+1)
    shift!(inds)
    shift!(ds)
    Neighborhood(inds, ds)
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