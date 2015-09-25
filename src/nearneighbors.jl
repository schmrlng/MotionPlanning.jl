import Base.length, Base.getindex
export NearNeighborCache, MetricNN, QuasiMetricNN, Neighborhood, NNDatastructureCache
export filter_neighborhood, addpoints, initcache, saveNN, loadNN!
export inball, knn, inball!, knn!, mutualknn!         # for metrics only (symmetric)
export inballF, knnF, inballF!, knnF!, mutualknnF!    # forwards
export inballB, knnB, inballB!, knnB!, mutualknnB!    # backwards

abstract NearNeighborCache
abstract MetricNN <: NearNeighborCache
abstract QuasiMetricNN <: NearNeighborCache

type Neighborhood{T<:AbstractFloat}
    inds::Vector{Int}
    ds::Vector{T}
end
filter_neighborhood(n::Neighborhood, f::BitVector) = Neighborhood(n.inds[f[n.inds]], n.ds[f[n.inds]])  # TODO: this is slow (for loop to get rid of double-scan, or perhaps best, use an iterator)

addpoints(NN::NearNeighborCache, W) = typeof(NN)([NN.V; W], NN.dist, NN.init)
length(NN::NearNeighborCache) = length(NN.V)
getindex(NN::NearNeighborCache, i::Int) = i > 0 ? NN.V[i] : NN.init
getindex(NN::NearNeighborCache, I::AbstractVector) = [NN[i] for i in I]

type NNDatastructureCache
    V
    DS
    US
end
NNDatastructureCache(NN::NearNeighborCache) = NNDatastructureCache(NN.V, NN.DS, NN.US)
function NNDatastructureCache(fname::AbstractString)
    open(fname, "r") do f
        deserialize(f)
    end
end
function saveNN(DC::NNDatastructureCache, fname::AbstractString)
    open(fname, "w") do f
        serialize(f, DC)
    end
end
saveNN(NN::NearNeighborCache) = saveNN(NNDatastructureCache(NN))

include("nearneighbors/metricNN.jl")
include("nearneighbors/quasimetricNN.jl")

# TODO: addpoints! function
# TODO: precompute_all function
# TODO: clean up huge amounts of copy/paste code with metaprogramming