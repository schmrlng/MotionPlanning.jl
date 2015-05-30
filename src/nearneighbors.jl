import Base.length, Base.getindex
export NearNeighborCache, MetricNN, QuasiMetricNN, Neighborhood, ControlInfo, NullControl, NNDatastructureCache
export filter_neighborhood, addpoints, initcache, controltype, saveNN, loadNN!
export inball, knn, inball!, knn!, mutualknn!         # for metrics only (symmetric)
export inballF, knnF, inballF!, knnF!, mutualknnF!    # forwards
export inballB, knnB, inballB!, knnB!, mutualknnB!    # backwards

abstract NearNeighborCache
abstract MetricNN <: NearNeighborCache
abstract QuasiMetricNN <: NearNeighborCache

abstract ControlInfo
type NullControl <: ControlInfo end

type Neighborhood{T<:FloatingPoint}
    inds::Vector{Int}
    ds::Vector{T}
end
filter_neighborhood(n::Neighborhood, f::BitVector) = Neighborhood(n.inds[f[n.inds]], n.ds[f[n.inds]])

controltype(d::PreMetric) = NullControl
addpoints(NN::NearNeighborCache, W) = typeof(NN)([NN.V, W], NN.dist)
length(NN::NearNeighborCache) = length(NN.V)
getindex(NN::NearNeighborCache, i) = NN.V[i]

type NNDatastructureCache
    V
    DS
    US
end
NNDatastructureCache(NN::NearNeighborCache) = NNDatastructureCache(NN.V, NN.DS, NN.US)
function NNDatastructureCache(fname::String)
    open(fname, "r") do f
        deserialize(f)
    end
end
function saveNN(DC::NNDatastructureCache, fname::String)
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