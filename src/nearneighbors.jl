import Base.length, Base.getindex
export NearNeighborCache, MetricNN, QuasiMetricNN, Neighborhood, NNDatastructureCache
export filter_neighborhood, addpoints, initcache, saveNN, loadNN!
export inball, knn, inball!, knn!, mutualknn!         # for metrics only (symmetric)
export inballF, knnF, inballF!, knnF!, mutualknnF!    # forwards
export inballB, knnB, inballB!, knnB!, mutualknnB!    # backwards

abstract SampleSet
abstract NearNeighborCache{T}
    immutable MutableNNC{T} <: NearNeighborCache{T}
        D::Vector{SparseVector{T,Int}}
        r::Vector{T}
    end
    MutableNNC(N::Int, T::DataType) = MutableNNC{T}(Array(SparseVector{T,Int}, N), zeros(T, N))
    immutable ImmutableNNC{T} <: NearNeighborCache{T}
        D::SparseMatrixCSC{T,Int}
        r::Vector{T}
    end

abstract DistanceDataStructure{T}
    immutable EmptyDistanceDS{T} <: DistanceDataStructure{T} end
    EmptyDistanceDS(T::DataType) = EmptyDistanceDS{T}()
    immutable BruteDistanceDS{A<:AbstractMatrix,T} <: DistanceDataStructure{T}
        Dmat::A
    end
    BruteDistanceDS{T}(Dmat::AbstractMatrix{T}) = BruteDistanceDS{typeof(Dmat),T}(Dmat)
    immutable TreeDistanceDS{NNT<:NearestNeighbors.NNTree,T} <: DistanceDataStructure{T}
        tree::NNT
    end
    TreeDistanceDS{T}(tree::NearestNeighbors.NNTree{T}) = TreeDistanceDS{typeof(tree),T}(tree)

abstract ControlDataStructure{U}
    immutable EmptyControlDS{U} <: ControlDataStructure{U} end
    EmptyControlDS(U::DataType) = EmptyControlDS{U}()
    EmptyControlDS() = EmptyControlDS{NullControl}()
    immutable BruteControlDS{U,A<:AbstractMatrix} <: ControlDataStructure{U}
        Umat::A
    end
    BruteControlDS{U}(Umat::AbstractMatrix{U}) = BruteControlDS{U,typeof(Umat)}(Umat)

helper_data_structures{S}(V::Vector{S}, dist::Metric) = (EmptyDistanceDS(eltype(S)), EmptyControlDS(controltype(dist)))

type MetricNN{M<:Metric,S<:State,N<:NearNeighborCache,D<:DistanceDataStructure,U<:ControlDataStructure} <: SampleSet
    V::Vector{S}
    dist::M
    init::S
    cache::N
    DS::D
    US::U
end
function MetricNN{S<:State}(V::Vector{S}, dist::Metric, init::S)
    N = length(V)
    T = eltype(S)
    MetricNN(V, dist, init, MutableNNC(N,T), helper_data_structures(V, dist)...)
end

type QuasiMetricNN{M<:QuasiMetric,S<:State,N<:NearNeighborCache,D<:DistanceDataStructure,U<:ControlDataStructure} <: SampleSet
    V::Vector{S}
    dist::M
    init::S
    cacheF::N   # Everything's just pointers, but F and B may be different (e.g. CSC and CSR), trading memory for speed
    cacheB::N
    DSF::D
    USF::U
    DSB::D
    USB::U
end
function QuasiMetricNN{S<:State}(V::Vector{S}, dist::QuasiMetric, init::S)
    N = length(V)
    T = eltype(S)
    QuasiMetricNN(V, dist, init, MutableNNC(N,T), MutableNNC(N,T), helper_data_structures(V, dist)...)
end

@inline function filter_neighborhood(n::AbstractSparseVector, f::BitVector)
    inds, ds = nonzeroinds(n), nonzeros(n)
    SparseVector(length(n), inds[f[inds]], ds[f[inds]])
end
addpoints(NN::MetricNN, W) = MetricNN([NN.V; W], NN.dist, NN.init)
addpoints(NN::QuasiMetricNN, W) = QuasiMetricNN([NN.V; W], NN.dist, NN.init)
length(NN::SampleSet) = length(NN.V)
getindex(NN::SampleSet, i::Int) = i > 0 ? NN.V[i] : NN.init
getindex(NN::SampleSet, I::AbstractVector) = [NN[i] for i in I]

## Computation of near neighbors

@inline inball!(NN::MetricNN, v::Int, r, f::BitVector) = filter_neighborhood(inball!(NN, v, r), f)
@inline inballF!(NN::QuasiMetricNN, v::Int, r, f::BitVector) = filter_neighborhood(inballF!(NN, v, r), f)
@inline inballB!(NN::QuasiMetricNN, v::Int, r, f::BitVector) = filter_neighborhood(inballB!(NN, v, r), f)
inball!{M,S,N<:ImmutableNNC}(NN::MetricNN{M,S,N}, v::Int, r) = getcol(NN.cache.D, v)
inballF!{M,S,N<:ImmutableNNC}(NN::QuasiMetricNN{M,S,N}, v::Int, r) = getcol(NN.cacheF.D, v)
inballB!{M,S,N<:ImmutableNNC}(NN::QuasiMetricNN{M,S,N}, v::Int, r) = getcol(NN.cacheB.D, v)
function inball!{M,S,N<:MutableNNC}(NN::MetricNN{M,S,N}, v::Int, r)
    if !isdefined(NN.cache.D, v)
        NN.cache.D[v] = inball(NN.V, NN.dist, NN.DS, v, r)
        NN.cache.r[v] = r
    end
    NN.cache.D[v]
end
function inballF!{M,S,N<:MutableNNC}(NN::QuasiMetricNN{M,S,N}, v::Int, r)
    if !isdefined(NN.cacheF.D, v)
        NN.cacheF.D[v] = inball(NN.V, NN.dist, NN.DSF, v, r)
        NN.cacheF.r[v] = r
    end
    NN.cacheF.D[v]
end
function inballB!{M,S,N<:MutableNNC}(NN::QuasiMetricNN{M,S,N}, v::Int, r)
    if !isdefined(NN.cacheB.D, v)
        NN.cacheB.D[v] = inball(NN.V, NN.dist, NN.DSB, v, r, false)
        NN.cacheB.r[v] = r
    end
    NN.cacheB.D[v]
end

function inball(V::Vector, dist::Metric, DS::EmptyDistanceDS, v::Int, r, forwards::Bool = true)
    ds = forwards ? colwise(dist, V[v], V) : colwise(dist, V, V[v])
    @devec nn_bool = ds .<= r
    nn_bool[v] = false
    inds = find(nn_bool)
    SparseVector(length(V), inds, ds[inds])
end

function inball{A<:Matrix}(V::Vector, dist::Metric, DS::BruteDistanceDS{A}, v::Int, r, forwards::Bool = true)
    @devec nn_bool = DS.dmat[:,v] .<= r
    nn_bool[v] = false
    inds = find(nn_bool)
    SparseVector(length(V), inds, DS.dmat[inds,v])
end

function inball{A<:SparseMatrixCSC}(V::Vector, dist::Metric, DS::BruteDistanceDS{A}, v::Int, r, forwards::Bool = true)
    vcol = getcol(DS.dmat, v)
    inds, ds = nonzeroinds(vcol), nonzeros(vcol)
    @devec nn_bool = ds .<= r
    SparseVector(length(V), inds[nn_bool], ds[nn_bool])
end

function inball{S}(V::Vector{S}, dist::Metric, DS::TreeDistanceDS, v::Int, r, forwards::Bool = true)
    inds = inrange(DS.tree, V[v], r, true)
    inds = deleteat!(inds, searchsortedfirst(inds, v))
    SparseVector(length(V), inds, forwards ? colwise(dist, V[v], DS.tree.data[:,inds]) : colwise(dist, DS.tree.data[:,inds], V[v]))
end

function inball{S}(V::Vector{S}, dist::ChoppedLowerBoundedPreMetric, DS::TreeDistanceDS, v::Int, r, forwards::Bool = true)
    inds = inrange(DS.tree, V[v], r, true)
    inds = deleteat!(inds, searchsortedfirst(inds, v))
    if S <: Vec
        ds = forwards ? colwise(dist, V[v], DS.tree.data[:,inds]) : colwise(dist, DS.tree.data[:,inds], V[v])
    elseif S <: SE2State
        ds = forwards ? colwise(dist, V[v], V[inds]) : colwise(dist, V[inds], V[v])
    else
        error("TODO: any other state types?")
    end
    @devec nn_bool = ds .<= r
    SparseVector(length(V), inds[nn_bool], ds[nn_bool])
end

# function inball(V::Vector, dist::Metric, DS::DistanceDataStructure, v::State, r)
#     ds = colwise(dist, v, V)
#     @devec nn_bool = ds .<= r
#     inds = find(nn_bool)
#     SparseVector(length(V), inds, ds[inds])
# end

# function inball(V::Vector, dist::Metric, DS::TreeDistanceDS, v::State, r)
#     inds = inrange(DS.tree, full(V[v]), r, true)
#     KDTrees.inball(NN.DS, convert(Vector{T}, v), r, true)
#     Neighborhood(inds, colwise(Euclidean(), v, hcat(NN[inds]...)))
# end

for ff in [x -> symbol("inball", x), x -> symbol("inball", x, "!")]
    @eval $(ff("F"))(NN::MetricNN, args...) = $(ff(""))(NN, args...)
    @eval $(ff("B"))(NN::MetricNN, args...) = $(ff(""))(NN, args...)
end



# abstract NearNeighborCache
# abstract MetricNN <: NearNeighborCache
# abstract QuasiMetricNN <: NearNeighborCache

# type Neighborhood{T<:AbstractFloat}
#     inds::Vector{Int}
#     ds::Vector{T}
# end
# filter_neighborhood(n::Neighborhood, f::BitVector) = Neighborhood(n.inds[f[n.inds]], n.ds[f[n.inds]])  # TODO: this is slow (for loop to get rid of double-scan, or perhaps best, use an iterator)

# addpoints(NN::NearNeighborCache, W) = typeof(NN)([NN.V; W], NN.dist, NN.init)
# length(NN::NearNeighborCache) = length(NN.V)
# getindex(NN::NearNeighborCache, i::Int) = i > 0 ? NN.V[i] : NN.init
# getindex(NN::NearNeighborCache, I::AbstractVector) = [NN[i] for i in I]

# type NNDatastructureCache
#     V
#     DS
#     US
# end
# NNDatastructureCache(NN::NearNeighborCache) = NNDatastructureCache(NN.V, NN.DS, NN.US)
# function NNDatastructureCache(fname::AbstractString)
#     open(fname, "r") do f
#         deserialize(f)
#     end
# end
# function saveNN(DC::NNDatastructureCache, fname::AbstractString)
#     open(fname, "w") do f
#         serialize(f, DC)
#     end
# end
# saveNN(NN::NearNeighborCache) = saveNN(NNDatastructureCache(NN))

# include("nearneighbors/metricNN.jl")
# include("nearneighbors/quasimetricNN.jl")

# # TODO: addpoints! function
# # TODO: precompute_all function
# # TODO: clean up huge amounts of copy/paste code with metaprogramming