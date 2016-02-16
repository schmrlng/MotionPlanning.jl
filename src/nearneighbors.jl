import Base: length, getindex
export SampleSet, MetricNN, QuasiMetricNN
export NearNeighborCache, MutableNNC, ImmutableNNC
export ControlCache, EmptyControlCache, MutableControlCache, ImmutableControlCache
export DistanceDataStructure, EmptyDistanceDS, BruteDistanceDS, TreeDistanceDS
export filter_neighborhood, addpoints
# export resetcache!
# export saveNN, loadNN!
export inball, knn, inball!, knn!, mutualknn!         # for metrics only (symmetric)
export inballF, knnF, inballF!, knnF!, mutualknnF!    # forwards
export inballB, knnB, inballB!, knnB!, mutualknnB!    # backwards

### Typedefs
abstract SampleSet
abstract NearNeighborCache{T}; begin
    ## NearNeighborCache implementations
    # TODO: consider EmptyNNC for online algorithms or streaming-batch
    immutable MutableNNC{T} <: NearNeighborCache{T}
        D::Vector{SparseVector{T,Int}}
        r::Vector{T}
    end
    MutableNNC{T}(N::Int, ::Type{T}) = MutableNNC{T}(Array(SparseVector{T,Int}, N), zeros(T, N))
    immutable ImmutableNNC{T} <: NearNeighborCache{T}
        D::SparseMatrixCSC{T,Int}
        r::Vector{T}
    end
    ImmutableNNC(M::MutableNNC) = ImmutableNNC(hcat(M.D...), M.r)
    MutableNNC{T}(I::ImmutableNNC{T}) = MutableNNC(size(I.D,1), T) # clear the cache
end

abstract ControlCache{U<:ControlInfo}; begin
    ## ControlCache implementations 
    immutable EmptyControlCache{U} <: ControlCache{U} end
    EmptyControlCache{U}(::Type{U}) = EmptyControlCache{U}()
    EmptyControlCache() = EmptyControlCache(NullControl)
    immutable MutableControlCache{U} <: ControlCache{U}
        US::Vector{SparseVector{U,Int}}
    end
    MutableControlCache{U}(N::Int, ::Type{U}) = MutableControlCache{U}(Array(SparseVector{U,Int}, N))
    immutable ImmutableControlCache{U} <: ControlCache{U}
        US::SparseMatrixCSC{U,Int}
    end
    ImmutableControlCache(M::MutableControlCache) = ImmutableControlCache(hcat(M.U...))
    MutableControlCache{U}(I::ImmutableControlCache{U}) = MutableControlCache(size(I.US,1), U) # clear the cache
end

abstract DistanceDataStructure{T}; begin
    ## DistanceDataStructure implementations
    immutable EmptyDistanceDS{T} <: DistanceDataStructure{T} end
    EmptyDistanceDS{T}(::Type{T}) = EmptyDistanceDS{T}()
    immutable BruteDistanceDS{A<:AbstractMatrix,T} <: DistanceDataStructure{T}
        Dmat::A
    end
    BruteDistanceDS{T}(Dmat::AbstractMatrix{T}) = BruteDistanceDS{typeof(Dmat),T}(Dmat)
    immutable TreeDistanceDS{NNT<:NearestNeighbors.NNTree,T} <: DistanceDataStructure{T}
        tree::NNT
    end
    TreeDistanceDS{T}(tree::NearestNeighbors.NNTree{T}) = TreeDistanceDS{typeof(tree),T}(tree)
end

### Construction
type MetricNN{M<:Metric,S<:State,N<:NearNeighborCache,D<:DistanceDataStructure,U<:ControlCache} <: SampleSet
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

type QuasiMetricNN{M<:QuasiMetric,S<:State,N<:NearNeighborCache,D<:DistanceDataStructure,U<:ControlCache} <: SampleSet
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

## Fallbacks (should usually be extended)
helper_data_structures{S}(V::Vector{S}, d::Metric) = (EmptyDistanceDS(eltype(S)), EmptyControlCache(controltype(d)))
function helper_data_structures{S}(V::Vector{S}, d::QuasiMetric)
    DS, US = EmptyDistanceDS(eltype(S)), EmptyControlCache(controltype(d))
    (DS, US, DS, US)
end

### Utilities

@inline function filter_neighborhood(n::AbstractSparseVector, f::BitVector)
    inds, ds = nonzeroinds(n), nonzeros(n)
    SparseVector(length(n), inds[f[inds]], ds[f[inds]])
end
addpoints(NN::MetricNN, W) = MetricNN([NN.V; W], NN.dist, NN.init)
addpoints(NN::QuasiMetricNN, W) = QuasiMetricNN([NN.V; W], NN.dist, NN.init)
length(NN::SampleSet) = length(NN.V)
getindex(NN::SampleSet, i::Int) = i > 0 ? NN.V[i] : NN.init
getindex(NN::SampleSet, I::AbstractVector) = [NN[i] for i in I]

### Serialization



### Near Neighbor Computation

@inline inball!(NN::MetricNN, v::Int, r, f::BitVector) = filter_neighborhood(inball!(NN, v, r), f)
@inline inballF!(NN::QuasiMetricNN, v::Int, r, f::BitVector) = filter_neighborhood(inballF!(NN, v, r), f)
@inline inballB!(NN::QuasiMetricNN, v::Int, r, f::BitVector) = filter_neighborhood(inballB!(NN, v, r), f)
for direction in ("", "F", "B")
    NNsym = (direction == "" ? :MetricNN : :QuasiMetricNN)
    inballD! = symbol(:inball, direction, "!")
    cacheD = symbol(:cache, direction)
    DSD = symbol(:DS, direction)
    @eval $inballD!{M,S,N<:ImmutableNNC}(NN::$NNsym{M,S,N}, v::Int, r) = subcol(NN.$cacheD.D, v)
    @eval function $inballD!{M,S,N<:MutableNNC}(NN::$NNsym{M,S,N}, v::Int, r)
        if !isdefined(NN.$cacheD.D, v)
            NN.$cacheD.D[v] = inball(NN.V, NN.dist, NN.$DSD, v, r, $(direction != "B"))
            NN.$cacheD.r[v] = r
        end
        NN.$cacheD.D[v]
    end
end

function inball(V::Vector, dist::Metric, DS::DistanceDataStructure, v::Int, r, forwards::Bool = true)
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
    vcol = subcol(DS.dmat, v)
    inds, ds = nonzeroinds(vcol), nonzeros(vcol)
    @devec nn_bool = ds .<= r
    SparseVector(length(V), inds[nn_bool], ds[nn_bool])
end

function inball{S}(V::Vector{S}, dist::Metric, DS::TreeDistanceDS, v::Int, r, forwards::Bool = true)
    inds = inrange(DS.tree, V[v], r, true)
    inds = deleteat!(inds, searchsortedfirst(inds, v))
    SparseVector(length(V), inds, forwards ? colwise(dist, V[v], DS.tree.data[:,inds]) : colwise(dist, DS.tree.data[:,inds], V[v]))
end

function inball{S}(V::Vector{S}, dist::ChoppedPreMetric, DS::TreeDistanceDS, v::Int, r, forwards::Bool = true)
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

for ff in [x -> symbol("inball", x), x -> symbol("inball", x, "!")]
    @eval $(ff("F"))(NN::MetricNN, args...) = $(ff(""))(NN, args...)
    @eval $(ff("B"))(NN::MetricNN, args...) = $(ff(""))(NN, args...)
end
