export TiledIndex, SampleSet, ExplicitSampleSet, TiledSampleSet
export zero_lattice_component, shift_lattice_component, indextype

# Sample Set Indices
struct TiledIndex{D}
    on_lattice::SVector{D,Int}
    off_lattice::Int
end
TiledIndex{D}(i::Int) where {D} = TiledIndex(zero(SVector{D,Int}), i)
(::Type{<:TiledIndex})(on, off) = (D = length(on); TiledIndex{D}(convert(SVector{D,Int}, on), off))
(::Type{<:TiledIndex})(tup::NTuple{N,Int}) where {N} = TiledIndex{N-1}(reverse(Base.tail(reverse(tup))), last(tup))
(::Type{<:TiledIndex})(ci::CartesianIndex) = TiledIndex(Tuple(ci))
Base.zero(::Type{TiledIndex{D}}) where {D} = TiledIndex{D}(0)
@inline Base.iterate(i::TiledIndex) = (i, nothing)
@inline Base.iterate(i::TiledIndex, ::Any) = nothing
@inline Base.eltype(::Type{TI}) where {TI<:TiledIndex} = TI
@inline Base.length(::TiledIndex) = 1
function Base.hash(x::TiledIndex, h::UInt=zero(UInt))
    hash(reduce((r, x) -> (1 << 19 - 1)*r + x, x.on_lattice .% UInt), x.off_lattice + h)
    # hash((x.on_lattice.data..., x.off_lattice), h)
    # reduce((r, x) -> 31*r + x, SVector(h, (x.on_lattice.data .% UInt)..., x.off_lattice % UInt))
end
zero_lattice_component(i::TiledIndex{D}) where {D} = TiledIndex(zero(SVector{D,Int}), i.off_lattice)
shift_lattice_component(i::TiledIndex, shift) = TiledIndex(i.on_lattice + shift, i.off_lattice)

const SampleIndex = Union{Int,TiledIndex}
@inline linear_index(i::Int) = i
@inline linear_index(i::TiledIndex) = i.off_lattice

# Sample Sets
abstract type SampleSet end
## ExplicitSampleSet
mutable struct ExplicitSampleSet{S<:State} <: SampleSet
    init::S
    goal_samples::Vector{S}
    V::Vector{S}
end
Base.eltype(nodes::ExplicitSampleSet{S}) where {S} = S
Base.length(nodes::ExplicitSampleSet) = length(nodes.V)
Base.getindex(nodes::ExplicitSampleSet, i::Int) = i > 0 ? nodes.V[i] : (i == 0 ? nodes.init : nodes.goal_samples[-i])
Base.eachindex(nodes::ExplicitSampleSet) = eachindex(nodes.V)
indextype(nodes::ExplicitSampleSet) = Int
function addstates!(nodes::ExplicitSampleSet, x::State)
    n = length(nodes.V)
    view(push!(nodes.V, x), n+1:n+1)
end
function addstates!(nodes::ExplicitSampleSet, xs::AbstractVector{<:State})
    n = length(nodes.V)
    view(append!(nodes.V, xs), n+1:n+length(xs))
end

## TiledSampleSet
mutable struct TiledSampleSet{S<:State,M} <: SampleSet
    init::S
    goal_samples::Vector{S}
    A::M
    b::S
    Voff::Vector{S}
end
Base.eltype(nodes::TiledSampleSet{S}) where {S} = S
Base.length(nodes::TiledSampleSet) = length(nodes.Voff)
function Base.getindex(nodes::TiledSampleSet, i::TiledIndex)
    if i.off_lattice > 0
        nodes.A*i.on_lattice + nodes.b + nodes.Voff[i.off_lattice]
    elseif i.off_lattice == 0
        nodes.init
    else
        nodes.goal_samples[-i.off_lattice]
    end
end
function Base.getindex(nodes::TiledSampleSet, i::Int)
    @assert i <= 0 "::Int TiledSampleSet indices must be <= 0 (corresponding to init or goal_samples)"
    i == 0 ? nodes.init : nodes.goal_samples[-i]
end
indextype(nodes::TiledSampleSet) = TiledIndex{size(nodes.A, 2)}

# Node Info Datastructures
struct LazyNodeInfoArray{T,A<:AbstractArray{Union{Missing,T}},F}
    main::A
    extras::Vector{Union{Missing,T}}
    initializer::F
end
indextype(::LazyNodeInfoArray{T,<:AbstractArray{Union{Missing,T},1}}) where {T} = Int
indextype(::LazyNodeInfoArray{T,<:AbstractArray{Union{Missing,T},N}}) where {T,N} = TiledIndex{N-1}
_flatten(i::Int) = (i,)
_flatten(i::TiledIndex) = (Tuple(i.on_lattice)..., i.off_lattice)
function _getindex(info::LazyNodeInfoArray, i)
    linear_index(i) <= 0 ? info.extras[-linear_index(i) + 1] : info.main[_flatten(i)...]
end
function Base.checkbounds(::Type{Bool}, info::LazyNodeInfoArray, i)
    linear_index(i) <= 0 ? linear_index(i) >= 1 - length(info.extras) : checkbounds(Bool, info.main, _flatten(i)...)
end
function Base.setindex!(info::LazyNodeInfoArray, x, i)
    linear_index(i) <= 0 ? setindex!(info.extras, x, -linear_index(i) + 1) : setindex!(info.main, x, _flatten(i)...)
end
function Base.getindex(info::LazyNodeInfoArray, i)
    !checkbounds(Bool, info, i) && return missing
    res = _getindex(info, i)
    res !== missing && return res    # infers better than ismissing (JuliaLang/julia#27681)
    info[i] = info.initializer(i)
end
function Base.push!(info::LazyNodeInfoArray, x)
    push!(info.main, x)
end
function Base.keys(info::LazyNodeInfoArray)
    X = indextype(info)
    Iterators.filter(i -> _getindex(info, i) !== missing,
                     Iterators.flatten((X.(keys(info.main)), X.(1-length(info.extras):0))))
end

@inline function node_info_datastructure(nodes::ExplicitSampleSet{S}, ::Type{T}, f=(x -> T());
                                         bounds=AxisAlignedBox(fill(-Inf, S), fill(Inf, S))) where {S,T}
    LazyNodeInfoArray(Vector{Union{Missing,T}}(missing, length(nodes)),
                      Vector{Union{Missing,T}}(missing, length(nodes.goal_samples) + 1),
                      i -> f(nodes[i]))
end
@inline function node_info_datastructure(nodes::TiledSampleSet{S}, ::Type{T}, f=(x -> T());
                                         bounds=AxisAlignedBox(fill(-Inf, S), fill(Inf, S))) where {S,T}
    if any(isinf.(bounds.lo)) || any(isinf.(bounds.hi))
        DefaultDict{indextype(nodes),T}(i -> f(nodes[i]), passkey=true)
    else
        surrounding_range((x, y)) = floor(Int, x):ceil(Int, y)
        update_offsets((x, y), z) = SVector(min(-z, x), max(-z, y))
        lo, hi = bounds.lo, bounds.hi
        A, b, Voff = nodes.A, nodes.b, nodes.Voff
        M, N = size(A)
        Q, R = qr(A)
        P = R\Q'
        Voff_idx_offsets = reduce((s, x) -> update_offsets.(s, x),
                                  (P*v for v in Voff),
                                  init=fill(SVector{2,eltype(P)}(Inf, -Inf), SVector{N}))
        inds = surrounding_range.(sum(SVector.(minmax.(P.*lo', P.*hi')), dims=2) .+ Voff_idx_offsets .- P*b)
        if prod(length.(inds)) > 1e7
            DefaultDict{indextype(nodes),T}(i -> f(nodes[i]), passkey=true)
        else
            LazyNodeInfoArray(OffsetArray{Union{Missing,T}}(missing, Tuple(inds)..., 1:length(Voff)),
                              Vector{Union{Missing,T}}(missing, length(nodes.goal_samples) + 1),
                              i -> f(nodes[i]))
        end
    end
end
