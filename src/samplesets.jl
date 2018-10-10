export TiledIndex, SampleSet, ExplicitSampleSet, TiledSampleSet
export zero_lattice_component, shift_lattice_component, indextype

# Sample Set Indices
struct TiledIndex{D}
    on_lattice::SVector{D,Int}
    off_lattice::Int
end
TiledIndex{D}(i::Int) where {D} = TiledIndex(zero(SVector{D,Int}), i)
Base.zero(::Type{TiledIndex{D}}) where {D} = TiledIndex{D}(0)
function Base.hash(x::TiledIndex, h::UInt=zero(UInt))
    hash(reduce((r, x) -> (1 << 19 - 1)*r + x, x.on_lattice .% UInt), x.off_lattice + h)
    # hash((x.on_lattice.data..., x.off_lattice), h)
    # reduce((r, x) -> 31*r + x, SVector(h, (x.on_lattice.data .% UInt)..., x.off_lattice % UInt))
end
zero_lattice_component(i::TiledIndex{D}) where {D} = TiledIndex(zero(SVector{D,Int}), i.off_lattice)
shift_lattice_component(i::TiledIndex, shift) = TiledIndex(i.on_lattice + shift, i.off_lattice)
const SampleIndex = Union{Int,TiledIndex}

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
# node_info_datastructure(nodes::ExplicitSampleSet, f) = [f(nodes[i]) for i in eachindex(nodes)]
@inline function node_info_datastructure(nodes::ExplicitSampleSet, ::Type{T}, f) where {T}
    DefaultDict{indextype(nodes),T}(i -> f(nodes[i]), passkey=true)
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
@inline function node_info_datastructure(nodes::TiledSampleSet, ::Type{T}, f) where {T}
    DefaultDict{indextype(nodes),T}(i -> f(nodes[i]), passkey=true)
end
