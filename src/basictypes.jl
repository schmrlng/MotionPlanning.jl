export StateSpace, BoundedStateSpace
export AbstractConfig, Config
export vecs2mat, mat2vecs

# StateSpace
abstract type StateSpace{S<:State} end
## BoundedStateSpace
struct BoundedStateSpace{S,B} <: StateSpace{S}
    bounded_set::B
end
BoundedStateSpace{S}(bounded_set::B) where {S,B} = BoundedStateSpace{S,B}(bounded_set)
BoundedStateSpace(bounded_set::B) where {B} = BoundedStateSpace{typeof(rand(bounded_set)),B}(bounded_set)
BoundedStateSpace(lo::S, hi::S) where {S<:StaticVector} = BoundedStateSpace{S}(AxisAlignedBox(lo, hi))
BoundedStateSpace(lo, hi) = BoundedStateSpace(AxisAlignedBox(lo, hi))
Base.in(x, SS::BoundedStateSpace) = x in SS.bounded_set
Base.rand(SS::BoundedStateSpace{S}) where {S} = convert(S, rand(SS.bounded_set))
boundingbox(SS::BoundedStateSpace) = boundingbox(SS.bounded_set)
@recipe function f(SS::BoundedStateSpace; dims=(1, 2)) #, statespace_color=:black, statespace_linecolor=statespace_color)
    dims      --> dims
    linecolor :=  :match
    fillalpha :=  0
    label     --> ""
    SS.bounded_set
end

# Config
abstract type AbstractConfig end
const Config = Union{AbstractConfig, AbstractVector{<:Number}}

# Matrix <--> Vector of (Column) Vectors Conversion
vecs2mat(V::AbstractVector{S}) where {S<:StaticArray}   = reshape(reinterpret(eltype(S), V), (length(S), length(V)))
vecs2mat(V::AbstractVector{S}) where {S<:AbstractArray} = reduce(hcat, V)
mat2vecs(::Type{S}, M::Matrix) where {S<:StaticArray}   = reshape(reinterpret(S, M), size(M, 2))
mat2vecs(::Type{S}, M::Matrix) where {S<:AbstractArray} = [convert(S, M[:,i]) for i in 1:size(M,2)]
