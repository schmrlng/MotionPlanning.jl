import Base: getindex, eltype, convert

### Planning Primitives
export AbstractState, State, Path

abstract AbstractState
typealias State Union(AbstractVector, AbstractState)
typealias Path{T<:State} Vector{T}

### State Typedefs
export SE2State

immutable SE2State{T<:FloatingPoint} <: AbstractState
    x::Vector2{T}
    t::T
end
SE2State{T<:FloatingPoint}(x::T, y::T, t::T) = SE2State(Vector2{T}(x,y), t)
function getindex(s::SE2State, d::Symbol)  # incurs a slight performance hit over direct element access
    if d == :x
        return s.x[1]
    elseif d == :y
        return s.x[2]
    elseif d == :t
        return s.t
    end
    throw(KeyError(d))
end
eltype{T}(::SE2State{T}) = T
eltype{T}(::Type{SE2State{T}}) = T
convert{T}(::Type{Vector{T}}, s::SE2State{T}) = convert(Vector{T}, s.x)  # TODO: a dangerous default; should be eliminated