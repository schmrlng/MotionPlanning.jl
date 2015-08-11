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