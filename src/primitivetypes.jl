import Base: getindex, eltype, convert, hcat

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
hcat(X::SE2State...) = hcat([X[i].x for i in 1:length(X)]...)             # TODO: also sort of weird, but very useful

### Abstract State Space Typedefs
export StateSpace, State2Workspace, QuasiMetric, ControlInfo, NullControl

abstract StateSpace
abstract State2Workspace
abstract QuasiMetric <: PreMetric   # positivity, positive definiteness, and triangle inequality (no symmetry)

abstract ControlInfo
type NullControl <: ControlInfo end