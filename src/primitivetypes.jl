import Base: getindex, eltype, convert, hcat

### Planning Primitives
export AbstractState, State, Path

abstract AbstractState
typealias State Union{AbstractVector, FixedVector, AbstractState}
typealias Path{T<:State} Vector{T}

### State Typedefs
export SE2State

immutable SE2State{T<:AbstractFloat} <: AbstractState
    x::Vec{2,T}
    t::T
end
SE2State{T<:AbstractFloat}(x::T, y::T, t::T) = SE2State(Vec{2,T}(x,y), t)
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
convert{T,S}(::Type{Vector{T}}, s::SE2State{S}) = convert(Vector{T}, s.x)       # TODO: a dangerous default; should be eliminated
function hcat{T}(X::SE2State{T}...)                                             # TODO: also sort of weird, but very useful
    result = Array(T, 2, length(X))
    @inbounds for i in 1:length(X), j in 1:2
        result[j,i] = X[i].x[j]
    end
    result
end

### Abstract State Space Typedefs
export StateSpace, State2Workspace, QuasiMetric, ControlInfo, NullControl

abstract StateSpace
abstract State2Workspace
abstract QuasiMetric <: PreMetric   # positivity, positive definiteness, and triangle inequality (no symmetry)

abstract ControlInfo
type NullControl <: ControlInfo end