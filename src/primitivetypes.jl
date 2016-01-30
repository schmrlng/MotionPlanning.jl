import Base: getindex, eltype, convert, hcat, full
import Distances: evaluate, colwise
import NearestNeighbors: inrange

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
SE2State{T<:AbstractFloat}(x::Vector{T}) = SE2State(x[1], x[2], x[3])
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
full{T}(s::SE2State{T}) = T[s.x[1], s.x[2], s.t]

### Abstract State Space Typedefs
export StateSpace, State2Workspace

abstract StateSpace{T<:AbstractFloat}
abstract State2Workspace

### Metrics and QuasiMetrics
export QuasiMetric, ChoppedLowerBoundedMetric, ChoppedLowerBoundedQuasiMetric, ChoppedLowerBoundedPreMetric

abstract QuasiMetric <: PreMetric   # positivity, positive definiteness, and triangle inequality (no symmetry)
colwise{S<:State}(d::PreMetric, v::S, W::Vector{S}) = eltype(S)[evaluate(d, v, w) for w in W]
colwise{S<:State}(d::PreMetric, W::Vector{S}, v::S) = eltype(S)[evaluate(d, w, v) for w in W]
evaluate(M::Euclidean, v::Vec, w::Vec) = norm(v - w)
evaluate(M::Euclidean, v::SE2State, w::SE2State) = norm(v.x - w.x)
inrange(tree::NearestNeighbors.NNTree, v::Vec, args...) = inrange(tree, full(v), args...)
inrange(tree::NearestNeighbors.NNTree, v::SE2State, args...) = inrange(tree, full(v.x), args...)
type ChoppedLowerBoundedMetric{M<:Metric,B<:Metric,T<:AbstractFloat} <: Metric
    m::M
    lowerbound::B
    chopval::T
end
type ChoppedLowerBoundedQuasiMetric{M<:QuasiMetric,B<:Metric,T<:AbstractFloat} <: QuasiMetric
    m::M
    lowerbound::B
    chopval::T
end
typealias ChoppedLowerBoundedPreMetric{M,B,T} Union{ChoppedLowerBoundedMetric{M,B,T}, ChoppedLowerBoundedQuasiMetric{M,B,T}}
function evaluate(clbm::ChoppedLowerBoundedPreMetric, v::State, w::State)
    lb = evaluate(clbm.lowerbound, v, w)
    lb >= clbm.chopval && return clbm.chopval
    min(evaluate(clbm.m, v, w), clbm.chopval)
end
# function NearestNeighbors.interpolate(clbm::ChoppedLowerBoundedPreMetric, v::State, w::State)
#     NearestNeighbors.interpolate(clbm.m, v, w)
# end

### Steering
export ControlInfo, NullControl

abstract ControlInfo
immutable NullControl <: ControlInfo end