import Base: getindex, eltype, convert, length, hcat
import Distances: evaluate, colwise
import NearestNeighbors: inrange

### State Typedefs
export AbstractState, State, Path, AbstractBitsState, BitsState
export statevec2mat, statemat2vec, changeprecision

"The `supertype` for all non-vector states."
abstract AbstractState
"Type union of all states, including user-defined (`AbstractState`) and mutable/immutable vectors."
typealias State Union{AbstractVector, FixedVector, AbstractState}
"Shorthand for a list (`Vector`) of states."
typealias Path{S<:State} Vector{S}
"The `supertype` for all non-`FixedVector` states that satisfy `isbits(S) == true`."
abstract AbstractBitsState <: AbstractState
"Type union encompassing all states that satisfy `isbits(S) == true`."
typealias BitsState Union{FixedVector, AbstractBitsState}

"Reinterprets a list of `BitsState`s as a `Matrix` without copying/re-allocating memory."
statevec2mat{S<:BitsState}(V::Vector{S}) =
    isbits(S) ? reinterpret(eltype(S), V, (length(S), length(V))) : error("eltype of BitsState may not be a bitstype?")
"Reinterprets a `Matrix` as a column-wise list of `BitsState`s."
statemat2vec{S<:BitsState,T}(::Type{S}, M::Matrix{T}) =
    isbits(S) && length(S) == size(M,1) ? reinterpret(S, M, (size(M,2),)) : error("State/Matrix size mismatch?")
"Convert a `State` to a `Vector`."
dense{S<:State}(x::S) = convert(Vector{eltype(S)}, x)
"Convert `x` to a similar type with numerical precision `T <: AbstractFloat` (e.g. `Float64`, `Float32`, etc.)."
changeprecision{T<:AbstractFloat,S<:State}(::Type{T}, x::S) = convert(changeprecision(T,S), x)
changeprecision{T<:AbstractFloat}(::Type{T}, x::Vector) = [changeprecision(T,i) for i in x]
changeprecision{T<:AbstractFloat}(::Type{T}, x::AbstractFloat) = T(x)

## SE2State
export SE2State

"State consisting of a position x ∈ R^2 and an angle θ ∈ R."
immutable SE2State{T<:AbstractFloat} <: AbstractBitsState
    x::Vec{2,T}
    θ::T
end
SE2State{T<:AbstractFloat}(x::T, y::T, θ::T) = SE2State(Vec{2,T}(x,y), θ)
SE2State{T<:AbstractFloat}(x::Vector{T}) = SE2State(x[1], x[2], x[3])
function getindex(s::SE2State, d::Symbol)  # incurs a slight performance hit over direct element access
    if d == :x
        return s.x[1]
    elseif d == :y
        return s.x[2]
    elseif d == :θ || d == :t
        return s.θ
    end
    throw(KeyError(d))
end
eltype{T}(::SE2State{T}) = T
eltype{T}(::Type{SE2State{T}}) = T
length{T}(::SE2State{T}) = 3
length{T}(::Type{SE2State{T}}) = 3
convert{T}(::Type{Vector{T}}, s::SE2State) = T[s.x[1], s.x[2], s.θ]
convert{T}(::Type{SE2State{T}}, s::SE2State) = SE2State(T(s.x[1]), T(s.x[2]), T(s.θ))
convert{T}(::Type{SE2State{T}}, x::Vector) = SE2State(T(x[1]), T(x[2]), T(x[3]))
changeprecision{T<:AbstractFloat,S}(::Type{T}, ::Type{SE2State{S}}) = SE2State{T}

## Vec (defined in FixedSizeArrays)
export Vec

function hcat{N,T}(X::Vec{N,T}...)
    result = Array(T, N, length(X))
    @inbounds for i in 1:length(X), j in 1:N
        result[j,i] = X[i][j]
    end
    result
end
changeprecision{N,T<:AbstractFloat,S}(::Type{T}, ::Type{Vec{N,S}}) = Vec{N,T}

## Vector
changeprecision{T<:AbstractFloat,S<:AbstractFloat}(::Type{T}, ::Type{Vector{S}}) = Vector{T}

### StateSpace Typedefs
export StateSpace, State2Workspace

"The `supertype` for all state spaces."
abstract StateSpace{T<:AbstractFloat}
"Encodes a transformation from the state space (dynamics) to the workspace (obstacles)."
abstract State2Workspace

### Metrics and QuasiMetrics
export QuasiMetric, ChoppedMetric, ChoppedQuasiMetric, ChoppedPreMetric

"A quasimetric satisfies positivity, positive definiteness, and the triangle inequality (no symmetry)."
abstract QuasiMetric <: PreMetric

## Fallbacks (should usually be extended)
colwise{S<:State}(d::PreMetric, v::S, W::Vector{S}) = eltype(S)[evaluate(d, v, w) for w in W]
colwise{S<:State}(d::PreMetric, W::Vector{S}, v::S) = eltype(S)[evaluate(d, w, v) for w in W]
changeprecision{T<:AbstractFloat}(::Type{T}, x::PreMetric) = x

## Chopped, Lower-Bounded Metrics
"""
Evaluates as `min(m(v, w), chopval)` for points `v` and `w`.\n
`lowerbound` should be easier to evaluate than `m` and must satisfy `lowerbound(v, w) ≤ m(v, w)` for all `v`, `w`.
"""
type ChoppedMetric{M<:Metric,B<:Metric,T<:AbstractFloat} <: Metric
    m::M
    lowerbound::B
    chopval::T
end
"""
Evaluates as `min(m(v, w), chopval)` for points `v` and `w`.\n
`lowerbound` should be easier to evaluate than `m` and must satisfy `lowerbound(v, w) ≤ m(v, w)` for all `v`, `w`.
"""
type ChoppedQuasiMetric{M<:QuasiMetric,B<:Metric,T<:AbstractFloat} <: QuasiMetric
    m::M
    lowerbound::B
    chopval::T
end
"Type union encompassing `ChoppedMetric` and `ChoppedQuasiMetric`."
typealias ChoppedPreMetric{M,B,T} Union{ChoppedMetric{M,B,T}, ChoppedQuasiMetric{M,B,T}}
function evaluate(clbm::ChoppedPreMetric, v::State, w::State)
    lb = evaluate(clbm.lowerbound, v, w)
    lb >= clbm.chopval && return oftype(clbm.chopval, Inf)
    d = evaluate(clbm.m, v, w)
    d <= clbm.chopval ? d : oftype(clbm.chopval, Inf)
end
for M in (:ChoppedMetric, :ChoppedQuasiMetric)
    @eval changeprecision{T<:AbstractFloat}(::Type{T}, x::$M) = $M(changeprecision(T, x.m),
                                                                   changeprecision(T, x.lowerbound),
                                                                   T(x.chopval))
end

### Steering
export ControlInfo, NullControl, StepControl, DurationAndTargetControl
export ControlSequence, ZeroOrderHoldControl, TimestampedTrajectoryControl
export duration, splitcontrol

abstract ControlInfo
immutable NullControl <: ControlInfo end
immutable StepControl{T<:AbstractFloat,N} <: ControlInfo
    t::T
    u::Vec{N,T}
end
immutable DurationAndTargetControl{T<:AbstractFloat,N} <: ControlInfo
    t::T
    x::Vec{N,T}
end
typealias ControlSequence{C<:ControlInfo} Vector{C}
typealias ZeroOrderHoldControl{C<:StepControl} Vector{C}
typealias TimestampedTrajectoryControl{C<:DurationAndTargetControl} Vector{C}

duration(x::StepControl) = x.t
duration(x::DurationAndTargetControl) = x.t
duration(x::ControlSequence) = sum(duration(s) for s in x)
function splitcontrol{T}(x::StepControl{T}, t::Real)
    if t < 0
        StepControl(T(0), x.u), x
    elseif t < x.t
        StepControl(T(t), x.u), StepControl(x.t - T(t), x.u)
    else
        x, StepControl(T(0), x.u)
    end
end
function splitcontrol{T,R<:Real}(x::StepControl{T}, t::AbstractVector{R})
    issorted(t) || error("Times should be sorted as input to splitcontrol.")
    d = duration(x)
    push!([StepControl(T(dt), x.u) for dt in diff(clamp(t, T(0), d))], StepControl(d-T(t[end]), x.u))
end
function splitcontrol{T<:AbstractFloat}(x0::ZeroOrderHoldControl, t::Union{T, AbstractVector{T}})
    issorted(t) || error("Times should be sorted as input to splitcontrol.")
    x = copy(x0)    # TODO: not the most memory efficient, but simple for now
    t = clamp(t, T(0), duration(x))
    splitu = [similar(x, 0) for i in 1:length(t)+1]
    d = T(0)
    uidx = 1
    for (i, ti) in enumerate(t)
        while ti > d + duration(x[uidx]) && uidx < length(x0)    # second clause necessary because of numerical error
            push!(splitu[i], x[uidx])
            d += duration(x[uidx])
            uidx += 1
        end
        u1, u2 = splitcontrol(x[uidx], ti - d)
        push!(splitu[i], u1)
        x[uidx] = u2
        d = ti
    end
    splitu[end] = x[uidx:end]
    splitu
end