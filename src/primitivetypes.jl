import Base: getindex, eltype, convert, length
import Distances: evaluate, colwise
import NearestNeighbors: inrange

### State Typedefs
export AbstractState, State, Path, AbstractBitsState, BitsState
export statevec2mat, statemat2vec, changeprecision, dense

"The `supertype` for all non-`AbstractVector` states."
abstract AbstractState
"Type union of all states, encompassing `AbstractVector`s and `AbstractState`s."
typealias State Union{AbstractVector, AbstractState}
"Shorthand for a list (`Vector`) of states."
typealias Path{S<:State} Vector{S}
"The `supertype` for all non-`SVector` states that satisfy `isbits(S) == true`."
abstract AbstractBitsState <: AbstractState
"Type union encompassing all states that may satisfy `isbits(S) == true`."
typealias BitsState Union{SVector, FieldVector, AbstractBitsState}

"Reinterprets a list of `BitsState`s as a `Matrix` without copying/re-allocating memory."    # TODO: methods for nonbits
function statevec2mat{S<:BitsState}(V::Vector{S})
    isbits(S) ? reinterpret(eltype(S), V, (length(S), length(V))) : error("$S must satisfy isbits($S) = true")
end
statevec2mat{S<:BitsState}(V::S...) = statevec2mat(collect(V))
statevec2mat{N,S<:BitsState}(V::NTuple{N,S}) = statevec2mat(collect(V))
"Reinterprets a `Matrix` as a column-wise list of `BitsState`s without copying/re-allocating memory."
function statemat2vec{S<:BitsState,T}(::Type{S}, M::Matrix{T})
    if isbits(S)
        length(S) == size(M,1) ? reinterpret(S, M, (size(M,2),)) : error("$S size mismatch with Matrix")
    else
        error("$S must satisfy isbits($S) = true")
    end
end
"Convert `x` to a similar type with numerical precision `T <: AbstractFloat` (e.g. `Float64`, `Float32`, etc.)."
changeprecision{T<:AbstractFloat,S}(::Type{T}, x::S) = convert(changeprecision(T,S), x)
changeprecision{T<:AbstractFloat,S,N}(::Type{T}, ::Type{Array{S,N}}) = Array{changeprecision(T,S),N}
changeprecision{T<:AbstractFloat,S<:Real,N}(::Type{T}, ::Type{Array{S,N}}) = Array{T,N}
changeprecision{T<:AbstractFloat,N,S<:Real}(::Type{T}, ::Type{SVector{N,S}}) = SVector{N,T}
changeprecision{T<:AbstractFloat,N1,N2,S<:Real,L}(::Type{T}, ::Type{SMatrix{N1,N2,S,L}}) = SMatrix{N1,N2,T,L}
changeprecision{T<:AbstractFloat,S<:Real}(::Type{T}, ::Type{S}) = T
changeprecision{T<:AbstractFloat,S}(::Type{T}, ::Type{S}) = S

## SE2State
export SE2State

"State consisting of a position (x,y) ∈ R^2 and an angle θ ∈ R."
immutable SE2State{T} <: FieldVector{T}
    x::T
    y::T
    θ::T
end
changeprecision{T<:AbstractFloat,S<:Real}(::Type{T}, ::Type{SE2State{S}}) = SE2State{T}


### StateSpace Typedefs
export StateSpace, State2Workspace

"The `supertype` for all state spaces; the type parameter is the state type."
abstract StateSpace{S<:State}
"Encodes a transformation from the state space (dynamics) to the workspace (obstacles)."
abstract State2Workspace


### Metrics and QuasiMetrics
export QuasiMetric, PreMetric, Metric, ChoppedMetric, ChoppedQuasiMetric, ChoppedPreMetric

"A quasimetric satisfies positivity, positive definiteness, and the triangle inequality (no symmetry)."
abstract QuasiMetric <: PreMetric

## Fallbacks (should usually be extended)
colwise{S<:State}(d::PreMetric, v::S, W::Vector{S}) = eltype(S)[evaluate(d, v, w) for w in W]
colwise{S<:State}(d::PreMetric, W::Vector{S}, v::S) = eltype(S)[evaluate(d, w, v) for w in W]

## Chopped, Lower-Bounded Metrics
"""
Evaluates as `m(v, w) <= chopval ? m(v, w) : Inf` for points `v` and `w`.\n
`lowerbound` should be easier to evaluate than `m` and must satisfy `lowerbound(v, w) ≤ m(v, w)` for all `v`, `w`.
"""
type ChoppedMetric{M<:Metric,B<:Metric,T<:AbstractFloat} <: PreMetric
    m::M
    lowerbound::B
    chopval::T
end
"""
Evaluates as `m(v, w) <= chopval ? m(v, w) : Inf` for points `v` and `w`.\n
`lowerbound` should be easier to evaluate than `m` and must satisfy `lowerbound(v, w) ≤ m(v, w)` for all `v`, `w`.
"""
type ChoppedQuasiMetric{M<:QuasiMetric,B<:Metric,T<:AbstractFloat} <: PreMetric
    m::M
    lowerbound::B
    chopval::T
end
"Type union encompassing `ChoppedMetric` and `ChoppedQuasiMetric`."
typealias ChoppedPreMetric{M,B,T} Union{ChoppedMetric{M,B,T}, ChoppedQuasiMetric{M,B,T}}
function evaluate(clbm::ChoppedPreMetric, v::State, w::State)
    lb = evaluate(clbm.lowerbound, v, w)
    lb > clbm.chopval && return oftype(clbm.chopval, Inf)   # Inf is a sentinel value more numerically robust than using
    d = evaluate(clbm.m, v, w)                              # chopval; technically breaks the metric inequality
    d <= clbm.chopval ? d : oftype(clbm.chopval, Inf)
end
for M in (:ChoppedMetric, :ChoppedQuasiMetric)
    @eval changeprecision{T<:AbstractFloat}(::Type{T}, x::$M) = $M(changeprecision(T, x.m),
                                                                   changeprecision(T, x.lowerbound),
                                                                   T(x.chopval))
end
typealias SymmetricDistance Union{Metric, ChoppedMetric}
typealias AsymmetricDistance Union{QuasiMetric, ChoppedQuasiMetric}


### Steering
export ControlInfo, NullControl, StepControl, DurationAndTargetControl
export ControlSequence, ZeroOrderHoldControl, TimestampedTrajectoryControl
export duration, splitcontrol

abstract ControlInfo
immutable NullControl <: ControlInfo end
immutable StepControl{T<:AbstractFloat,N} <: ControlInfo
    t::T
    u::SVector{N,T}
end
StepControl{T}(t::T, u::AbstractVector{T}) = StepControl(t, SVector(u))
immutable DurationAndTargetControl{T<:AbstractFloat,N} <: ControlInfo
    t::T
    x::SVector{N,T}
end
DurationAndTargetControl{T}(t::T, x::AbstractVector{T}) = DurationAndTargetControl(t, SVector(x))
typealias ControlSequence{C<:ControlInfo} Vector{C}
typealias ZeroOrderHoldControl{C<:StepControl} Vector{C}
typealias TimestampedTrajectoryControl{C<:DurationAndTargetControl} Vector{C}

duration(x::ControlInfo) = x.t
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