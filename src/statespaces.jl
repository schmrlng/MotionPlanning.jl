export RealVectorStateSpace, SE2StateSpace
export Identity, VectorView, First2Vector2, Select2Vector2, OutputMatrix, ExtractVector
export sample_space, state2workspace, volume, dim, propagate, collision_waypoints, waypoints
export in_state_space, is_free_state, is_free_motion, is_free_path, defaultNN, setup_steering, controltype

#====================== Implementation Notes ======================
TODO: the information below is outdated.

To specify a new state space (with associated distance function M,
a Metric or QuasiMetric) the user must specify:
1. evaluate(M, v, w)
2. waypoints(M, v, w)
3. time_waypoint(M, v, w, t)    (if time-based discretization is desired)
4. cost_waypoint(M, v, w, c)    (if cost-based discretization is desired)

Optionally, to improve performance, the user may specify: ## TODO: OUTDATED BELOW
1. defaultNN(M, init)
2. workspace_waypoints(v, w, M, s2w)
3. collision_waypoints(v, w, M, s2w)    (i.e. a sparser version of workspace_waypoints)

In rare circumstances, the user might need to specify:
1. setup_steering(SS, r)
2. controltype(M)

==================================================================#

### State Space Definitions
immutable RealVectorStateSpace{N, T, M<:PreMetric, W<:State2Workspace} <: StateSpace{T}
    lo::Vec{N,T}
    hi::Vec{N,T}
    dist::M
    s2w::W
end

immutable SE2StateSpace{T, M<:PreMetric, W<:State2Workspace} <: StateSpace{T}
    lo::Vec{2,T}    # workspace only; 0 <= theta <= 2pi is assumed
    hi::Vec{2,T}
    dist::M
    s2w::W
end

### Utilities
changeprecision{T<:AbstractFloat}(::Type{T}, SS::RealVectorStateSpace) =
    RealVectorStateSpace(map(x -> changeprecision(T,x), (SS.lo, SS.hi, SS.dist, SS.s2w))...)
changeprecision{T<:AbstractFloat}(::Type{T}, SS::SE2StateSpace) =
    SE2StateSpace(map(x -> changeprecision(T,x), (SS.lo, SS.hi, SS.dist, SS.s2w))...)
changeprecision{T<:AbstractFloat}(::Type{T}, s2w::State2Workspace) = s2w

### Sampling
sample_space(SS::RealVectorStateSpace) = (SS.lo + rand(typeof(SS.lo)).*(SS.hi-SS.lo))
sample_space{T,M,W}(SS::SE2StateSpace{T,M,W}) = SE2State(SS.lo + rand(Vec{2,T}).*(SS.hi-SS.lo), convert(T, 2pi*rand(T)))
function volume(SS::StateSpace)
    # !(isa(SS, RealVectorStateSpace) && isa(SS.dist, Euclidean)) && warn("Volume not yet implemented for non-Euclidean metrics!")
    prod(SS.hi-SS.lo)
end
dim{N}(SS::RealVectorStateSpace{N}) = N # TODO: this should really be a function of the metric dist
dim(SS::SE2StateSpace) = 3

### State2Workspace (configuration space -> collision checking workspace)
immutable Identity <: State2Workspace end
immutable VectorView <: State2Workspace
    r::UnitRange{Int}
end
immutable First2Vector2 <: State2Workspace end
immutable Select2Vector2 <: State2Workspace
    a::Int
    b::Int
end
immutable OutputMatrix{M,N,T<:AbstractFloat} <: State2Workspace
    C::Mat{M,N,T}
end
OutputMatrix(C::Matrix) = OutputMatrix(Mat(C))
immutable ExtractVector <: State2Workspace end
changeprecision{T<:AbstractFloat}(::Type{T}, s2w::OutputMatrix) = OutputMatrix(changeprecision(T, s2w.C))

@unfix state2workspace(v::State, SS::StateSpace) = state2workspace(v, SS.s2w)
@unfix state2workspace(v::State, s2w::Identity) = v
@unfix state2workspace(v::Vec, s2w::VectorView) = Vec(v[s2w.r])
@unfix state2workspace(v::Vec, s2w::First2Vector2) = Vec(v[1], v[2])
@unfix state2workspace(v::Vec, s2w::Select2Vector2) = Vec(v[s2w.a], v[s2w.b])
@unfix state2workspace{M,N}(v::Vec{N}, s2w::OutputMatrix{M,N}) = s2w.C*v
@unfix state2workspace(v::SE2State, s2w::ExtractVector) = v.x

@unfix workspace2state(w::Vec, SS::StateSpace) = workspace2state(w, SS, SS.s2w)
@unfix workspace2state(w::Vec, SS::StateSpace, s2w::Identity) = w
@unfix workspace2state(w::Vec, SS::StateSpace, s2w::State2Workspace) = workspace2state(w, sample_space(SS), s2w)
@unfix function workspace2state(w::Vec, v::Vec, s2w::VectorView)
    for i in 1:length(s2w.r)
        v = setindex(v, w[i], s2w.r[i])
    end
    v
end
@unfix workspace2state(w::Vec, v::Vec, s2w::First2Vector2) = setindex(setindex(v, w[1], 1), w[2], 2)
@unfix workspace2state(w::Vec, v::Vec, s2w::Select2Vector2) = setindex(setindex(v, w[1], s2w.a), w[2], s2w.b)
@unfix workspace2state(w::Vec, v::Vec, s2w::OutputMatrix) = v + Vec(dense(s2w.C)\dense(w - s2w.C*v))
@unfix workspace2state(w::Vec{2}, v::SE2State, s2w::ExtractVector) = SE2State(w, v.t)

### Propagation and Waypoints
setup_steering(SS::StateSpace, r) = setup_steering(SS.dist, r)
setup_steering(d::PreMetric, r) = nothing
setup_steering(d::ChoppedPreMetric, r) = (d.chopval = r)
controltype(d::PreMetric) = NullControl  # TODO: not sure if I even use this anymore, but will be relevant for caching

# function propagate{T<:AbstractFloat}(d::PreMetric,
#                                      v::State, u::Union{ZeroOrderHoldControl, StepControl}, t::AbstractVector{T})
#     path = typeof(v)[]
#     for ui in splitcontrol(u, t)[1:end-1]
#         v = propagate(d, v, ui)
#         push!(path, v)
#     end
#     path
# end
function propagate{T}(d::PreMetric, v::State, u::StepControl{T}, s::AbstractFloat)
    s <= 0 ? v : s >= duration(u) ? propagate(d, v, u) : propagate(d, v, StepControl(T(s),u.u))
end

function propagate(d::PreMetric, v::State, us::ControlSequence)
    for u in us
        v = propagate(d, v, u)
    end
    v
end

function propagate{T<:AbstractFloat}(d::PreMetric, v::State, u::ControlInfo, s::AbstractVector{T})
    [propagate(d, v, u, t) for t in s]
end
function propagate{T<:AbstractFloat}(d::PreMetric, v::State, u::ControlSequence, s::AbstractVector{T})
    issorted(s) || error("Times should be sorted as input to propagate.")
    tf = duration(u)
    s = clamp(s, zero(tf), tf)
    path = typeof(v)[]
    t0, i = zero(tf), 1
    for t in s
        while t >= t0 + duration(u[i]) && i < length(u)    # second clause necessary because of numerical error
            v = propagate(d, v, u[i])
            t0 += duration(u[i])
            i += 1
        end
        push!(path, propagate(d, v, u[i], t-t0))
    end
    path
end

function collision_waypoints(d::PreMetric, v::State, us::ControlSequence)
    path = typeof(v)[]
    for u in us
        append!(path, collision_waypoints(d, v, u))
        v = propagate(d, v, u)
    end
end

## Waypoints
function waypoints(d::PreMetric, v::State, w::State, n::Int)
    u = steering_control(d, v, w)
    t = duration(u)
    propagate(d, v, u, linspace(typeof(t)(0), t, n))
end
function waypoints(d::PreMetric, v::State, w::State, dt::AbstractFloat)
    u = steering_control(d, v, w)
    t = duration(u)
    propagate(d, v, u, [typeof(t)(0):dt:t; t])
end
collision_waypoints(d::PreMetric, v::State, w::State) = push!(collision_waypoints(d, v, steering_control(d,v,w)), w)

## StateSpace and CLBM
for f in (:propagate, :collision_waypoints, :waypoints, :steering_control)
    @eval $f(SS::StateSpace, args...) = $f(SS.dist, args...)
    @eval $f(CLBM::ChoppedPreMetric, args...) = $f(CLBM.m, args...)
end

### Validity Checking
in_state_space{N}(v::Vec{N}, SS::StateSpace) = @all [SS.lo[i] <= v[i] <= SS.hi[i] for i in 1:N]
in_state_space(v::SE2State, SS::SE2StateSpace) = (SS.lo[1] < v.x[1] < SS.hi[1] && SS.lo[2] < v.x[2] < SS.hi[2])
is_free_state(v::State, CC::CollisionChecker, SS::StateSpace) =
    in_state_space(v, SS) && is_free_state(state2workspace(v, SS), CC)
function is_free_motion(v::State, w::State, CC::CollisionChecker, SS::StateSpace)
    wps = collision_waypoints(SS, v, w)
    @all [(in_state_space(wps[i], SS) &&
           is_free_motion(state2workspace(wps[i], SS.s2w), state2workspace(wps[i+1], SS.s2w), CC))
           for i in 1:length(wps)-1]
end
is_free_path(p::Path, CC::CollisionChecker, SS::StateSpace) =
    @all [is_free_motion(p[i], p[i+1], CC, SS) for i in 1:length(p)-1]

### Near Neighbor Setup (should be extended for best performance)
defaultNN(SS::StateSpace, init) = defaultNN(SS.dist, init)
defaultNN(d::Metric, init) = MetricNN(typeof(init)[init], d, init)
defaultNN(d::QuasiMetric, init) = QuasiMetricNN(typeof(init)[init], d, init)
function defaultNN{M}(d::ChoppedPreMetric{M}, init)
    M <: Metric && return MetricNN(typeof(init)[init], d, init)
    M <: QuasiMetric && return QuasiMetricNN(typeof(init)[init], d, init)
    error("Unsupported ChoppedPreMetric")
end
# pairwise_distances(d::PreMetric, V) = # TODO: implement this in a reasonable manner with ControlInfo

### Plotting
plot(SS::StateSpace) = plot_bounds(SS.lo, SS.hi)
plot_waypoints(SS::StateSpace, v::State, w::State) =
    map(x -> state2workspace(x, SS.s2w), isa(SS.dist, Euclidean) ? collision_waypoints(SS,v,w) : waypoints(SS,v,w,10))
function plot_path(p::Path, SS::StateSpace; kwargs...)
    plt.scatter(zip([state2workspace(x, SS)[1:2] for x in p]...)...; kwargs...)
    wps = hcat([hcat(plot_waypoints(SS, p[i], p[i+1])...) for i in 1:length(p)-1]...)
    plot_path(wps; kwargs...)
end
function plot_tree(V::Path, A::Vector{Int}, SS::StateSpace; kwargs...)    # V is a vector of states, i.e. a Path
    VW = [state2workspace(v, SS) for v in V]
    pts = hcat(VW[find(A)]...)
    plt.scatter(pts[1,:], pts[2,:], zorder=1; kwargs...)
    X = vcat([[hcat(plot_waypoints(SS, V[A[w]], V[w])...)[1,:]; nothing] for w in find(A)]...)
    Y = vcat([[hcat(plot_waypoints(SS, V[A[w]], V[w])...)[2,:]; nothing] for w in find(A)]...)
    plt.plot(X, Y, linewidth=.5, linestyle="-", zorder=1; kwargs...)
end

include("statespaces/geometric.jl")
include("statespaces/linearquadratic.jl")
include("statespaces/simplecars.jl")