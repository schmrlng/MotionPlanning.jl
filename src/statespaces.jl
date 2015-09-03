export RealVectorStateSpace, SE2StateSpace
export Identity, VectorView, First2Vector2, Select2Vector2, ExtractVector
export sample_space, state2workspace, waypoints, workspace_waypoints, collision_waypoints, time_waypoint, cost_waypoint
export inbounds, is_free_state, is_free_motion, is_free_path, defaultNN, setup_steering
import Distances: evaluate

#====================== Implementation Notes ======================
To specify a new state space (with associated distance function M,
a Metric or QuasiMetric) the user must specify:
1. evaluate(M, v, w)
2. waypoints(v, w, M)           # TODO: consider swapping order to match Distances.jl?
3. time_waypoint(v, w, M, t)    (if time-based discretization is desired)
4. cost_waypoint(v, w, N, t)    (if cost-based discretization is desired)

Optionally, to improve performance, the user may specify:
1. defaultNN(M, init)
2. workspace_waypoints(v, w, M, s2w)
3. collision_waypoints(v, w, M, s2w)    (i.e. a sparser version of workspace_waypoints)

In rare circumstances, the user might need to specify:
1. setup_steering

==================================================================#

### State Space Definitions
immutable RealVectorStateSpace{T<:FloatingPoint, M<:PreMetric, W<:State2Workspace} <: StateSpace
    dim::Int
    lo::Vector{T}           # workspace only; 0 <= theta <= 2pi
    hi::Vector{T}
    dist::M
    s2w::W
end

immutable SE2StateSpace{T<:FloatingPoint, M<:PreMetric, W<:State2Workspace} <: StateSpace
    dim::Int
    lo::Vector2{T}           # workspace only; 0 <= theta <= 2pi
    hi::Vector2{T}
    dist::M
    s2w::W
end

### Sampling
sample_space(SS::RealVectorStateSpace) = (SS.lo + rand(SS.dim).*(SS.hi-SS.lo))   # TODO: @devec
sample_space{T,M,W}(SS::SE2StateSpace{T,M,W}) = SE2State(SS.lo + Vector2{T}(rand(T), rand(T)).*(SS.hi-SS.lo), convert(T, 2pi*rand(T)))   # TODO: @devec
function volume(SS::StateSpace)
    !(isa(SS, RealVectorStateSpace) && isa(SS.dist, Euclidean)) && warn("Volume not yet implemented for non-Euclidean metrics!")
    prod(SS.hi-SS.lo)
end

### State2Workspace (configuration space -> collision checking workspace)
type Identity <: State2Workspace end
type VectorView <: State2Workspace
    r::UnitRange{Int}
end
type First2Vector2 <: State2Workspace end
type Select2Vector2 <: State2Workspace
    a::Int
    b::Int
end
type ExtractVector <: State2Workspace end
state2workspace(v, SS::StateSpace) = state2workspace(v, SS.s2w)
state2workspace(v, s2w::Identity) = v
state2workspace(v, s2w::VectorView) = view(v, s2w.r)
state2workspace(v, s2w::First2Vector2) = Vector2(v[1], v[2])    # without staged functions, this is faster than Select2Vector2
state2workspace(v, s2w::Select2Vector2) = Vector2(v[s2w.a], v[s2w.b])
state2workspace(v::SE2State, s2w::ExtractVector) = v.x

workspace2state(w, SS::StateSpace) = workspace2state(w, SS, SS.s2w)
workspace2state(w, SS::StateSpace, s2w::Identity) = w
function workspace2state(w, SS::StateSpace, s2w::State2Workspace)
    v = sample_space(SS)
    workspace2state(w, v, s2w)
end
workspace2state(w::AbstractVector, v::AbstractVector, s2w::VectorView) = (v[s2w.r] = w; v)
workspace2state(w::AbstractVector, v::AbstractVector, s2w::First2Vector2) = (v[1] = w[1]; v[2] = w[2]; v)
workspace2state(w::AbstractVector, v::AbstractVector, s2w::Select2Vector2) = (v[s2w.a] = w[1]; v[s2w.b] = w[2]; v)
workspace2state(w::AbstractVector, v::SE2State, s2w::ExtractVector) = SE2State(Vector2(w), v.t)

### Waypoints
setup_steering(SS::StateSpace, r) = nothing

waypoints{S<:State}(v::S, w::S, SS::StateSpace) = waypoints(v, w, SS.dist)
time_waypoint{S<:State}(v::S, w::S, SS::StateSpace, t) = time_waypoint(v, w, SS.dist, t)
cost_waypoint{S<:State}(v::S, w::S, SS::StateSpace, c) = cost_waypoint(v, w, SS.dist, c)
workspace_waypoints{S<:State}(v::S, w::S, SS::StateSpace) = workspace_waypoints(v, w, SS.dist, SS.s2w)
collision_waypoints{S<:State}(v::S, w::S, SS::StateSpace) = collision_waypoints(v, w, SS.dist, SS.s2w)
# workspace_waypoints(v::State, w::State, d::PreMetric, ::Identity) = waypoints(v, w, d)

## Defaults (should be extended for best performance)
workspace_waypoints(v::State, w::State, d::PreMetric, s2w::State2Workspace) = [state2workspace(x, s2w) for x in waypoints(v, w, d)]
collision_waypoints(v::State, w::State, d::PreMetric, s2w::State2Workspace) = workspace_waypoints(v, w, d, s2w)

### Validity Checking
function inbounds(v::AbstractVector, SS::StateSpace)
    for i in 1:length(v)
        (SS.lo[i] > v[i] || v[i] > SS.hi[i]) && return false
    end
    true
end
inbounds(v::SE2State, SS::SE2StateSpace) = (SS.lo[1] < v.x[1] < SS.hi[1] && SS.lo[2] < v.x[2] < SS.hi[2])
is_free_state(v::State, CC::CollisionChecker, SS::StateSpace) = inbounds(v, SS) && is_free_state(state2workspace(v, SS), CC)
function is_free_motion(v::State, w::State, CC::CollisionChecker, SS::StateSpace)
    wps = collision_waypoints(v, w, SS)
    for i in 1:length(wps)-1
        (!inbounds(wps[i], SS) || !is_free_motion(wps[i], wps[i+1], CC)) && return false
    end
    true
end
function is_free_path(p::Path, CC::CollisionChecker, SS::StateSpace)
    for i in 1:length(p)-1
        !is_free_motion(p[i], p[i+1], CC, SS) && return false
    end
    true
end

### Near Neighbor Setup (should be extended for best performance)
defaultNN(SS::StateSpace, init) = defaultNN(SS.dist, init)
defaultNN(d::Metric, init) = MetricNN_BruteForce(typeof(init)[init], d)
defaultNN(d::QuasiMetric, init) = QuasiMetricNN_BruteForce(typeof(init)[init], d)

### Plotting
plot(SS::StateSpace) = plot_bounds(SS.lo, SS.hi)
function plot_path(p::Path, SS::StateSpace; kwargs...)
    wps = hcat([hcat(workspace_waypoints(p[i], p[i+1], SS)...) for i in 1:length(p)-1]...)
    plot_path(wps; kwargs...)
end
function plot_tree(V::Path, A::Vector{Int}, SS::StateSpace; kwargs...)    # V is a vector of states, i.e. a Path
    VW = [state2workspace(v, SS) for v in V]
    pts = hcat(VW[find(A)]...)
    scatter(pts[1,:], pts[2,:], zorder=1; kwargs...)
    X = vcat([[hcat(workspace_waypoints(V[A[w]], V[w], SS)...)[1,:]', nothing] for w in find(A)]...)
    Y = vcat([[hcat(workspace_waypoints(V[A[w]], V[w], SS)...)[2,:]', nothing] for w in find(A)]...)
    plt.plot(X, Y, linewidth=.5, linestyle="-", zorder=1; kwargs...)
end

include("statespaces/geometric.jl")
# include("statespaces/linearquadratic.jl")
include("statespaces/simplecars.jl")