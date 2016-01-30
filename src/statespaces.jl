export RealVectorStateSpace, SE2StateSpace
export Identity, VectorView, First2Vector2, Select2Vector2, OutputMatrix, ExtractVector
export sample_space, state2workspace, volume, dim, waypoints, workspace_waypoints, collision_waypoints, time_waypoint, cost_waypoint
export in_state_space, is_free_state, is_free_motion, is_free_path, defaultNN, setup_steering, controltype

#====================== Implementation Notes ======================
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
immutable ExtractVector <: State2Workspace end
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
@unfix workspace2state(w::Vec, v::Vec, s2w::OutputMatrix) = v + pinv(s2w.C)*(w - s2w.C*v)
@unfix workspace2state(w::Vec{2}, v::SE2State, s2w::ExtractVector) = SE2State(w, v.t)

### Waypoints
setup_steering(SS::StateSpace, r) = setup_steering(SS.dist, r)
setup_steering(d::PreMetric, r) = nothing
setup_steering(d::ChoppedLowerBoundedPreMetric, r) = (d.chopval = r + eps(r))
controltype(d::PreMetric) = NullControl

waypoints{S<:State}(SS::StateSpace, v::S, w::S) = waypoints(SS.dist, v, w)
time_waypoint{S<:State}(SS::StateSpace, v::S, w::S, t) = time_waypoint(SS.dist, v, w, t)
cost_waypoint{S<:State}(SS::StateSpace, v::S, w::S, c) = cost_waypoint(SS.dist, v, w, c)
workspace_waypoints{S<:State}(SS::StateSpace, v::S, w::S) = workspace_waypoints(SS.dist, v, w, SS.s2w)
collision_waypoints{S<:State}(SS::StateSpace, v::S, w::S) = collision_waypoints(SS.dist, v, w, SS.s2w)
# workspace_waypoints(d::PreMetric, v::State, w::State, ::Identity) = waypoints(d, v, w)

## CLBM
waypoints(CLBM::ChoppedLowerBoundedPreMetric, v, w) = waypoints(CLBM.m, v, w)
time_waypoint(CLBM::ChoppedLowerBoundedPreMetric, v, w, t) = time_waypoint(CLBM.m, v, w, t)
cost_waypoint(CLBM::ChoppedLowerBoundedPreMetric, v, w, c) = cost_waypoint(CLBM.m, v, w, c)
workspace_waypoints(CLBM::ChoppedLowerBoundedPreMetric, v, w) = workspace_waypoints(CLBM.m, v, w, SS.s2w)
collision_waypoints(CLBM::ChoppedLowerBoundedPreMetric, v, w) = collision_waypoints(CLBM.m, v, w, SS.s2w)

## Defaults (should be extended for best performance)
workspace_waypoints(d::PreMetric, v::State, w::State, s2w::State2Workspace) = [state2workspace(x, s2w) for x in waypoints(d, v, w)]
collision_waypoints(d::PreMetric, v::State, w::State, s2w::State2Workspace) = workspace_waypoints(d, v, w, s2w)

### Validity Checking
in_state_space{N}(v::Vec{N}, SS::StateSpace) = @all [SS.lo[i] <= v[i] <= SS.hi[i] for i in 1:N]
in_state_space(v::SE2State, SS::SE2StateSpace) = (SS.lo[1] < v.x[1] < SS.hi[1] && SS.lo[2] < v.x[2] < SS.hi[2])
is_free_state(v::State, CC::CollisionChecker, SS::StateSpace) = in_state_space(v, SS) && is_free_state(state2workspace(v, SS), CC)
function is_free_motion(v::State, w::State, CC::CollisionChecker, SS::StateSpace)
    wps = collision_waypoints(SS, v, w)
    @all [(in_state_space(wps[i], SS) && is_free_motion(wps[i], wps[i+1], CC)) for i in 1:length(wps)-1]
end
is_free_path(p::Path, CC::CollisionChecker, SS::StateSpace) = @all [is_free_motion(p[i], p[i+1], CC, SS) for i in 1:length(p)-1]

### Near Neighbor Setup (should be extended for best performance)
defaultNN(SS::StateSpace, init) = defaultNN(SS.dist, init)
defaultNN(d::Metric, init) = MetricNN(typeof(init)[init], d, init)
defaultNN(d::QuasiMetric, init) = QuasiMetricNN(typeof(init)[init], d, init)
function defaultNN{M}(d::ChoppedLowerBoundedPreMetric{M}, init)
    M <: Metric && return MetricNN(typeof(init)[init], d, init)
    M <: QuasiMetric && return QuasiMetricNN(typeof(init)[init], d, init)
    error("Unsupported ChoppedLowerBoundedPreMetric")
end
# pairwise_distances(d::PreMetric, V) = # TODO: implement this in a reasonable manner with ControlInfo

### Plotting
plot(SS::StateSpace) = plot_bounds(SS.lo, SS.hi)
function plot_path(p::Path, SS::StateSpace; kwargs...)
    wps = hcat([hcat(workspace_waypoints(SS, p[i], p[i+1])...) for i in 1:length(p)-1]...)
    plot_path(wps; kwargs...)
end
function plot_tree(V::Path, A::Vector{Int}, SS::StateSpace; kwargs...)    # V is a vector of states, i.e. a Path
    VW = [state2workspace(v, SS) for v in V]
    pts = hcat(VW[find(A)]...)
    scatter(pts[1,:], pts[2,:], zorder=1; kwargs...)
    X = vcat([[hcat(workspace_waypoints(SS, V[A[w]], V[w])...)[1,:]'; nothing] for w in find(A)]...)
    Y = vcat([[hcat(workspace_waypoints(SS, V[A[w]], V[w])...)[2,:]'; nothing] for w in find(A)]...)
    plt[:plot](X, Y, linewidth=.5, linestyle="-", zorder=1; kwargs...)
end

include("statespaces/geometric.jl")
# include("statespaces/linearquadratic.jl")
include("statespaces/simplecars.jl")