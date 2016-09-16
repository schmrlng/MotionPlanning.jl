export BoundedStateSpace
export Identity, VectorView, OutputMatrix
export sample_space, state2workspace, volume, dim
export steering_control, propagate, collision_waypoints, waypoints, plot_waypoints
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
immutable BoundedStateSpace{S<:State,M<:PreMetric,W<:State2Workspace} <: StateSpace{S}
    lo::S
    hi::S
    dist::M
    s2w::W
end
changeprecision{T<:AbstractFloat}(::Type{T}, SS::BoundedStateSpace) =
    BoundedStateSpace(map(x -> changeprecision(T,x), (SS.lo, SS.hi, SS.dist, SS.s2w))...)
changeprecision{T<:AbstractFloat}(::Type{T}, s2w::State2Workspace) = s2w

### Sampling
sample_space(SS::BoundedStateSpace) = (SS.lo + rand(typeof(SS.lo)).*(SS.hi-SS.lo))
volume(SS::BoundedStateSpace) = prod(SS.hi-SS.lo)    # TODO: only makes sense for a Euclidean metric
dim{S}(SS::BoundedStateSpace{S}) = length(S)         # TODO: this should really be a function of the metric dist

### State2Workspace (configuration space -> collision checking workspace)
immutable Identity <: State2Workspace end
immutable VectorView{N} <: State2Workspace
    inds::NTuple{N,Int}
end
VectorView(inds::Int...) = VectorView(inds)
VectorView(inds) = VectorView(Tuple(inds))
immutable OutputMatrix{M,N,T<:AbstractFloat,L} <: State2Workspace
    C::SMatrix{M,N,T,L}
end
OutputMatrix(C::Matrix) = OutputMatrix(SMatrix{size(C,1),size(C,2)}(C))
changeprecision{T<:AbstractFloat}(::Type{T}, s2w::OutputMatrix) = OutputMatrix(changeprecision(T, s2w.C))

state2workspace(v::State, SS::StateSpace) = state2workspace(v, SS.s2w)
state2workspace(v::State, s2w::Identity) = v
state2workspace(v::AbstractVector, s2w::VectorView) = v[s2w.inds]
state2workspace(v::AbstractVector, s2w::OutputMatrix) = s2w.C*SVector(v)

workspace2state(w::AbstractVector, SS::StateSpace) = workspace2state(w, SS, SS.s2w)
workspace2state(w::AbstractVector, SS::StateSpace, s2w::Identity) = w
workspace2state(w::AbstractVector, SS::StateSpace, s2w::State2Workspace) = workspace2state(w, sample_space(SS), s2w)
function workspace2state(w::AbstractVector, v::AbstractVector, s2w::VectorView)
    mv = MVector(v)
    mv[s2w.inds] = w
    typeof(v)(mv)
end
workspace2state(w::AbstractVector, v::AbstractVector, s2w::OutputMatrix) = v + typeof(v)(Matrix(s2w.C) \ (w - s2w.C*v))

### Propagation and Waypoints
setup_steering(SS::StateSpace, r) = setup_steering(SS.dist, r)
setup_steering(d::PreMetric, r) = nothing
setup_steering(d::ChoppedPreMetric, r) = (d.chopval = r)
controltype(d::PreMetric) = NullControl  # TODO: not sure if I even use this anymore, but will be relevant for caching

## Generic propagation methods
function propagate{T}(d::Union{Metric, QuasiMetric}, v::State, u::StepControl{T}, s::AbstractFloat)
    s <= 0 ? v : s >= duration(u) ? propagate(d, v, u) : propagate(d, v, StepControl(T(s),u.u))
end
function propagate(d::Union{Metric, QuasiMetric}, v::State, us::ControlSequence)
    for u in us
        v = propagate(d, v, u)
    end
    v
end
function propagate(d::Union{Metric, QuasiMetric}, v::State, us::ControlSequence, s::AbstractFloat)
    s <= 0 && return v
    t = zero(s)
    for u in us
        if s >= t + duration(u)
            v = propagate(d, v, u)
            t += duration(u)
        else
            return propagate(d, v, u, s - t)
        end
    end
    v
end
function propagate{T<:AbstractFloat}(d::Union{Metric, QuasiMetric}, v::State, u::ControlInfo, s::AbstractVector{T})
    [propagate(d, v, u, t) for t in s]
end
function propagate{T<:AbstractFloat}(d::Union{Metric, QuasiMetric}, v::State, u::ControlSequence, s::AbstractVector{T})
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

## Generic waypoints methods
function waypoints(d::Union{Metric, QuasiMetric}, v::State, w::State, n::Int)
    u = steering_control(d, v, w)
    t = duration(u)
    propagate(d, v, u, linspace(typeof(t)(0), t, n))
end
function waypoints(d::Union{Metric, QuasiMetric}, v::State, w::State, dt::AbstractFloat)
    u = steering_control(d, v, w)
    t = duration(u)
    propagate(d, v, u, [typeof(t)(0):dt:t; t])
end
function collision_waypoints(d::Union{Metric, QuasiMetric}, v::State, us::ControlSequence)
    path = typeof(v)[]
    for u in us
        append!(path, collision_waypoints(d, v, u))
        v = propagate(d, v, u)
    end
    path
end
collision_waypoints(d::Union{Metric, QuasiMetric}, v::State, w::State) = push!(collision_waypoints(d, v, steering_control(d,v,w)), w)

## StateSpace and CLBM
for f in (:propagate, :collision_waypoints, :waypoints, :steering_control)
    @eval $f(SS::StateSpace, args...) = $f(SS.dist, args...)
    @eval $f(CLBM::ChoppedPreMetric, args...) = $f(CLBM.m, args...)
end
steering_control(SS::StateSpace, V::State...) = vcat([steering_control(SS.dist,V[i],V[i+1]) for i in 1:length(V)-1]...)

### Validity Checking
in_state_space(v::AbstractVector, SS::StateSpace) = @all [SS.lo[i] <= v[i] <= SS.hi[i] for i in 1:length(v)]
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
    wps = hcat([statevec2mat(plot_waypoints(SS, p[i], p[i+1])) for i in 1:length(p)-1]...)
    plot_path(wps; kwargs...)
end
function plot_tree(V::Path, A::Vector{Int}, SS::StateSpace; kwargs...)    # V is a vector of states, i.e. a Path
    VW = [state2workspace(v, SS) for v in V]
    pts = statevec2mat(VW[find(A)])
    plt.scatter(pts[1,:], pts[2,:], zorder=1; kwargs...)
    X = vcat([[statevec2mat(plot_waypoints(SS, V[A[w]], V[w]))[1,:]; nothing] for w in find(A)]...)
    Y = vcat([[statevec2mat(plot_waypoints(SS, V[A[w]], V[w]))[2,:]; nothing] for w in find(A)]...)
    plt.plot(X, Y, linewidth=.5, linestyle="-", zorder=1; kwargs...)
end

include("statespaces/geometric.jl")
include("statespaces/linearquadratic.jl")
include("statespaces/simplecars.jl")
