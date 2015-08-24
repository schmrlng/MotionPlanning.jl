export StateSpace, GeometricStateSpace, RealVectorStateSpace, DifferentialStateSpace
export vector_to_state, sample_space, volume, defaultNN, pairwise_distances, setup_steering

### State Space Definitions
abstract StateSpace

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
sample_space{T}(SS::ReedsSheppStateSpace{T}) = RSState(SS.lo + Vector2{T}(rand(T), rand(T)).*(SS.hi-SS.lo), convert(T, 2pi*rand(T)))   # TODO: @devec

### State2Workspace (configuration space -> collision checking workspace)
abstract State2Workspace
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

### Validity Checking
inbounds(v::AbstractVector, SS::StateSpace) = (SS.lo[1] < v[1] < SS.hi[1] && SS.lo[2] < v[2] < SS.hi[2])
inbounds(v::SE2State, SS::SE2StateSpace) = (SS.lo[1] < v.x[1] < SS.hi[1] && SS.lo[2] < v.x[2] < SS.hi[2])
is_free_state(v::State, CC::CollisionChecker, SS::StateSpace) = inbounds(v, SS) && is_free_state(state2workspace(v, SS.s2w), CC)
function is_free_motion(v::State, w::State, CC::CollisionChecker, SS::StateSpace)
    wps = workspace_waypoints(v, w, SS)
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

### Plotting
plot(SS::StateSpace) = plot_bounds(SS.lo, SS.hi)
function plot_path(p::Path, SS::StateSpace)
    wps = hcat([hcat(workspace_waypoints(p[i], p[i+1], SS.dist, SS.s2w)...) for i in 1:length(p)-1]...)
end


### NearNeighborCache Stuff (TODO: could be a macro, perhaps)
plot_path(SS::StateSpace)

### DEFAULTS - essentially this is all just for GeometricStateSpace; having general defaults may dangerous

is_free_state(v, CC::CollisionChecker, SS::StateSpace) = is_free_state(v, CC)
is_free_motion(v, w, CC::CollisionChecker, SS::StateSpace) = is_free_motion(v, w, CC)
is_free_path(path, CC::CollisionChecker, SS::StateSpace) = is_free_path(path, CC)
plot(SS::StateSpace) = plot_bounds(SS.lo, SS.hi)
plot_path(SS::StateSpace, NN::NearNeighborCache, args...; kwargs...) = plot_path(NN.V, args...; kwargs...)
plot_tree(SS::StateSpace, NN::NearNeighborCache, args...; kwargs...) = plot_tree(NN.V, args...; kwargs...)
setup_steering(SS::StateSpace, r) = nothing

# include("statespaces/geometric.jl")
# include("statespaces/linearquadratic.jl")
include("statespaces/simplecars.jl")
# include("statespaces/reedsshepp.jl")