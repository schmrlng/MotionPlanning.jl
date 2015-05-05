export StateSpace, GeometricStateSpace, RealVectorStateSpace, DifferentialStateSpace
export vector_to_state, sample_space, volume, defaultNN, pairwise_distances, setup_steering

abstract StateSpace

abstract GeometricStateSpace <: StateSpace  # straight-line connections
abstract RealVectorStateSpace <: GeometricStateSpace  # state space is a hyperrectangle

abstract DifferentialStateSpace <: StateSpace  # connections subject to differential constraints

### DEFAULTS - essentially this is all just for GeometricStateSpace; having general defaults may dangerous

is_free_state(v, CC::CollisionChecker, SS::StateSpace) = is_free_state(v, CC)
is_free_motion(v, w, CC::CollisionChecker, SS::StateSpace) = is_free_motion(v, w, CC)
is_free_path(path, CC::CollisionChecker, SS::StateSpace) = is_free_path(path, CC)
plot(SS::StateSpace) = plot_bounds(SS.lo, SS.hi)
plot_path(SS::StateSpace, NN::NearNeighborCache, args...; kwargs...) = plot_path(NN.V, args...; kwargs...)
plot_tree(SS::StateSpace, NN::NearNeighborCache, args...; kwargs...) = plot_tree(NN.V, args...; kwargs...)
setup_steering(SS::StateSpace, r) = nothing

include("statespaces/geometric.jl")
include("statespaces/linearquadratic.jl")
include("statespaces/reedsshepp.jl")