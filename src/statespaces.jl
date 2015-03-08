export StateSpace, GeometricStateSpace, RealVectorStateSpace, DifferentialStateSpace
export vector_to_state, sample_space

abstract StateSpace

abstract GeometricStateSpace <: StateSpace  # straight-line connections
abstract RealVectorStateSpace <: GeometricStateSpace  # state space is a hyperrectangle

abstract DifferentialStateSpace <: StateSpace  # connections subject to differential constraints

### DEFAULTS

## For state spaces with states that are Vectors

is_free_state(v, CC::CollisionChecker, SS::GeometricStateSpace) = is_free_state(v, CC)
is_free_motion(v, w, CC::CollisionChecker, SS::GeometricStateSpace) = is_free_motion(v, w, CC)
is_free_path(path, CC::CollisionChecker, SS::GeometricStateSpace) = is_free_path(path, CC)
plot(SS::StateSpace) = plot_bounds(SS.lo, SS.hi)

include("statespaces/geometric.jl")
# include("statespaces/linearquadratic.jl")