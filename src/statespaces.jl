export StateSpace, GeometricStateSpace, RealVectorStateSpace, DifferentialStateSpace
export vector_to_state, sample_space

abstract StateSpace

abstract GeometricStateSpace <: StateSpace  # straight-line connections
abstract RealVectorStateSpace <: GeometricStateSpace  # state space is a hyperrectangle

abstract DifferentialStateSpace <: StateSpace  # connections subject to differential constraints

### DEFAULTS

## For state spaces with states that are Vectors

vector_to_state(v::Vector, SS::StateSpace) = v  # for goal sampling
sample_space(SS::StateSpace) = (SS.lo + rand(SS.dim).*(SS.hi-SS.lo))
is_free_state(v, CC::CollisionChecker, SS::GeometricStateSpace) = is_free_state(v, CC)
is_free_motion(v, w, CC::CollisionChecker, SS::GeometricStateSpace) = is_free_motion(v, w, CC)
is_free_path(path, CC::CollisionChecker, SS::GeometricStateSpace) = is_free_path(path, CC)

## GeometricStateSpace path checking

# is_free_motion(v::Vector, w::Vector, obs::ObstacleSet, SS::GeometricStateSpace) = is_free_motion(v, w, obs)
# is_free_path(path::Matrix, obs::ObstacleSet, SS::GeometricStateSpace) = is_free_path(path, obs)

# include("RealVectorMetricSpaces.jl")
# include("BoundedEuclideanStateSpace.jl")
# include("BoundedEuclideanStateSpaceCost.jl")
# include("BESSCuboidUnionRobot.jl")
# include("ReedsSheppStateSpace.jl")
# include("LQDMPStateSpace.jl")
# include("ReedsSheppStateSpace.jl")
include("statespaces/geometric.jl")
# include("statespaces/linearquadratic.jl")