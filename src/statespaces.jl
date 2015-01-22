abstract StateSpace

abstract GeometricStateSpace <: StateSpace  # straight-line connections
abstract RealVectorStateSpace <: GeometricStateSpace  # state space is a hyperrectangle

abstract DifferentialStateSpace <: StateSpace  # connections subject to differential constraints

abstract NearNeighborCache

abstract ProblemSetup

### DEFAULTS

## For state spaces with states that are Vectors

vector_to_state(v::Vector, SS::StateSpace) = v
is_valid_state(v::Vector, SS::StateSpace) = true # all(SS.lo .<= v .<= SS.hi) - with the current sampling scheme, this will never be false?
is_free_state(v::Vector, obs::ObstacleSet, SS::StateSpace) = (is_valid_state(v, SS) && is_free_pt(v, obs))
sample_space(SS::StateSpace) = (SS.lo + rand(SS.dim).*(SS.hi-SS.lo))

## GeometricStateSpace path checking

is_free_motion(v::Vector, w::Vector, obs::ObstacleSet, SS::GeometricStateSpace) = is_free_motion(v, w, obs)
is_free_path(path::Matrix, obs::ObstacleSet, SS::GeometricStateSpace) = is_free_path(path, obs)

# include("RealVectorMetricSpaces.jl")
# include("BoundedEuclideanStateSpace.jl")
# include("BoundedEuclideanStateSpaceCost.jl")
# include("BESSCuboidUnionRobot.jl")
# include("ReedsSheppStateSpace.jl")
# include("LQDMPStateSpace.jl")
# include("ReedsSheppStateSpace.jl")
include("statespaces/geometric.jl")
include("statespaces/linearquadratic.jl")