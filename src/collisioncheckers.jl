export Obstacle, StateSpaceObstacle, ConfigSpaceObstacle, WorkspaceObstacle
export RobotCollisionBody, RobotCollisionComponent, RobotWorkspaceComponent, WorkspaceRobotBody, PointRobot
export CollisionChecker, ContinuousCollisionChecker, DiscreteCollisionChecker
export is_free_state, is_free_motion, is_free_config, is_free_sweep, rand_free_state, intersecting, sweep_intersecting

# Obstacle Sets
struct Obstacle{T,S}
    set::S
end
const StateSpaceObstacle{S}  = Obstacle{:statespace,S}
const ConfigSpaceObstacle{S} = Obstacle{:configspace,S}
const WorkspaceObstacle{S}   = Obstacle{:workspace,S}
StateSpaceObstacle(set::S) where {S}  = Obstacle{:statespace,S}(set)
ConfigSpaceObstacle(set::S) where {S} = Obstacle{:configspace,S}(set)
WorkspaceObstacle(set::S) where {S}   = Obstacle{:workspace,S}(set)
intersecting(o::Obstacle, x) = intersecting(o.set, x)
sweep_intersecting(o::Obstacle, x0, xf) = sweep_intersecting(o.set, x0, xf)
@recipe function f(o::Obstacle; dims=(1, 2))#, obstacle_color=:red, obstacle_alpha=1)
    dims  --> dims
    # color :=  obstacle_color
    # alpha :=  obstacle_alpha
    # label --> ""
    o.set
end

# Abstract Robot Types (speculative for now)
abstract type RobotCollisionBody end
abstract type RobotCollisionComponent end
abstract type RobotWorkspaceComponent <: RobotCollisionComponent end

# Collision Checker
struct CollisionChecker{C,O,R<:RobotCollisionBody,S2C}
    obstacles::O
    robot::R
    state2config::S2C
    motion_count::Base.RefValue{Int}
    edge_count::Base.RefValue{Int}
end
const ContinuousCollisionChecker{O,R,S2C} = CollisionChecker{true,O,R,S2C}
const DiscreteCollisionChecker{O,R,S2C}   = CollisionChecker{false,O,R,S2C}
reset!(CC::CollisionChecker) = (CC.motion_count[] = CC.edge_count[] = 0)
@recipe function f(CC::CollisionChecker; dims=(1, 2))#, obstacle_color=:red, obstacle_alpha=1)
    dims           --> dims
    # obstacle_color --> obstacle_color
    # obstacle_alpha --> obstacle_alpha
    foreach(o -> @series(begin o end), CC.obstacles)
end

# General Methods for Dispatch
## is_free_state
is_free_state(CC::CollisionChecker, x) = all(o -> is_free_state(o, CC.robot, CC.state2config, x), CC.obstacles)
is_free_state(o::Obstacle, robot, s2c, x) = is_free_config(o, robot, s2c(x))
is_free_state(o::StateSpaceObstacle,  robot, s2c, x) = !intersecting(o, x)

## is_free_motion
function is_free_motion(CC::CollisionChecker, x0, xf)
    CC.motion_count[] += 1
    all(o -> is_free_motion(o, CC.robot, CC.state2config, x0, xf), CC.obstacles)
end
is_free_motion(o::Obstacle, robot, s2c, x0, xf) = is_free_sweep(o, robot, s2c(x0), s2c(xf))
is_free_motion(o::StateSpaceObstacle,  robot, s2c, x0, xf) = !sweep_intersecting(o, x0, xf)

## is_free_config
is_free_config(o::ConfigSpaceObstacle, robot, q) = !intersecting(o, q)

## is_free_sweep
is_free_sweep(o::ConfigSpaceObstacle, robot, q0, qf) = !sweep_intersecting(o, q0, qf)

## rand_free_state
rand_free_state(CC::CollisionChecker, X) = (x = rand(X); is_free_state(CC, x) ? x : rand_free_state(CC, X))

# WorkspaceRobotBody
struct WorkspaceRobotBody{C2W} <: RobotCollisionBody
    config2workspace::C2W
end
is_free_config(o::WorkspaceObstacle, robot::WorkspaceRobotBody, q) = !intersecting(o, robot.config2workspace(q))
function is_free_sweep(o::WorkspaceObstacle, robot::WorkspaceRobotBody, q0, qf)
    !sweep_intersecting(o, robot.config2workspace(q0), robot.config2workspace(qf))
end

PointRobot() = WorkspaceRobotBody(identity)

# SingleLinkRobot
struct SingleLinkRobot{L,C2T} <: RobotCollisionBody
    link::L
    config2transform::C2T
end
function is_free_config(o::WorkspaceObstacle, robot::SingleLinkRobot, q)
    !intersecting(o, robot.link, robot.config2transform(q))
end
function is_free_sweep(o::WorkspaceObstacle, robot::SingleLinkRobot, q0, qf)
    !sweep_intersecting(o, robot.link, robot.config2transform(q0), robot.config2transform(qf))
end

## Convenience Constructors
_type_sort(obstacles) = obstacles
_type_sort(obstacles::Vector{T}) where {T} = isconcretetype(T) ? obstacles : TypeSortedCollection(obstacles)
function ContinuousCollisionChecker(_obstacles, robot::R=PointRobot(), state2config::S2C=identity) where {R,S2C}
    obstacles = _type_sort(_obstacles)
    CollisionChecker{true,typeof(obstacles),R,S2C}(obstacles, robot, state2config, Ref(0), Ref(0))
end
function DiscreteCollisionChecker(_obstacles, robot::R=PointRobot(), state2config::S2C=identity) where {R,S2C}
    obstacles = _type_sort(_obstacles)
    CollisionChecker{false,typeof(obstacles),R,S2C}(obstacles, robot, state2config, Ref(0), Ref(0))
end

# Super basic collision checking; `import SeparatingAxisTheorem2D` or ... to load more advanced capabilities
intersecting(X::Union{SimpleConvexSet,SetComplement}, x) = x in X

function sweep_intersecting(B::AxisAlignedBox{N}, x0, xf) where {N}
    all(B.lo[i] <= max(x0[i], xf[i]) && min(x0[i], xf[i]) <= B.hi[i] for i in 1:N) && _sweep_intersecting(B, x0, xf)
end
function _sweep_intersecting(B::AxisAlignedBox{N}, x0, xf) where {N}
    any((λ = (ifelse(x0[i] < B.lo[i], B.lo[i], B.hi[i]) - x0[i])/(xf[i] - x0[i]);
        0 <= λ <= 1 && all(i == j || B.lo[j] <= x0[j] + (xf[j] - x0[j])*λ <= B.hi[j] for j in 1:N)) for i in 1:N)
end

function sweep_intersecting(B::Ball, x0, xf)
    all(B.c[i] - B.r <= max(x0[i], xf[i]) &&
        min(x0[i], xf[i]) <= B.c[i] + B.r for i in 1:N) && _sweep_intersecting(B, x0, xf)
end
function _sweep_intersecting(B::Ball, x0, xf)
    x0 in B && return true
    xf in B && return true
    X = xf - x0
    Y = B.c - x0
    a = dot(X, X)
    b = dot(Y, Y)
    c = dot(X, Y)
    a*C.r*C.r < a*b - c*c && return false
    0 <= c <= a
end

sweep_intersecting(S::SetComplement{<:SimpleConvexSet}, x0, xf) = x0 in S || xf in S    # S.set is convex
