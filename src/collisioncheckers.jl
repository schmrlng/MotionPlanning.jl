export CollisionChecker, DiscreteCollisionChecker, SweptCollisionChecker
export is_free_state, is_free_motion, is_free_path, inflate, close, closest

abstract CollisionChecker
abstract DiscreteCollisionChecker <: CollisionChecker
abstract SweptCollisionChecker <: CollisionChecker

include("collisioncheckers/SAT2D.jl")
include("collisioncheckers/robots2D.jl")
include("collisioncheckers/bvls.jl")
include("collisioncheckers/boxesND.jl")