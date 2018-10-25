# __precompile__()

module MotionPlanning

using LinearAlgebra
using Statistics

using Reexport
@reexport using StaticArrays
@reexport using DifferentialDynamicsModels
using TypeSortedCollections
using OffsetArrays
using DataStructures
using Requires
using RecipesBase

include("simplesets.jl")
include("basictypes.jl")
include("collisioncheckers.jl")
include("goals.jl")
include("samplesets.jl")
include("connections.jl")
include("nearneighbors.jl")
include("graphs.jl")
include("problems.jl")
# include("sampling.jl")
include("planners.jl")
# include("postprocessors.jl")
# include("plotting.jl")

# TODO: test suite!
# TODO: T<:AbstractFloat - either ensure type stability (at least for states), or just go with Float64 everything
# TODO: change up lazy vcat [.,.] syntax in anticipation of v0.4; also consider push!/append!
# TODO: half-baked controls caching - do it or kill it
# TODO: "prototype" (.h-ish) files to specify required methods; julia type system makes OOP a pain
# TODO: FMT re-solving is wonky with different sample sizes
# TODO: Think about parameterizing MPProblem on the StateSpace type (for ease of postprocessing)
# TODO: Requires.jl for conditional dependency loading

# MODERN TODO: broadcast API for ObstacleList
# MODERN TODO: abstract `StateSet`, `ConfigSet` (and maybe `WorkspaceSet`?) for use in `CollisionChecker`s and `Goal`s

function __init__()
    # CollisionChecker glue code
    @require SeparatingAxisTheorem2D="33f35c42-577e-5c66-8569-bbf1771cbf32" include("collisioncheckers/separatingaxistheorem2D.jl")
    # @require GeometryTypes
    # @require Bullet3

    # NearNeighborDataStructure glue code
    @require NearestNeighbors="b8a86587-4115-5ab1-83bc-aa920d37bbce" include("nearneighbors/nearestneighbors.jl")
    @require FLANN="4ef67f76-e0de-5105-ac01-03b6482fb4f8" include("nearneighbors/flann.jl")
end

end # module
