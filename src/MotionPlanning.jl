module MotionPlanning

using PyPlot
using ArrayViews
using Devectorize
using Iterators
using GZip
using Grid
# using HypothesisTests
using Distances
using KDTrees
# using Graphs
using ImmutableArrays

### Planning Primitives
export AbstractState, State, Path

abstract AbstractState
typealias State Union(AbstractVector, AbstractState)
typealias Path{T<:State} Vector{T}

### Includes
include("utilities/collections.jl")
include("collisioncheckers.jl")
include("nearneighbors.jl")
include("goals.jl")
include("statespaces.jl")
include("problems.jl")
include("sampling.jl")
include("planners.jl")
include("postprocessors.jl")
include("plotting.jl")

# TODO: test suite!
# TODO: T<:FloatingPoint - either ensure type stability (at least for states), or just go with Float64 everything
# TODO: change up lazy vcat [.,.] syntax in anticipation of v0.4; also consider push!/append!
# TODO: half-baked controls caching - do it or kill it
# TODO: "prototype" (.h-ish) files to specify required methods; julia type system makes OOP a pain
# TODO: FMT re-solving is wonky with different sample sizes 
# TODO: Think about parameterizing MPProblem on the StateSpace type (for ease of postprocessing)

end # module
