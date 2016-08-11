# __precompile__()

module MotionPlanning

using Reexport
@reexport using StaticArrays
@reexport using Distances
import PyPlot; const plt = PyPlot
using Distributions; import Base.eltype; eltype{T}(::Dirichlet{T}) = T    # TODO: submit Distributions PR
using SCS
using Devectorize
using Iterators
using NearestNeighbors

### Includes
include("primitivetypes.jl")
include("utilities/utils.jl")
include("collisioncheckers.jl")
include("goals.jl")
include("nearneighbors.jl")
include("statespaces.jl")
include("problems.jl")
include("sampling.jl")
include("planners.jl")
include("postprocessors.jl")
include("plotting.jl")

# TODO: test suite!
# TODO: T<:AbstractFloat - either ensure type stability (at least for states), or just go with Float64 everything
# TODO: change up lazy vcat [.,.] syntax in anticipation of v0.4; also consider push!/append!
# TODO: half-baked controls caching - do it or kill it
# TODO: "prototype" (.h-ish) files to specify required methods; julia type system makes OOP a pain
# TODO: FMT re-solving is wonky with different sample sizes 
# TODO: Think about parameterizing MPProblem on the StateSpace type (for ease of postprocessing)
# TODO: Requires.jl for conditional dependency loading

end # module
