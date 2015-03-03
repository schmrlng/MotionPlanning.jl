module MotionPlanning

# using PyPlot
using ArrayViews
using Devectorize
using Iterators
# using HypothesisTests
using Distances
using NearestNeighbors
using Graphs
using FLANN
using ImmutableArrays

include("problems.jl")
include("nearneighbors.jl")
# include("obstacles.jl")
include("goals.jl")
include("statespaces.jl")
include("sampling.jl")
include("planners.jl")
include("plotting.jl")

end # module
