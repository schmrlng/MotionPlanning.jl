module MotionPlanning

using PyPlot
using ArrayViews
using Devectorize
using Iterators
using HypothesisTests
using Distances
using NearestNeighbors

export
    # Obstacles
    ObstacleSet,
    ObstacleList,
    AABoxes,
    volume,
    is_free_pt,
    is_free_motion,
    is_free_path,
    closest_obs_pt,
    plot_obstacles,

    # Goal Regions
    RectangleGoal,
    BallGoal,
    PointGoal,
    plot_goal,
    is_goal_pt,
    sample_goal,

    # State Spaces
    StateSpace,
    GeometricStateSpace,
    RealVectorStateSpace,
    DifferentialStateSpace,
    NearNeighborCache,
    ProblemSetup,
    vector_to_state,
    is_valid_state,
    is_free_state,
    sample_space,
    steer,

        ## Geometric State Spaces
        RealVectorMetricSpace,
        GeometricProblem,
        BoundedEuclideanStateSpace,
        UnitHypercube,
        shortcut,
        adaptive_shortcut,
        MetricNN,
        MetricNNKDTree,
        Neighborhood,
        filter_neighborhood,
        nearestk,
        nearestr,
        mutualnearestk,

        ## Quasi Metric State Spaces
        LinearQuadraticStateSpace,
        QuasiMetricProblem,
        QuasiMetricNN,
        pairwise_distances,
        pairwise_distances_approx_opt,
        NNCache,
        nearRF,
        nearRB,
        waypoints,
        plot_solution,

    # Sampling
    sample_free_goal,
    sample_free,

    # Planners
    fmtstar,

    # Plotting
    plot_rectangle,
    plot_circle,
    plot_ellipse,
    plot_bounds,
    plot_graph,
    plot_tree,
    plot_path,
    plot_problem_setup

include("obstacles.jl")
include("goals.jl")
include("statespaces.jl")
include("sampling.jl")
include("planners.jl")
include("plotting.jl")

end # module
