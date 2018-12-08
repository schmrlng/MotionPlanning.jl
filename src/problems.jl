export MPProblem, MPSolution

struct MPSolution{T}
    status::Symbol
    cost::T
    elapsed::Float64
    metadata::Dict{Symbol,Any}
end

mutable struct MPProblem
    state_space::StateSpace
    bvp::SteeringBVP
    init::State
    goal::Goal
    collision_checker::CollisionChecker
    status::Symbol
    graph::NearNeighborGraph
    solution::MPSolution

    function MPProblem(state_space::StateSpace,
                       bvp::SteeringBVP,
                       init::State,
                       goal::Goal,
                       collision_checker::CollisionChecker)
        new(state_space, bvp, init, goal, collision_checker, :unsolved)
    end
end

@recipe function f(P::MPProblem; dims=(1, 2),
                   statespace_color=:black, statespace_alpha=1,
                   goal_color=:green, goal_alpha=1, goal_markershape=:star, goal_markersize=10,
                   obstacle_color=:red, obstacle_alpha=1,
                   show_graph=false, graph_color=:grey, graph_alpha=0.5, graph_markersize=2, graph_edge_waypoints=10,
                   solution_color=:blue, solution_alpha=1, solution_markersize=2, solution_edge_waypoints=10)
    dims  --> dims
    label --> ""
    x, y = dims
    @series begin
        color --> statespace_color
        alpha --> statespace_alpha
        P.state_space
    end
    @series begin
        color       --> goal_color
        alpha       --> goal_alpha
        markershape --> goal_markershape
        markersize  --> goal_markersize
        P.goal
    end
    @series begin
        color --> obstacle_color
        alpha --> obstacle_alpha
        P.collision_checker
    end
    if show_graph && isdefined(P, :solution) && :graph in keys(P.solution.metadata)
        @series begin
            color          --> graph_color
            alpha          --> graph_alpha
            markersize     --> graph_markersize
            num_waypoints  --> graph_edge_waypoints
            state2config   --> P.collision_checker.state2config
            # config2viz     --> P.collision_checker.config2viz    # TODO
            plot_endpoints --> true
            plot_x0        --> false
            plot_xf        --> true
            graph = P.solution.metadata[:graph]
            graph isa ForwardAdjacencyDict ? [P.graph[i, j] for (i, js) in P.solution.metadata[:graph].d for j in js] :
                                             [P.graph[i, j] for (j, is) in P.solution.metadata[:graph].d for i in is]
        end
    end
    if P.status === :solved
        @series begin
            color          --> solution_color
            alpha          --> solution_alpha
            markersize     --> solution_markersize
            num_waypoints  --> solution_edge_waypoints
            state2config   --> P.collision_checker.state2config
            # config2viz     --> P.collision_checker.config2viz    # TODO
            plot_endpoints --> true
            [P.graph[i, j] for (i, j) in pairwise(P.solution.metadata[:solution_nodes])]
        end
    end
end
