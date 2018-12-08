export rrt!

function rrt!(P::MPProblem; max_sample_count=Inf,
                            time_limit=Inf,
                            steering_eps=0.1,    # TODO: compute from P.state_space dims?
                            goal_bias=0.05,
                            compute_full_metadata=true)
    metadata = @standard_setup!(P)
    @assert(isfinite(max_sample_count) || isfinite(time_limit),
            "At least one of max_sample_count and time_limit must be finite.")
    metadata[:planner] = :RRT
    metadata[:max_sample_count] = max_sample_count
    metadata[:time_limit] = time_limit
    metadata[:steering_eps] = steering_eps
    metadata[:goal_bias] = goal_bias

    # Graph Setup
    S = typeof(P.init)
    nodes = ExplicitSampleSet(P.init, S[], [P.init])    # nodes.V = [P.init] to work around FLANN empty index segfault
    P.graph = NearNeighborGraph(nodes, P.bvp, near_style=Val(:variable))

    # Solve
    rrt!(P.state_space, P.bvp, P.init, P.goal, P.collision_checker, P.graph, metadata, P.solution.elapsed)

    # Post-Processing
    if compute_full_metadata
        record_tree!(metadata, metadata[:node_info])
    end

    standard_wrapup!(P)
    P.solution
end

function rrt!(state_space::StateSpace,
              bvp::SteeringBVP,
              init::State,
              goal::Goal,
              collision_checker::CollisionChecker,
              graph::NearNeighborGraph{NeighborInfo{X,D,U}},
              metadata::Dict{Symbol,Any},
              planner_start_time=time(),
              max_sample_count=metadata[:max_sample_count],
              time_limit=metadata[:time_limit],
              steering_eps=D(metadata[:steering_eps]),
              goal_bias=metadata[:goal_bias],
              node_info=node_info_datastructure(graph.nodes, TreeNodeInfo{X,D})) where {X,D,U}
    i = 1
    while i < max_sample_count && time() - planner_start_time < time_limit
        x_rand = rand_free_state(collision_checker, rand() < goal_bias ? goal : state_space)
        v_near, cost, controls = first(neighbors(graph, x_rand, Nearest(1), dir=Val(:B)))
        x_near = graph[v_near]
        x_new, cost, controls = steer_towards(bvp, x_near, x_rand, steering_eps, cost, controls)
        if is_free_edge(collision_checker, bvp, x_near, x_new, controls)
            i += 1
            addstates!(graph, x_new)
            push!(node_info, (parent=v_near, cost_to_come=node_info[v_near].cost_to_come + cost))
        end
        x_new in goal && break
    end

    metadata[:solved] = graph[i] in goal
    record_solution!(metadata, node_info, i, 0)
    nothing
end
