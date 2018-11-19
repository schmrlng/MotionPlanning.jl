export rrtstar!

function update_descendant_costs_to_come!(node_info, v, Δc)
    v_info = node_info[v]
    node_info[v] = (v_info..., cost_to_come=v_info.cost_to_come + Δc)
    foreach(w -> update_descendant_costs_to_come!(node_info, w, Δc), v_info.children)
end

function rrtstar!(P::MPProblem; max_sample_count=Inf,
                                time_limit=Inf,
                                steering_eps=0.1,    # TODO: compute from P.state_space dims?
                                radius_scale_factor=1,
                                d=length(P.init),
                                gamma=(2*(1 + 1/d)*volume(P.state_space)/volume(Ball{d}))^(1/d),
                                r=(i -> radius_scale_factor*gamma*(log(i)/i)^(1/d)),
                                goal_bias=0.05,
                                compute_full_metadata=true)
    metadata = @standard_setup!(P)
    @assert(isfinite(max_sample_count) || isfinite(time_limit),
            "At least one of max_sample_count and time_limit must be finite.")
    metadata[:planner] = :RRTstar
    metadata[:max_sample_count] = max_sample_count
    metadata[:time_limit] = time_limit
    metadata[:steering_eps] = steering_eps
    metadata[:radius_scale_factor] = radius_scale_factor
    metadata[:d] = d
    metadata[:gamma] = gamma
    metadata[:r] = r
    metadata[:goal_bias] = goal_bias

    # Graph Setup
    S = typeof(P.init)
    nodes = ExplicitSampleSet(P.init, S[], [P.init])    # nodes.V = [P.init] to work around FLANN empty index segfault
    P.graph = NearNeighborGraph(nodes, P.bvp, near_style=Val(:variable))

    # Solve
    rrtstar!(P.state_space, P.bvp, P.init, P.goal, P.collision_checker, P.graph, metadata, P.solution.elapsed)

    # Post-Processing
    if compute_full_metadata
        record_tree!(metadata, metadata[:node_info])
    end

    standard_wrapup!(P)
    P.solution
end

function rrtstar!(state_space::StateSpace,
                  bvp::SteeringBVP,
                  init::State,
                  goal::Goal,
                  collision_checker::CollisionChecker,
                  graph::NearNeighborGraph{NeighborInfo{X,D,U}},
                  metadata::Dict{Symbol,Any},
                  planner_start_time=time(),
                  r::R=metadata[:r],
                  max_sample_count=metadata[:max_sample_count],
                  time_limit=metadata[:time_limit],
                  steering_eps=D(metadata[:steering_eps]),
                  goal_bias=metadata[:goal_bias],
                  node_info=node_info_datastructure(graph.nodes, FullTreeNodeInfo{X,D})) where {X,D,U,R}
    i = 1
    while i < max_sample_count && time() - planner_start_time < time_limit
        x_rand = rand_free_state(collision_checker, rand() < goal_bias ? goal : state_space)
        v_near = first(neighbors(graph, x_rand, Nearest(1), dir=Val(:B))).index    # TODO: include_controls=Val(true)
        x_near = graph[v_near]
        x_new  = steer_towards(bvp, x_near, x_rand, steering_eps)
        if is_free_edge(collision_checker, bvp, x_near, x_new)    # more TODO
            ri = min(D(r(i)), steering_eps)
            i += 1
            addstates!(graph, x_new)
            Δc_min = bvp(x_near, x_new).cost
            v_min, c_min = v_near, node_info[v_near].cost_to_come + Δc_min
            for (v, Δc, u) in neighbors(graph, i, ri, dir=Val(:B))
                c = node_info[v].cost_to_come + Δc
                if c < c_min && is_free_edge(collision_checker, bvp, graph[v], x_new, u)
                    v_min, c_min, Δc_min = v, c, Δc
                end
            end

            i_info = FullTreeNodeInfo{X,D}(parent=v_min, cost_to_come=c_min)
            push!(node_info, i_info)
            push!(node_info[v_min].children, i)

            for (v, Δc, u) in neighbors(graph, i, ri, dir=Val(:F))
                v_info = node_info[v]
                if c_min + Δc < v_info.cost_to_come && is_free_edge(collision_checker, bvp, x_new, graph[v], u)
                    p_info = node_info[v_info.parent]
                    deleteat!(p_info.children, searchsortedfirst(p_info.children, v))
                    push!(i_info.children, v)
                    node_info[v] = (v_info..., parent=i)
                    update_descendant_costs_to_come!(node_info, v, c_min + Δc - v_info.cost_to_come)
                end
            end
            sort!(i_info.children)
        end
    end

    v_opt = 0
    c_opt = D(Inf)
    solved = false
    for v in eachindex(graph.nodes)
        if graph[v] in goal
            solved = true
            c = node_info[v].cost_to_come
            if c < c_opt
                v_opt = v
                c_opt = c
            end
        end
    end
    metadata[:solved] = solved
    record_tree_solution!(metadata, node_info, v_opt, 0)
    nothing
end
