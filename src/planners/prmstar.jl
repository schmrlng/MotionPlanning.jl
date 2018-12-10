export prmstar!

# status = 'U', 'O', or 'C' (unvisited, open, closed)
const PRMNodeInfo{X,D} = NamedTuple{(:is_free, :status, :forward_neighbors, :parent, :cost_to_come),
                                    Tuple{Bool,Char,Vector{X},X,D}}
function PRMNodeInfo{X,D}(; is_free=false, status='U',
                            forward_neighbors=X[], parent=zero(X), cost_to_come=D(Inf)) where {X,D}
    (is_free=is_free, status=status, forward_neighbors=forward_neighbors, parent=parent, cost_to_come=cost_to_come)
end

function prmstar_default_radius_function(state_space::StateSpace, bvp::GeometricSteering, radius_scale_factor=1)
    d = dimension(state_space)
    γ = radius_scale_factor*2*((1/d)*volume(state_space)/volume(Ball{d}))^(1/d)
    N -> γ*(log(N)/N)^(1/d)
end
function prmstar_default_radius_function(state_space::StateSpace, bvp::SteeringBVP)
    error("prmstar_default_radius_function not implemented yet for SteeringBVP type $(typeof(bvp))")
end

function prmstar!(P::MPProblem; use_current_graph=false,
                                sample_count=0,
                                only_sample_free=true,
                                radius_scale_factor=1,
                                r=rrtstar_default_radius_function(P.state_space, P.bvp, radius_scale_factor),
                                ensure_goal_count=0,
                                compute_full_metadata=true)
    metadata = @standard_setup!(P)
    metadata[:planner] = :PRMstar
    metadata[:r] = r

    # Graph Setup
    if !use_current_graph || !isdefined(P, :graph)
        S = typeof(P.init)
        samples = only_sample_free ? S[rand(P.state_space) for i in 1:sample_count] :
                                     S[rand_free_state(P.collision_checker, P.state_space) for i in 1:sample_count]
        nodes = ExplicitSampleSet(P.init, S[], samples)
        P.graph = NearNeighborGraph(nodes, P.bvp)
    end
    setradius!(P.graph, r isa Function ? r(length(P.graph)) : r)    # TODO: some brittleness with `TiledSampleSet`s
    setinit!(P.graph, P.init)
    setgoal!(P.graph, () -> rand_free_state(P.collision_checker, P.goal), sample_count=ensure_goal_count)
    compute_init_and_goal_neighbors!(P.graph)

    # Solve
    prmstar!(P.state_space, P.bvp, P.init, P.goal, P.collision_checker, P.graph, metadata)

    # Post-Processing
    if compute_full_metadata
        record_graph!(metadata, metadata[:node_info], omit=(((k, info),) -> info.status == 'U'))
    end

    standard_wrapup!(P)
    P.solution
end

function prmstar!(state_space::StateSpace,
                  bvp::SteeringBVP,
                  init::State,
                  goal::Goal,
                  collision_checker::CollisionChecker,
                  graph::NearNeighborGraph{NeighborInfo{X,D,U}},
                  metadata::Dict{Symbol,Any},
                  node_info=node_info_datastructure(
                      graph.nodes, PRMNodeInfo{X,D},
                      x -> PRMNodeInfo{X,D}(is_free=is_free_state(collision_checker, x)),
                      bounds=boundingbox(state_space))
                  ) where {X,D,U}
    is_free = n -> (ni = node_info[n.index]; ni !== missing && ni.is_free)

    open_queue = PriorityQueue{X,D}()
    z = X(0)
    node_info[z] = (node_info[z]..., status='O', cost_to_come=zero(D))

    while !(graph[z] in goal)    # uniform-cost search (to accommodate infinite `TiledSampleSet`s)
        z_info = node_info[z]
        for (x, c, u) in neighbors(is_free, graph, z, dir=Val(:F))    # includes status == 'C', in the spirit of PRM
            if is_free_edge(collision_checker, bvp, graph[z], graph[x], u)
                push!(z_info.forward_neighbors, x)
                x_info = node_info[x]
                x_info.status == 'C' && continue
                cost_through_z = z_info.cost_to_come + c
                if cost_through_z < x_info.cost_to_come
                    node_info[x] = (x_info..., status='O', parent=z, cost_to_come=cost_through_z)
                    open_queue[x] = cost_through_z
                end
            end
        end
        node_info[z] = (node_info[z]..., status='C')
        isempty(open_queue) ? break : z = dequeue!(open_queue)
    end

    metadata[:solved] = graph[z] in goal
    record_solution!(metadata, node_info, z, X(0))
    nothing
end
