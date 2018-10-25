export fmtstar!

# status = 'U', 'O', or 'C' (unvisited, open, closed)
const FMTNodeInfo{X,D} = NamedTuple{(:is_free, :status, :parent, :cost_to_come),Tuple{Bool,Char,X,D}}
function FMTNodeInfo{X,D}(; is_free=false, status='U', parent=zero(X), cost_to_come=zero(D)) where {X,D}
    (is_free=is_free, status=status, parent=parent, cost_to_come=cost_to_come)
end

function fmtstar!(P::MPProblem; r, ensure_goal_count=0, compute_full_metadata=true)
    # Graph Node Setup
    !isdefined(P, :graph) && error("TODO")
    setradius!(P.graph, r)
    setinit!(P.graph, P.init)
    setgoal!(P.graph, () -> rand_free_state(P.collision_checker, P.goal), sample_count=ensure_goal_count)
    compute_init_and_goal_neighbors!(P.graph)

    # Solve
    P.status, P.solution = fmtstar!(P.state_space, P.bvp, P.init, P.goal, P.collision_checker,
                                    P.graph, r)

    # Post-Processing
    if P.status === :solved && compute_full_metadata
        metadata = P.solution.metadata

        X = indextype(P.graph)
        node_info = metadata[:node_info]
        metadata[:tree] = Dict{X,X}(k => node_info[k].parent for k in keys(node_info) if node_info[k].status != 'U')
    end
    P.status, P.solution
end

function fmtstar!(state_space::StateSpace,
                  bvp::SteeringBVP,
                  init::State,
                  goal::Goal,
                  collision_checker::CollisionChecker,
                  graph::NearNeighborGraph{NeighborInfo{X,D,U}},
                  r,
                  node_info=node_info_datastructure(
                      graph.nodes, FMTNodeInfo{X,D},
                      x -> FMTNodeInfo{X,D}(is_free=is_free_state(collision_checker, x)),
                      bounds=boundingbox(state_space))
                  ) where {X,D,U}
    tic = time()
    collision_checker.motion_count[] = collision_checker.edge_count[] = 0    # reset!(collision_checker)

    if !is_free_state(collision_checker, init)
        # @warn("Initial state is infeasible!")
        status = :infeasible
        solution = MPSolution(status, D(Inf), time() - tic, Dict{Symbol,Any}())
        return status, solution
    end

    is_unvisited = n -> (ni = node_info[n.index]; ni !== missing && ni.is_free && ni.status == 'U')
    is_open      = n -> (ni = node_info[n.index]; ni !== missing && ni.status == 'O')

    open_queue = PriorityQueue{X,D}()
    z = X(0)
    node_info[z] = (node_info[z]..., status='O')

    while !(graph[z] in goal)
        for (x, _, _) in neighbors(is_unvisited, graph, z, dir=Val(:F))
            # @show (z, x)
            y_min, c_min, u_min = mapreduce(n -> (index=n.index,
                                                  cost=node_info[n.index].cost_to_come + n.cost,
                                                  controls=n.controls),
                                            (a, b) -> ifelse(a.cost < b.cost, a, b),
                                            neighbors(is_open, graph, x, dir=Val(:B)))
            if is_free_edge(collision_checker, bvp, graph[y_min], graph[x], u_min)
                node_info[x] = (is_free=true, status='O', parent=y_min, cost_to_come=c_min)
                open_queue[x] = c_min
            end
        end
        node_info[z] = (node_info[z]..., status='C')
        isempty(open_queue) ? break : z = dequeue!(open_queue)
    end

    solution_nodes = [z]
    costs_to_come = [node_info[z].cost_to_come]
    while linear_index(solution_nodes[1]) != 0
        pushfirst!(solution_nodes, node_info[solution_nodes[1]].parent)
        pushfirst!(costs_to_come, node_info[solution_nodes[1]].cost_to_come)
    end

    status = graph[z] in goal ? :solved : :failed
    solution_metadata = Dict{Symbol,Any}(
        :r => r,
        :motion_checks => collision_checker.motion_count[],
        :edge_checks => collision_checker.edge_count[],
        :cost => costs_to_come[end],
        :costs_to_come => costs_to_come,
        :planner => :FMTstar,
        :solved => status === :solved,
        :solution_nodes => solution_nodes,
        :node_info => node_info
    )
    solution = MPSolution(status, costs_to_come[end], time() - tic, solution_metadata)
    status, solution
end
