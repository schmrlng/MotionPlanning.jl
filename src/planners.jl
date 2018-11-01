function standard_setup!(P)
    tic = time()
    reset!(P.collision_checker)
    P.status = is_free_state(P.collision_checker, P.init) ? :inprogress :
                                                            (@warn "Initial state is infeasible!"; :infeasible)
    P.solution = MPSolution(P.status, nothing, tic, Dict{Symbol,Any}())
    nothing
end
macro standard_setup!(P)
    P = :($(esc(P)))
    quote
        standard_setup!($P)
        $P.status === :infeasible && return $P.solution
        $P.solution.metadata
    end
end

function backtrace_tree_path(node_info, z, terminal_index)
    path_nodes = [z]
    costs_to_come  = [node_info[z].cost_to_come]
    while path_nodes[1] != terminal_index
        pushfirst!(path_nodes, node_info[path_nodes[1]].parent)
        pushfirst!(costs_to_come,  node_info[path_nodes[1]].cost_to_come)
    end
    path_nodes, costs_to_come
end
function record_tree_solution!(metadata, node_info, z, terminal_index)
    metadata[:node_info] = node_info
    if metadata[:solved]
        metadata[:solution_nodes], metadata[:costs_to_come] = backtrace_tree_path(node_info, z, terminal_index)
        metadata[:cost] = metadata[:costs_to_come][end]
    end
end

function record_tree!(metadata, node_info)
    metadata[:tree] = Dict(k => node_info[k].parent for k in keys(node_info))
end

function standard_wrapup!(P)
    metadata = P.solution.metadata
    metadata[:motion_checks] = P.collision_checker.motion_count[]
    metadata[:edge_checks] = P.collision_checker.edge_count[]
    P.status = metadata[:solved] ? :solved : :failed
    P.solution = MPSolution(P.status, get(metadata, :cost, missing), time() - P.solution.elapsed, metadata)
end

include("planners/fmtstar.jl")
# include("planners/rrt.jl")
# include("planners/rrtstar.jl")
