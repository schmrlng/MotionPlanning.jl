export shortcut, adaptive_shortcut, adaptive_shortcut!
export smooth_solution!, time_discretize_solution!, time_space_solution!

### ADAPTIVE-SHORTCUT (Hsu 2000)

function shortcut(path::Path, CC::CollisionChecker)
    N = length(path)
    if N == 2
        return path
    end
    if is_free_motion(path[1], path[end], CC)
        return path[[1,end]]
    end
    mid = ceil(Int, N/2)
    return [shortcut(path[1:mid], CC)[1:end-1]; shortcut(path[mid:end], CC)]
end

function cut_corner(v1::State, v2::State, v3::State, CC::CollisionChecker)
    m1 = (v1 + v2)/2
    m2 = (v3 + v2)/2
    while !is_free_motion(m1, m2, CC)
        m1 = (m1 + v2)/2
        m2 = (m2 + v2)/2
    end
    return typeof(v1)[v1, m1, m2, v3]
end

function adaptive_shortcut(path::Path, CC::CollisionChecker, iterations::Int = 10)
    while (short_path = shortcut(path, CC)) != path
        path = short_path
    end
    for i in 1:iterations
        path = [path[1:1]; vcat([cut_corner(path[j-1:j+1]..., CC)[2:3] for j in 2:length(path)-1]...); path[end:end]]
        while (short_path = shortcut(path, CC)) != path
            path = short_path
        end
    end
    return path, cumsum([0; map(norm, diff(path))])
end

function adaptive_shortcut!(P::MPProblem, iterations::Int = 10)
    P.status == :solved || error("Cannot post-process unsolved problem! (adaptive-shortcut)")
    isa(P.SS.dist, Euclidean) || error("Adaptive-shortcut requires Euclidean SS")
    S = P.solution
    smoothed_path, smoothed_cumcost = adaptive_shortcut(P.V.V[S.metadata["path"]], P.CC, iterations)
    S.metadata["smoothed_path"] = smoothed_path
    S.metadata["smoothed_cumcost"] = smoothed_cumcost
    S.metadata["smoothed_cost"] = smoothed_cumcost[end]
    smoothed_cumcost[end]
end

### Smoothing

function smooth_solution!(P::MPProblem)
    P.status != :solved && error("Cannot post-process unsolved problem! (adaptive-shortcut)")
    isa(P.SS.dist, Euclidean) && return adaptive_shortcut!(P)
end

### Path discretization

function time_discretize_solution!(P::MPProblem, dt::Real, usesmoothed = true)
    if usesmoothed && haskey(P.solution.metadata, "smoothed_path")
        path = P.solution.metadata["smoothed_path"]
        costs = P.solution.metadata["smoothed_cumcost"]
    else
        path = P.V[P.solution.metadata["path"]]
        costs = P.solution.metadata["cumcost"]
    end
    controlseq = steering_control(P.SS, path...)
    P.solution.metadata["discretized_path"] = propagate(P.SS, path[1], controlseq, 0:dt:duration(controlseq))
end

function time_space_solution!(P::MPProblem, n::Integer, usesmoothed = true)
    if usesmoothed && haskey(P.solution.metadata, "smoothed_path")
        path = P.solution.metadata["smoothed_path"]
        costs = P.solution.metadata["smoothed_cumcost"]
    else
        path = P.V[P.solution.metadata["path"]]
        costs = P.solution.metadata["cumcost"]
    end
    controlseq = steering_control(P.SS, path...)
    P.solution.metadata["discretized_path"] = propagate(P.SS, path[1], controlseq, linspace(0., duration(controlseq), n))
end
