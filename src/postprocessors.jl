export shortcut, adaptive_shortcut, adaptive_shortcut!

### ADAPTIVE-SHORTCUT (Hsu 2000)

function shortcut(path::Path, CC::CollisionChecker)
    N = length(path)
    if N == 2
        return path
    end
    if is_free_motion(path[1], path[end], CC)
        return path[[1,end]]
    end
    mid = iceil(N/2)
    return [shortcut(path[1:mid], CC)[1:end-1], shortcut(path[mid:end], CC)]
end

function cut_corner(v1::AbstractVector, v2::AbstractVector, v3::AbstractVector, CC::CollisionChecker)
    m1 = (v1 + v2)/2
    m2 = (v3 + v2)/2
    while ~is_free_motion(m1, m2, CC)
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
        path = [path[1:1], vcat([cut_corner(path[j-1:j+1]..., CC)[2:3] for j in 2:length(path)-1]...), path[end:end]]
        while (short_path = shortcut(path, CC)) != path
            path = short_path
        end
    end
    return path, sum(mapslices(norm, diff(hcat(path...), 2), 1))
end

function adaptive_shortcut!(P::MPProblem)
    P.status != :solved && error("Cannot post-process unsovled problem! (adaptive-shortcut)")
    (!isa(P.SS, RealVectorMetricSpace) || P.SS.dist != Euclidean()) && error("Adaptive-shortcut requires Euclidean SS")
    S = P.solution
    smoothed_path, smoothed_cost = adaptive_shortcut(P.V.V[S.metadata["path"]], P.CC)
    S.metadata["smoothed_path"] = smoothed_path
    S.metadata["smoothed_cost"] = smoothed_cost
    smoothed_cost
end