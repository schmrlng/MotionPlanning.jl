export shortcut, adaptive_shortcut, adaptive_shortcut!, discretize_path
export cost_discretize_solution!, cost_space_solution!, time_discretize_solution!

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
    P.status != :solved && error("Cannot post-process unsolved problem! (adaptive-shortcut)")
    (!isa(P.SS, RealVectorStateSpace) || P.SS.dist != Euclidean()) && error("Adaptive-shortcut requires Euclidean SS")
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
    isa(P.SS, RealVectorStateSpace) && isa(P.SS.dist, Euclidean) && return adaptive_shortcut!(P)
    # isa(P.SS, SE2StateSpace) && isa(P.SS.dist, SimpleCarMetric) && return 
end


### Path discretization (Euclidean)

# function discretize_path(path0::Path, dx)
#     path = path0[1:1]
#     for i in 2:length(path0)
#          norm(path0[i] - path[end]) > 2dx/3 && push!(path, path0[i])  # cut out small waypoint steps
#     end
#     dpath = path[1:1]
#     for i in 2:length(path)
#         segment_length = norm(path[i] - path[i-1])
#         M = ceil(segment_length / dx)
#         append!(dpath, [path[i-1] + (j/M)*(path[i] - path[i-1]) for j in 1:M])
#     end
#     map(dense, dpath)
# end

### Reeds-Shepp fix inexact steering (okay not really) and discretize

# function smooth_waypoints{T}(v::RSState{T}, s::RSSegment{T}, r::T, dx)
#     s.t == 0 && return Vector2{T}[v.x + a*Vector2(cos(v.t), sin(v.t)) for a in linspace(0, abs(s.d), iceil(abs(s.d)/dx)+1)]
#     center = v.x + sign(s.t)*Vector2(-r*sin(v.t), r*cos(v.t))
#     turnpts = [r*Vector2(cos(x), sin(x)) for x in linspace(0, abs(s.d), iceil(r*abs(s.d)/dx)+1)]
#     if s.t*s.d < 0
#         for i in 1:length(turnpts)      # manual @devec
#             turnpts[i] = Vector2(turnpts[i][1], -turnpts[i][2])
#         end
#     end
#     [(center + sign(s.t)*rotate(p, v.t-pi/2)) for p in turnpts]
# end
# function smooth_waypoints{T}(v::RSState{T}, w::RSState{T}, SS::ReedsSheppStateSpace, dx)
#     pts = Array(Vector2{T}, 0)
#     for s in SS.dist.paths[RSvec2sub(v, w, SS.dist)...]
#         s_pts = smooth_waypoints(v, s, SS.r, dx)
#         append!(pts, s_pts[1:end-1])
#         v = RSState(s_pts[end], v.t + s.t*s.d)
#     end
#     scale_factor = (w.x - pts[1]) ./ (v.x - pts[1])
#     pts = [(scale_factor.*(p - pts[1]) + pts[1]) for p in pts]
#     push!(pts, w.x)
# end
# smooth_waypoints(i::Int, j::Int, NN::NearNeighborCache, SS::ReedsSheppStateSpace, dx) = smooth_waypoints(NN[i], NN[j], SS, dx)

# function RSsmooth_and_discretize!(P::MPProblem, dx)
#     (!isa(P.SS, ReedsSheppStateSpace) || P.status != :solved) && error("RSsmooth_and_discretize only works for solved Reeds-Shepp problems!")
#     sol = P.solution.metadata["path"]
#     dpath = vcat([smooth_waypoints(sol[i], sol[i+1], P.V, P.SS, dx)[1:end-1] for i in 1:length(sol)-1]...)
#     push!(dpath, P.V[sol[end]].x)
#     dpath = dpath[[true, map(norm, diff(dpath)) .> dx / 4]] # cut out tiny waypoint steps
#     P.solution.metadata["smoothed_path"] = dpath
# end

function time_discretize_solution!(P::MPProblem, dt)
    if haskey(P.solution.metadata, "smoothed_path")
        path = P.solution.metadata["smoothed_path"]
        costs = P.solution.metadata["smoothed_cumcost"]
    else
        path = P.V[P.solution.metadata["path"]]
        costs = P.solution.metadata["cumcost"]
    end
    controlseq = steering_control(P.SS, path...)
    P.solution.metadata["discretized_path"] = propagate(P.SS, path[1], controlseq, 0:dt:duration(controlseq))
end

function time_space_solution!(P::MPProblem, n)
    if haskey(P.solution.metadata, "smoothed_path")
        path = P.solution.metadata["smoothed_path"]
        costs = P.solution.metadata["smoothed_cumcost"]
    else
        path = P.V[P.solution.metadata["path"]]
        costs = P.solution.metadata["cumcost"]
    end
    controlseq = steering_control(P.SS, path...)
    P.solution.metadata["discretized_path"] = propagate(P.SS, path[1], controlseq, linspace(0., duration(controlseq), n))
end

# function cost_discretize_solution!(P::MPProblem, dc)  # TODO: write time_discretize_solution! after standardizing StateSpace definitions
#     if haskey(P.solution.metadata, "smoothed_path")
#         path = P.solution.metadata["smoothed_path"]
#         costs = P.solution.metadata["smoothed_cumcost"]
#     else
#         path = P.V[P.solution.metadata["path"]]
#         costs = P.solution.metadata["cumcost"]
#     end
#     x0 = path[1]
#     dpath = typeof(x0)[x0]
#     c = dc
#     for i in 1:length(path)-1
#         while c < costs[i+1]
#             push!(dpath, cost_waypoint(P.SS, path[i], path[i+1], c - costs[i]))
#             c = c + dc
#         end
#     end
#     push!(dpath, path[end])
#     P.solution.metadata["discretized_path"] = dpath
# end

# function cost_space_solution!(P::MPProblem, n)
#     for k in n-1:n+10
#         cost_discretize_solution!(P, P.solution.cost / k)
#         if length(P.solution.metadata["discretized_path"]) >= n
#             P.solution.metadata["discretized_path"] = P.solution.metadata["discretized_path"][1:n]
#             return P.solution.metadata["discretized_path"]
#         end
#     end
#     error("Something must be seriously wrong with costs to get here...")
# end

# time_discretize_solution!(P::MPProblem, dt) = cost_discretize_solution!(P::MPProblem, dt) # TODO: TEMP

### Linear Quadratic Discretization

# function discretize_path(pidx, dt, NN::NearNeighborCache, SS::LinearQuadraticStateSpace)
#     dpath = NN.V[[pidx[1]]]
#     for i in 2:length(pidx)
#         segment_length = NN.US[pidx[i-1],pidx[i]].t
#         M = iceil(segment_length / dt) + 1
#         wps = statepoints(pidx[i-1], pidx[i], NN, SS, M)
#         append!(dpath, wps[2:end])
#     end
#     dpath
# end

### General (should perhaps not live in this package?)

# function discretize_path(P::MPProblem, dt)
#     if P.SS.dist == Euclidean()
#         adaptive_shortcut!(P)
#         P.solution.metadata["discretized_path"] = discretize_path(P.solution.metadata["smoothed_path"], dt)
#         return P.solution.metadata["discretized_path"]   # TODO: get method for MPSolution type
#     end
#     P.solution.metadata["discretized_path"] = discretize_path(P.solution.metadata["path"], dt, P.V, P.SS)
# end