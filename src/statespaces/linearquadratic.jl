immutable LinearQuadraticStateSpace <: DifferentialStateSpace
    dim::Int
    lo::Vector
    hi::Vector

    A::Matrix
    B::Matrix
    c::Vector
    R::Matrix
    G::Function
    Ginv::Function
    expAt::Function
    cdrift::Function

    xbar::Function
    cost::Function
    cost_deriv::Function
    x::Function
end

immutable QuasiMetricProblem <: ProblemSetup
    init::Vector{Float64}
    goal::Goal
    obs::ObstacleSet
    V0::Vector{Vector{Float64}}
    SS::LinearQuadraticStateSpace
    config_name::String
end

immutable QuasiMetricNN <: NearNeighborCache
    V::Vector{Vector{Float64}}
    D::Matrix
    TT::Matrix
    NNF::Vector{Vector{Int64}}
    NNB::Vector{Vector{Int64}}
end

function LinearQuadraticStateSpace(dim, lo, hi, A, B, c, R, G, Ginv, expAt, cdrift)
    xbar(x0, t) = expAt(t)*x0 + cdrift(t)
    cost(x0, x1, t) = ( t + (x1 - xbar(x0,t))'*Ginv(t)*(x1 - xbar(x0,t)) )[1]
    function cost_deriv(x0, x1, t)
        d = Ginv(t)*(x1 - xbar(x0,t))
        return ( 1 - 2*(A*x1 + c)'*d - d'*B*inv(R)*B'*d )[1]
    end
    x(x0, x1, t, s) = xbar(x0,s) + G(s)*expAt(t-s)'*Ginv(t)*(x1 - xbar(x0,t))

    return LinearQuadraticStateSpace(dim, lo, hi, A, B, c, R, G, Ginv, expAt, cdrift, xbar, cost, cost_deriv, x)
end

volume(SS::LinearQuadraticStateSpace) = prod(SS.hi-SS.lo)

function pairwise_distances{T}(V::Vector{Vector{T}}, SS::LinearQuadraticStateSpace, t_bound::Float64)
    return pairwise_distances_approx_opt(V, SS, t_bound)
end

function pairwise_distances_approx_opt{T}(V::Vector{Vector{T}}, SS::LinearQuadraticStateSpace, t_bound::Float64, res::Int64 = 10)
    N = length(V)
    V1 = hcat(V...)
    t = t_bound/2.
    V0bar = SS.expAt(t)*V1 .+ SS.cdrift(t)
    return (t + pairwise(SqMahalanobis(SS.Ginv(t)), V0bar, V1)), fill(t, N, N)
end

# function pairwise_distances_naive_opt{T}(V::Vector{Vector{T}}, SS::LinearQuadraticStateSpace, t_bound::Float64)
#     N = length(V)
#     D = fill(Inf, N, N)
#     TT = fill(Inf, N, N)
#     for j = 1:N
#         for i = 1:N
#             if i != j
#                 t_opt = try fzero(t -> SS.cost_deriv(V[i],V[j],t), 0.001, t_bound) catch; continue end
#                 @inbounds D[i,j] = SS.cost(V[i], V[j], t_opt)
#                 @inbounds TT[i,j] = t_opt
#             end
#         end
#     end
#     return D, TT
# end

function NNCache{T}(V::Vector{Vector{T}}, SS::LinearQuadraticStateSpace, t_bound::Float64)        # consider pre-filling NNF, NNB here?
    return QuasiMetricNN(V, pairwise_distances(V, SS, t_bound)..., fill(Int64[], length(V)), fill(Int64[], length(V)))
end

function nearRF(NN::QuasiMetricNN, v::Int64, r::Float64)
    if isempty(NN.NNF[v])
        nn_bool = NN.D[v,:] .< r
        nn_bool[v] = false
        nn_idx = find(nn_bool)
        NN.NNF[v] = nn_idx
    else
        nn_idx = NN.NNF[v]
    end
    return nn_idx, vec(NN.D[v, nn_idx])
end

function nearRB(NN::QuasiMetricNN, v::Int64, r::Float64)
    if isempty(NN.NNB[v])
        nn_bool = NN.D[:,v] .< r
        nn_bool[v] = false
        nn_idx = find(nn_bool)
        NN.NNB[v] = nn_idx
    else
        nn_idx = NN.NNB[v]
    end
    return nn_idx, NN.D[nn_idx, v]
end

function nearRF(NN::NearNeighborCache, v::Int64, r::Float64, filter::BitVector)
    nn_idx, D = nearRF(NN, v, r)
    return nn_idx[filter[nn_idx]], D[filter[nn_idx]]
end

function nearRB(NN::NearNeighborCache, v::Int64, r::Float64, filter::BitVector)
    nn_idx, D = nearRB(NN, v, r)
    return nn_idx[filter[nn_idx]], D[filter[nn_idx]]
end

function waypoints(i::Int64, j::Int64, NN::QuasiMetricNN, SS::LinearQuadraticStateSpace, res=5)       # make arg ordering consistent NN to front
    return hcat([SS.x(NN.V[i], NN.V[j], NN.TT[i,j], s) for s in linspace(0, NN.TT[i,j], res)]...)
end

function is_free_motion(i::Int64, j::Int64, NN::QuasiMetricNN, obs::ObstacleSet, SS::LinearQuadraticStateSpace)
    wps = waypoints(i, j, NN, SS)
    return all(SS.lo .<= wps .<= SS.hi) && is_free_path(Vector{Float64}[wps[:,c] for c in 1:size(wps,2)], obs)
end

function plot_tree(NN::QuasiMetricNN, A, SS::LinearQuadraticStateSpace, col="black"; kwargs...)
    pts = hcat(NN.V[find(A)]...)
    scatter(pts[1,:], pts[2,:], zorder=1; kwargs...)
    X = vcat([[waypoints(A[v], v, NN, SS, 20)[1,:]', nothing] for v in find(A)]...)
    Y = vcat([[waypoints(A[v], v, NN, SS, 20)[2,:]', nothing] for v in find(A)]...)
    # @show X
    plt.plot(X, Y, color=col, linewidth=.5, linestyle="-", zorder=1; kwargs...)
end

function plot_solution(NN::QuasiMetricNN, sol, SS::LinearQuadraticStateSpace, disc = 20, col="black"; kwargs...)
    pts = hcat(NN.V[sol]...)
    scatter(pts[1,:], pts[2,:], zorder=2, color=col, s=25; kwargs...)
    if length(sol) > 1
        wps = hcat([waypoints(sol[i], sol[i+1], NN, SS, disc) for i in 1:length(sol)-1]...)
        plot_path(wps, color=col, linewidth=2.0; kwargs...)
        plt.quiver([wps[row,1:3:end]' for row in 1:4]..., zorder=5, width=.003, headwidth=8)
        # print(maximum(abs(wps),2))
    end
end

function fmtstar(P::QuasiMetricProblem, N::Int64, r::Float64)
    tic()
    collision_checks = 0

    V, free_volume_ub = sample_free(P, N)
    NN = NNCache(V, P.SS, r)

    tic()
    A = zeros(Int64,N)
    W = trues(N)
    W[1] = false
    H = falses(N)
    H[1] = true
    C = zeros(N)
    z = 1

    while ~is_goal_pt(V[z], P.goal)
        H_new = Int64[]
        for x in nearRF(NN, z, r, W)[1]
            Y_near, D_near = nearRB(NN, x, r, H)
            c_min, y_idx = findmin(C[Y_near] + D_near)
            y_min = Y_near[y_idx]
            if (collision_checks = collision_checks + 1; is_free_motion(y_min, x, NN, P.obs, P.SS))   # doesn't matter for other algs, but should make "direction" consistent elsewhere
                A[x] = y_min
                C[x] = c_min
                push!(H_new, x)
                W[x] = false
            end
        end
        H[H_new] = true
        H[z] = false
        if any(H)
            z = find(H)[indmin(C[H])]
        else
            break
            # return Inf, Int64[], V, A
        end
    end

    sol = [z]
    while sol[1] != 1
        unshift!(sol, A[sol[1]])
    end

    runtime_elapsed = toq()
    tot_elapsed = toq()

    solution_metadata = {
        "radius_multiplier" => rm,
        "collision_checks" => collision_checks,
        "elapsed" => runtime_elapsed,
        "total_elapsed" => tot_elapsed,
        "num_samples" => N,
        "cost" => is_goal_pt(V[z], P.goal) ? C[z] : Inf,
        "test_config" => P.config_name,
        "planner" => "FMT",
        "solved" => is_goal_pt(V[z], P.goal)
    }

    return solution_metadata, sol, NN, A
end