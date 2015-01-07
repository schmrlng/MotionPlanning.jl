using Graphs

function fmtstar(P::ProblemSetup, N::Int, rm::Float64)
    tic()
    collision_checks = 0

    V, free_volume_ub = sample_free(P, N)
    NN = MetricNNKDTree(V, P.SS.dist)

    A = zeros(Int,N)
    W = trues(N)
    W[1] = false
    H = falses(N)
    H[1] = true
    C = zeros(N)
    z = 1

    k = min(iceil((2*rm)^P.SS.dim*(e/P.SS.dim)*log(N)), N-1)

    while ~is_goal_pt(V[z], P.goal)
        H_new = Int[]
        for x in mutualnearestk(NN, z, k, W).inds
            neighborhood = nearestk(NN, x, k, H)
            c_min, y_idx = findmin(C[neighborhood.inds] + neighborhood.ds)
            y_min = neighborhood.inds[y_idx]
            if (collision_checks = collision_checks + 1; is_free_motion(V[x], V[y_min], P.obs, P.SS))
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
        end
    end

    sol = [z]
    while sol[1] != 1
        unshift!(sol, A[sol[1]])
    end

    solution_metadata = {
        "radius_multiplier" => rm,
        "collision_checks" => collision_checks,
        "elapsed" => toq(),
        "num_samples" => N,
        "cost" => C[z],
        "test_config" => P.config_name,
        "planner" => "FMTstar",
        "solved" => is_goal_pt(V[z], P.goal)
    }

    return solution_metadata, sol, V, A
end