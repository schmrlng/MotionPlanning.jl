export fmtstar!

fmtstar!(P::MPProblem; kwargs...) = fmtstar!(P, length(P.V); kwargs...)
function fmtstar!(P::MPProblem, N::Int; rm::Float64 = 1.0,
                                        connections::Symbol = :R,
                                        k = min(iceil((2*rm)^P.SS.dim*(e/P.SS.dim)*log(N)), N-1),
                                        r = 0.,
                                        checkpts = true)  # TODO: bleh, prefer false
    tic()
    P.CC.count = 0

    # TODO: staged functions (Julia v0.4) for knn vs ball... or something clever in v0.3
    if connections == :R
        nearF = inballF!
        nearB = inballB!
    elseif connections == :K
        nearF = mutualknnF!
        nearB = knnB!
    else
        error("Connection type must be radial (:R) or k-nearest (:K)")
    end
    r > 0 && setup_steering(P.SS, r)
    free_volume_ub = sample_free!(P, N - length(P.V))  # TODO: clean this logic up
    if checkpts
        F = trues(N)
        for i in 1:N
            F[i] = is_free_state(P.V[i], P.CC, P.SS)
        end
    end
    dim = P.SS.dim
    if r == 0.
        r = rm*2*(1/dim*free_volume_ub/(pi^(dim/2)/gamma(dim/2+1))*log(N)/N)^(1/dim)
        setup_steering(P.SS, r)
    end

    A = zeros(Int,N)
    W = trues(N); W[1] = false
    H = falses(N); H[1] = true
    C = zeros(Float64,N)
    HHeap = CollectionsJ4.PriorityQueue([1], [0.])
    z = CollectionsJ4.dequeue!(HHeap)    # i.e. z = 1

    while ~is_goal_pt(P.V[z], P.goal)
        H_new = Int[]
        for x in (connections == :R ? nearF(P.V, z, r, W).inds : nearF(P.V, z, k, W).inds)
            checkpts && !F[x] && continue
            neighborhood = (connections == :R ? nearB(P.V, x, r, H) : nearB(P.V, x, k, H))
            c_min, y_idx = findmin(C[neighborhood.inds] + neighborhood.ds)
            y_min = neighborhood.inds[y_idx]
            if is_free_motion(P.V[y_min], P.V[x], P.CC, P.SS)
                A[x] = y_min
                C[x] = c_min
                HHeap[x] = c_min
                push!(H_new, x)
                W[x] = false
            end
        end
        H[H_new] = true
        H[z] = false
        if !isempty(HHeap)
            z = CollectionsJ4.dequeue!(HHeap)
        else
            break
        end
    end

    sol = [z]
    costs = [C[z]]
    while sol[1] != 1
        unshift!(sol, A[sol[1]])
        unshift!(costs, C[sol[1]])
    end

    P.status = is_goal_pt(P.V[z], P.goal) ? :solved : :failed
    solution_metadata = {
        "radius_multiplier" => rm,
        "collision_checks" => P.CC.count,
        "num_samples" => N,
        "cost" => C[z],
        "cumcost" => costs,
        "planner" => "FMTstar",
        "solved" => is_goal_pt(P.V[z], P.goal),
        "tree" => A,
        "path" => sol
    }
    connections == :R && (solution_metadata["r"] = r)
    connections == :K && (solution_metadata["k"] = k)
    P.solution = MPSolution(P.status, C[z], toq(), solution_metadata)
    C[z]
end