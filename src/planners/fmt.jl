export fmtstar!

fmtstar!(P::MPProblem; kwargs...) = fmtstar!(P, length(P.V); kwargs...)
function fmtstar!(P::MPProblem, N::Int; rm::Float64 = 1.0,
                                        connections::Symbol = :R,
                                        k = min(ceil(Int, (2*rm)^P.SS.dim*(e/P.SS.dim)*log(N)), N-1),
                                        r = 0.,
                                        ensure_goal_ct = 1,
                                        init_idx = 1,
                                        checkpts = true)  # TODO: bleh, prefer false
    tic()
    P.CC.count = 0

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
    if !is_free_state(P.init, P.CC, P.SS)
        warn("Initial state is infeasible!")
        P.status = :failed
        P.solution = MPSolution(P.status, Inf, toq(), Dict())
        return Inf
    end
    free_volume_ub = sample_free!(P, N - length(P.V), ensure_goal_ct = ensure_goal_ct)  # TODO: clean this logic up
    if checkpts
        F = trues(N)
        for i in 1:N
            F[i] = is_free_state(P.V[i], P.CC, P.SS)
        end
    end
    if r == 0.             # TODO: this default connection radius only really applies for geometric planning
        dim = P.SS.dim     # in particular, state space dimension here should come from the metric instead
        r = rm*2*(1/dim*free_volume_ub/(pi^(dim/2)/gamma(dim/2+1))*log(N)/N)^(1/dim)
        setup_steering(P.SS, r)
    end

    A = zeros(Int,N)
    W = trues(N)
    H = falses(N)
    C = zeros(Float64,N)
    P.V.init = P.init
    if P.V[init_idx] == P.init
        W[init_idx] = false
        H[init_idx] = true
        HHeap = Collections.PriorityQueue([init_idx], [0.])
    else    # special casing the first expansion round of FMT if P.init is not in the sample set
        HHeap = Collections.PriorityQueue(Int[], Float64[])
        neighborhood = (connections == :R ? inballF(P.V, P.init, r) : knnF(P.V, P.init, r))
        for ii in 1:length(neighborhood.inds)
            x, c = neighborhood.inds[ii], neighborhood.ds[ii]
            if is_free_motion(P.init, P.V[x], P.CC, P.SS)
                A[x] = 0
                C[x] = c
                HHeap[x] = c
                H[x] = true
                W[x] = false
            end
        end
    end
    z = Collections.dequeue!(HHeap)    # i.e. z = init_idx

    while !is_goal_pt(P.V[z], P.goal, P.SS)
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
            z = Collections.dequeue!(HHeap)
        else
            break
        end
    end

    sol = [z]
    costs = [C[z]]
    while sol[1] != 1
        unshift!(sol, A[sol[1]])
        if sol[1] == 0
            unshift!(costs, 0.)
            break
        end
        unshift!(costs, C[sol[1]])
    end

    P.status = is_goal_pt(P.V[z], P.goal, P.SS) ? :solved : :failed
    solution_metadata = Dict(
        "radius_multiplier" => rm,
        "collision_checks" => P.CC.count,
        "num_samples" => N,
        "cost" => C[z],
        "cumcost" => costs,
        "planner" => "FMTstar",
        "solved" => is_goal_pt(P.V[z], P.goal, P.SS),
        "tree" => A,
        "path" => sol
    )
    connections == :R && (solution_metadata["r"] = r)
    connections == :K && (solution_metadata["k"] = k)
    P.solution = MPSolution(P.status, C[z], toq(), solution_metadata)
    C[z]
end