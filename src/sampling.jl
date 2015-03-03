function sample_free_goal(P::MPProblem)
    v = vector_to_state(sample_goal(P.goal), P.SS)
    while ~is_free_state(v, P.obs, P.SS)
        v = vector_to_state(sample_goal(P.goal), P.SS)
    end
    return v
end

function sample_free!(P::MPProblem, N::Integer, ensure_goal::Bool = true, goal_bias::FloatingPoint = 0.0)
    V = P.V.V
    length(V) > 0 && V[0] != P.init
    fill(P.init, N)
    sample_count = 1
    attempts = 0
    for v in chain(P.V0, iterate(x -> sample_space(P.SS), sample_space(P.SS))) # TODO: Iterators.repeatedly
        attempts += 1
        if is_free_state(v, P.obs, P.SS)
            sample_count = sample_count + 1
            if 0.0 < goal_bias && rand() < goal_bias
                V[sample_count] = sample_free_goal(P)
            else
                V[sample_count] = v
            end
            if sample_count == N
                break
            end
        end
    end
    if ensure_goal && ~any([is_goal_pt(v, P.goal) for v in V])
        V[N] = sample_free_goal(P)
    end
    return V, 1. #ci(BinomialTest(sample_count, attempts), .05, tail = :left)[2] * volume(P.SS)
end