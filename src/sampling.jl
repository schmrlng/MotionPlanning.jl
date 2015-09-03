export sample_free_goal, sample_free!

function sample_free_goal(P::MPProblem)
    v = sample_goal(P.goal, P.SS)
    while !is_free_state(v, P.CC, P.SS)
        v = sample_goal(P.goal, P.SS)
    end
    return v
end

function sample_free!(P::MPProblem, N::Integer, ensure_goal::Bool = true, goal_bias = 0.0)
    N == 0 && return volume(P.SS)
    V = P.V.V
    W = Array(eltype(V), N)
    if length(V) > 0 && V[1] == P.init
        sample_count = 0
    else
        sample_count = 1
        W[1] = P.init
    end
    attempts = 0
    successes = 0
    for v in repeatedly(() -> sample_space(P.SS))
        attempts += 1
        if is_free_state(v, P.CC, P.SS)
            sample_count += 1
            successes += 1
            if 0.0 < goal_bias && rand() < goal_bias
                W[sample_count] = sample_free_goal(P)
            else
                W[sample_count] = v
            end
            if sample_count == N
                break
            end
        end
    end
    ensure_goal && (W[N] = sample_free_goal(P))
    P.V = addpoints(P.V, W)
    return volume(P.SS) #ci(BinomialTest(successes, attempts), .05, tail = :left)[2] * volume(P.SS)
end