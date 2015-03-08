export MPProblem, MPSolution

type MPSolution
    status::Symbol
    cost::Float64
    elapsed::Float64
    metadata::Dict
end

type MPProblem
    SS::StateSpace
    init::State
    goal::Goal
    CC::CollisionChecker
    V::NearNeighborCache
    config_name::String
    status::Symbol
    solution::MPSolution

    function MPProblem(SS::StateSpace,
                       init::State,
                       goal::Goal,
                       CC::CollisionChecker,
                       V::NearNeighborCache,
                       config_name::String="$(SS.dim)D $(typeof(SS))")
        new(SS, init, goal, CC, V, config_name, "not yet solved")
    end
end

MPProblem(SS::StateSpace, init::State, goal::Goal, CC::CollisionChecker) = MPProblem(SS, init, goal, CC, defaultNN(SS, init))

function plot(P::MPProblem; SS=true, CC=true, goal=true, meta=false, sol=true)
    SS && plot(P.SS)
    CC && plot(P.CC)
    goal && plot(P.goal)
    if isdefined(P, :solution)
        S = P.solution
        if meta
            haskey(S.metadata, "tree") && plot_tree(P.V.V, S.metadata["tree"], color="gray", alpha=0.5)
            # TODO: graph (PRM)
        end
        sol && plot_path(P.V.V, S.metadata["path"], color="blue")
    end
end