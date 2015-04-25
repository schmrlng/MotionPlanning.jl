import Base.copy
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
function copy(P::MPProblem)
    Pcopy = MPProblem(P.SS, P.init, P.goal, P.CC, P.V, P.config_name)
    Pcopy.status = P.status
    Pcopy.solution = P.solution
    Pcopy
end

function plot(P::MPProblem; SS=true, CC=true, goal=true, meta=false, sol=true, smoothed=false)
    SS && plot(P.SS)
    CC && plot(P.CC, P.SS.lo, P.SS.hi)
    goal && plot(P.goal)
    if isdefined(P, :solution)
        S = P.solution
        if meta
            haskey(S.metadata, "tree") && plot_tree(P.SS, P.V, S.metadata["tree"], color="gray", alpha=0.5)
            # TODO: graph (PRM)
        end
        sol && plot_path(P.SS, P.V, S.metadata["path"], color="blue")
        smoothed && haskey(S.metadata, "smoothed_path") && plot_path(S.metadata["smoothed_path"], color="orange")
    end
end