export Goal, StateSpaceGoal, GoalState, ConfigSpaceGoal, GoalConfig

abstract type Goal end
# StateSpaceGoal
struct StateSpaceGoal{S} <: Goal
    set::S
end
const GoalState{S<:State} = StateSpaceGoal{SVector{1,S}}
GoalState(x::State) = StateSpaceGoal(SVector((x,)))
Base.in(x::State, G::StateSpaceGoal) = x in G.set
Base.rand(G::StateSpaceGoal) = rand(G.set)
@recipe function f(G::StateSpaceGoal; dims=(1, 2))
    if G.set isa AbstractArray{<:State}
        seriestype  --> :scatter
        markercolor :=  :match
        x, y = dims
        [q[x] for q in G.set], [q[y] for q in G.set]
    else
        dims --> dims
        G.set
    end
end

# ConfigSpaceGoal
struct ConfigSpaceGoal{C,S2C,C2S} <: Goal
    projection::C
    state2config::S2C
    config2state::C2S
end
const GoalConfig{C<:Config,S2C,C2S} = ConfigSpaceGoal{SVector{1,C},S2C,C2S}
GoalConfig(q::Config, state2config, config2state) = ConfigSpaceGoal(SVector((q,)), state2config, config2state)
Base.in(x::State, G::ConfigSpaceGoal) = G.state2config(x) in G.projection
Base.rand(G::ConfigSpaceGoal) = G.config2state(rand(G.projection))
@recipe function f(G::ConfigSpaceGoal; dims=(1, 2))
    dims --> dims
    G.projection
end
