export Goal, RectangleGoal, BallGoal, PointGoal, StateGoal, is_goal_pt, sample_goal

abstract Goal
abstract WorkspaceGoal <: Goal
abstract StateSpaceGoal <: Goal

immutable RectangleGoal{T<:AbstractFloat} <: WorkspaceGoal
    bounds::Matrix{T}
end

immutable BallGoal{T<:AbstractFloat} <: WorkspaceGoal
    center::AbstractVector{T}
    radius::T
end

immutable PointGoal{T<:AbstractFloat} <: WorkspaceGoal
    pt::AbstractVector{T}
end

immutable StateGoal{S<:State} <: StateSpaceGoal
    st::S
end

plot(G::RectangleGoal, SS::StateSpace) = plot_rectangle(G.bounds, color = "green")
plot(G::BallGoal, SS::StateSpace) = plot_circle(G.center, G.radius, color = "green")
plot(G::PointGoal, SS::StateSpace) = scatter(G.pt[1:2]..., color = "green", zorder=5)
plot(G::StateGoal, SS::StateSpace) = scatter(state2workspace(G.st, SS)[1:2]..., color = "green", zorder=5)

sample_goal(G::WorkspaceGoal, SS::StateSpace) = workspace2state(sample_goal(G), SS)
sample_goal(G::StateSpaceGoal, SS::StateSpace) = sample_goal(G)

## Rectangle Goal
is_goal_pt(v::State, G::RectangleGoal, SS::StateSpace) = all(G.bounds[:,1] .<= state2workspace(v, SS) .<= G.bounds[:,2])
sample_goal(G::RectangleGoal) = (G.bounds[:,1] + (G.bounds[:,2] - G.bounds[:,1]).*rand(size(G.bounds,1)))

## Ball Goal
is_goal_pt(v::State, G::BallGoal, SS::StateSpace) = (norm(convert(typeof(G.center), state2workspace(v, SS)) - G.center) <= G.radius)
function sample_goal(G::BallGoal)
    while true
        v = G.center + 2*G.radius*(convert(typeof(G.center), rand(length(G.center))) - .5)
        if is_goal_pt(v, G)
            return v
        end
    end
end

## Point Goal
is_goal_pt(v::State, G::PointGoal, SS::StateSpace) = (state2workspace(v, SS) == G.pt)
sample_goal(G::PointGoal) = G.pt

## State Goal
is_goal_pt(v::State, G::StateGoal, SS::StateSpace) = (v == G.st)
sample_goal(G::StateGoal) = G.st