export Goal, RectangleGoal, BallGoal, PointGoal, is_goal_pt, sample_goal

abstract Goal

immutable RectangleGoal{T<:FloatingPoint} <: Goal
    bounds::Matrix{T}
end

immutable BallGoal{T<:FloatingPoint} <: Goal
    center::AbstractVector{T}
    radius::T
end

immutable PointGoal <: Goal
    pt::State
end

plot(G::RectangleGoal, SS::StateSpace) = plot_rectangle(G.bounds, color = "green")
plot(G::BallGoal, SS::StateSpace) = plot_circle(G.center, G.radius, color = "green")
plot(G::PointGoal, SS::StateSpace) = scatter(state2workspace(G.pt, SS)[1:2]..., color = "green", zorder=5)

## Rectangle Goal
is_goal_pt(v::State, G::RectangleGoal, SS::StateSpace) = all(G.bounds[:,1] .<= state2workspace(v, SS) .<= G.bounds[:,2])
sample_goal(G::RectangleGoal) = (G.bounds[:,1] + (G.bounds[:,2] - G.bounds[:,1]).*rand(size(G.bounds,1)))  # TODO: workspace2state

## Ball Goal
is_goal_pt(v::State, G::BallGoal, SS::StateSpace) = (norm(convert(typeof(G.center), state2workspace(v, SS)) - G.center) <= G.radius)
function sample_goal(G::BallGoal)    # TODO: workspace2state
    while true
        v = G.center + 2*G.radius*(convert(typeof(G.center), rand(length(G.center))) - .5)
        if is_goal_pt(v, G)
            return v
        end
    end
end

## Point Goal
is_goal_pt(v::State, G::PointGoal) = (v == G.pt)
sample_goal(G::PointGoal, SS::StateSpace) = G.pt