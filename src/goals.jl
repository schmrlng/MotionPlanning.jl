export Goal, RectangleGoal, BallGoal, PointGoal, is_goal_pt, sample_goal

abstract Goal

immutable RectangleGoal{T<:FloatingPoint} <: Goal
    bounds::Matrix{T}
end

immutable BallGoal{T<:FloatingPoint} <: Goal
    center::AbstractVector{T}
    radius::T
end

immutable PointGoal{T<:FloatingPoint} <: Goal
    pt::AbstractVector{T}
end

plot(G::RectangleGoal) = plot_rectangle(G.bounds, color = "green")
plot(G::BallGoal) = plot_circle(G.center, G.radius, xmax = 1, ymax = 1, color = "green")
plot(G::PointGoal) = scatter(G.pt[1], G.pt[2], color = "green", zorder=5)

## Rectangle Goal

is_goal_pt(v::AbstractVector, G::RectangleGoal) = all(G.bounds[:,1] .<= v .<= G.bounds[:,2])
sample_goal(G::RectangleGoal) = (G.bounds[:,1] + (G.bounds[:,2] - G.bounds[:,1]).*rand(size(G.bounds,1)))

## Ball Goal

is_goal_pt(v::AbstractVector, G::BallGoal) = (norm(v - G.center) <= G.radius)
function sample_goal(G::BallGoal)
    while true
        v = G.center + 2*G.radius*(convert(typeof(G.center), rand(length(G.center))) - .5)
        if is_goal_pt(v, G)
            return v
        end
    end
end

## Point Goal

is_goal_pt(v::AbstractVector, G::PointGoal) = (v == G.pt)
sample_goal(G::PointGoal) = G.pt