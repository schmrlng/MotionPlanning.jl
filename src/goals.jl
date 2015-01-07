abstract Goal

immutable RectangleGoal{T<:FloatingPoint} <: Goal
    bounds::Matrix{T}
end

immutable BallGoal{T<:FloatingPoint} <: Goal
    center::Vector{T}
    radius::T
end

immutable PointGoal{T<:FloatingPoint} <: Goal
    pt::Vector{T}
end

### RectangleGoal

plot_goal(G::RectangleGoal) = plot_rectangle(G.bounds, color = "green")
plot_goal(G::BallGoal) = plot_circle(G.center, G.radius, xmax = 1, ymax = 1, color = "green")
plot_goal(G::PointGoal) = scatter(G.pt[1], G.pt[2], color = "green", zorder=5)

## Rectangle Goal

is_goal_pt(v::Vector, G::RectangleGoal) = all(G.bounds[:,1] .<= v .<= G.bounds[:,2])
sample_goal(G::RectangleGoal) = (G.bounds[:,1] + (G.bounds[:,2] - G.bounds[:,1]).*rand(size(G.bounds,1)))

## Ball Goal

is_goal_pt(v::Vector, G::BallGoal) = (norm(v - G.center) <= G.radius)
function sample_goal(G::BallGoal)
    while true
        v = G.center + 2*G.radius*(rand(length(G.center)) - .5)
        if is_goal_pt(v, G)
            return v
        end
    end
end

## Point Goal

is_goal_pt(v::Vector, G::PointGoal) = (v == G.pt)
sample_goal(G::PointGoal) = G.pt