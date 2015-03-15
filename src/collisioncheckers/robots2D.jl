export PointRobot2D

# ---------- Point Robot ----------

type PointRobot2D <: SweptCollisionChecker
    obstacles::Shape2D      # pretty much always Compound2D
    count::Int
end
PointRobot2D(obstacles) = PointRobot2D(obstacles, 0)

is_free_state(v::Vector2, CC::PointRobot2D) = !colliding(v, CC.obstacles)
is_free_motion(v::Vector2, w::Vector2, CC::PointRobot2D) = (CC.count += 1; !colliding(Line(v,w), CC.obstacles))
function is_free_path(path::Path, CC::PointRobot2D)
    for i in 1:length(path)-1
        !is_free_motion(path[i], path[i+1], CC) && return false
    end
    true
end
inflate(CC::PointRobot2D, eps) = eps > 0 ? PointRobot2D(inflate(CC.obstacles, eps)) : CC
plot(CC::PointRobot2D, lo = zeros(2), hi = zeros(2); kwargs...) = plot(CC.obstacles, color = "red", edgecolor = "none",
                                                                       xmin = lo[1], xmax = hi[1], ymin = lo[1], ymax = hi[1]; kwargs...)
composable(CC::PointRobot2D, opts = (stroke("black"),)) = composable(CC.obstacles, opts)