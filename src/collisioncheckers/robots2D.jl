export PointRobot2D

# ---------- Point Robot ----------

type PointRobot2D <: SweptCollisionChecker
    obstacles::Shape2D      # pretty much always Compound2D
end

is_free_state(v::Vector2, CC::PointRobot2D) = !colliding(v, CC.obstacles)
is_free_motion(v::Vector2, w::Vector2, CC::PointRobot2D) = !colliding(Line(v,w), CC.obstacles)
function is_free_path(path::Path, CC::PointRobot2D)
    for i in 1:length(path)-1
        !is_free_motion(path[i], path[i+1], CC) && return false
    end
    true
end
plot(CC::PointRobot2D; kwargs...) = plot(CC.obstacles, color = "red"; kwargs...)
composable(CC::PointRobot2D, opts = (stroke("black"),)) = composable(CC.obstacles, opts)