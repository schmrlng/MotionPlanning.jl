export PointRobot2D

# ---------- Point Robot ----------

type PointRobot2D{S<:Shape2D} <: SweptCollisionChecker
    obstacles::S      # pretty much always Compound2D
    count::Int
end
PointRobot2D(obstacles) = PointRobot2D(obstacles, 0)

@unfix is_free_state(v::Vec{2}, CC::PointRobot2D) = !colliding(v, CC.obstacles)
@unfix is_free_motion(v::Vec{2}, w::Vec{2}, CC::PointRobot2D) = (CC.count += 1; !colliding_ends_free(Line(v,w), CC.obstacles))
function is_free_path(path::Path, CC::PointRobot2D)
    for i in 1:length(path)-1
        !is_free_motion(path[i], path[i+1], CC) && return false
    end
    true
end
inflate(CC::PointRobot2D, eps; roundcorners = true) = eps > 0 ? PointRobot2D(inflate(CC.obstacles, eps, roundcorners = roundcorners)) : CC
addobstacle(CC::PointRobot2D, o::Shape2D) = PointRobot2D(Compound2D(CC.obstacles, o))
@unfix addblocker(CC::PointRobot2D, p::Vec{2}, r) = addobstacle(CC, Circle(p, r))
@unfix closest(p::Vec{2}, CC::PointRobot2D, W::Mat{2,2}) = closest(p, CC.obstacles, W)
@unfix closeR(p::Vec{2}, CC::PointRobot2D, W::Mat{2,2}, r2) = closeR(p, CC.obstacles, W, r2)

plot(CC::PointRobot2D, lo = zeros(2), hi = ones(2); kwargs...) = plot(CC.obstacles, color = "red", edgecolor = "none",
                                                                      xmin = lo[1], xmax = hi[1], ymin = lo[2], ymax = hi[2]; kwargs...)
composable(CC::PointRobot2D, opts = (stroke("black"),)) = composable(CC.obstacles, opts)