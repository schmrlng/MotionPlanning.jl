export PointRobot2D

# ---------- Point Robot ----------

type PointRobot2D{S<:Shape2D} <: SweptCollisionChecker
    obstacles::S      # pretty much always Compound2D
    count::Int
end
PointRobot2D(obstacles) = PointRobot2D(obstacles, 0)
changeprecision{T<:AbstractFloat}(::Type{T}, CC::PointRobot2D) = PointRobot2D(changeprecision(T, CC.obstacles))

is_free_state(v::AbstractVector, CC::PointRobot2D) = !colliding(v, CC.obstacles)
is_free_motion(v::AbstractVector, w::AbstractVector, CC::PointRobot2D) = (CC.count += 1;
                                                                          !colliding(Line(v,w), CC.obstacles))
function is_free_path(path::Path, CC::PointRobot2D)
    for i in 1:length(path)-1
        !is_free_motion(path[i], path[i+1], CC) && return false
    end
    true
end
inflate(CC::PointRobot2D, ε; roundcorners = true) =
    ε > 0 ? PointRobot2D(inflate(CC.obstacles, ε, roundcorners = roundcorners)) : CC
addobstacle(CC::PointRobot2D, o::Shape2D) = PointRobot2D(Compound2D(CC.obstacles, o))
addblocker(CC::PointRobot2D, p::AbstractVector, r) = addobstacle(CC, Circle(p, r))
closest(p::AbstractVector, CC::PointRobot2D, W::AbstractMatrix) = closest(p, CC.obstacles, W)
closeR(p::AbstractVector, CC::PointRobot2D, W::AbstractMatrix, r2) = closeR(p, CC.obstacles, W, r2)

function plot(CC::PointRobot2D, lo = zeros(2), hi = ones(2); kwargs...)
    plot(CC.obstacles, color = "red", edgecolor = "none",
         xmin = lo[1], xmax = hi[1], ymin = lo[2], ymax = hi[2]; kwargs...)
end
# composable(CC::PointRobot2D, opts = (stroke("black"),)) = composable(CC.obstacles, opts)