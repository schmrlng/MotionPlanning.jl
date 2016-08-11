export Goal, RectangleGoal, BallGoal, ConvexHullWorkspaceGoal, ConvexHullStateSpaceGoal
export PointGoal, StateGoal, is_goal_pt, sample_goal

abstract Goal
abstract WorkspaceGoal <: Goal
abstract StateSpaceGoal <: Goal

immutable RectangleGoal{N,T<:AbstractFloat} <: WorkspaceGoal
    lo::SVector{N,T}
    hi::SVector{N,T}
end
RectangleGoal(lo::AbstractVector, hi::AbstractVector) = RectangleGoal(SVector(lo), SVector(hi))
changeprecision{T<:AbstractFloat}(::Type{T}, G::RectangleGoal) = RectangleGoal(changeprecision(T, G.lo),
                                                                               changeprecision(T, G.hi))


immutable BallGoal{N,T<:AbstractFloat} <: WorkspaceGoal
    center::SVector{N,T}
    radius::T
end
BallGoal(center::AbstractVector, radius) = BallGoal(SVector(center), radius)
changeprecision{T<:AbstractFloat}(::Type{T}, G::BallGoal) = BallGoal(changeprecision(T, G.center), T(radius))

immutable ConvexHullWorkspaceGoal{N,T<:AbstractFloat} <: WorkspaceGoal
    pts::Vector{SVector{N,T}}
    A::Matrix{Float64}    # for use with the SCS LP solver
    b::Vector{Float64}    # for use with the SCS LP solver
    c::Vector{Float64}    # for use with the SCS LP solver
    m::Int
    n::Int

    function ConvexHullWorkspaceGoal(pts::Vector{SVector{N,T}})
        m, n = N, length(pts)
        A = n > 1 ? convert(Matrix{Float64}, [statevec2mat(pts); ones(1,n); -eye(n)]) : Array(Float64,0,0)
        b = zeros(m + n + 1); b[m+1] = 1
        c = zeros(n)
        new(pts, A, b, c, m, n)
    end
end
ConvexHullWorkspaceGoal{N,T<:AbstractFloat}(pts::Vector{SVector{N,T}}) = ConvexHullWorkspaceGoal{N,T}(pts)
ConvexHullWorkspaceGoal{V<:AbstractVector}(pts::Vector{V}) = ConvexHullWorkspaceGoal([SVector(p) for p in pts])
PointGoal(pt::AbstractVector) = ConvexHullWorkspaceGoal([pt])
changeprecision{T<:AbstractFloat}(::Type{T}, G::ConvexHullWorkspaceGoal) =
    ConvexHullWorkspaceGoal(changeprecision(T, G.pts))

immutable ConvexHullStateSpaceGoal{S<:State} <: StateSpaceGoal
    sts::Vector{S}
    A::Matrix{Float64}    # for use with the SCS LP solver
    b::Vector{Float64}    # for use with the SCS LP solver
    c::Vector{Float64}    # for use with the SCS LP solver
    m::Int
    n::Int

    function ConvexHullStateSpaceGoal(sts::Vector{S})
        m, n = length(sts[1]), length(sts)
        A = n > 1 ? convert(Matrix{Float64}, [statevec2mat(sts); ones(1,n); -eye(n)]) : Array(Float64,0,0)
        b = zeros(m + n + 1); b[m+1] = 1
        c = zeros(n)
        new(sts, A, b, c, m, n)
    end
end
ConvexHullStateSpaceGoal{S<:State}(sts::Vector{S}) = ConvexHullStateSpaceGoal{S}(sts)
ConvexHullStateSpaceGoal{S<:SVector}(sts::Vector{S}) = ConvexHullStateSpaceGoal{S}(sts)
ConvexHullStateSpaceGoal{S<:FieldVector}(sts::Vector{S}) = ConvexHullStateSpaceGoal{S}(sts)
ConvexHullStateSpaceGoal{V<:AbstractVector}(sts::Vector{V}) = ConvexHullStateSpaceGoal([SVector(s) for s in sts])
StateGoal(s::AbstractVector) = ConvexHullStateSpaceGoal([s])
changeprecision{T<:AbstractFloat}(::Type{T}, G::ConvexHullStateSpaceGoal) =
    ConvexHullStateSpaceGoal(changeprecision(T, G.sts))

plot(G::RectangleGoal, SS::StateSpace) = plot_rectangle(G.lo, G.hi, color = "green")
plot(G::BallGoal, SS::StateSpace) = plot_circle(G.center, G.radius, color = "green")
function plot(G::ConvexHullWorkspaceGoal, SS::StateSpace)
    if G.n == 1    # PointGoal
        plt.scatter(G.pts[1][1], G.pts[1][2], color = "green", zorder=5)
    elseif G.n == 2
        plot_line_segments([G.pts[1]], [G.pts[2]], color = "green", zorder=5, linewidth=2.5)
    else
        plot_polygon(G.pts, color = "green", zorder=5)
    end
end
function plot(G::ConvexHullStateSpaceGoal, SS::StateSpace)
    pts = [state2workspace(s, SS) for s in G.sts]
    if G.n == 1    # StateGoal
        plt.scatter(pts[1][1], pts[1][2], color = "green", zorder=5)
    elseif G.n == 2
        plot_line_segments([pts[1]], [pts[2]], color = "green", zorder=5, linewidth=2.5)
    else
        plot_polygon(pts, color = "green", zorder=5)
    end
end

sample_goal(G::WorkspaceGoal, SS::StateSpace) = workspace2state(sample_goal(G), SS)
sample_goal(G::StateSpaceGoal, SS::StateSpace) = sample_goal(G)

## Rectangle Goal
is_goal_pt(v::State, G::RectangleGoal, SS::StateSpace) = (reduce(&, G.lo .<= state2workspace(v, SS) .<= G.hi) != 0)
sample_goal(G::RectangleGoal) = G.lo + (G.hi - G.lo).*rand(typeof(G.lo))

## Ball Goal
is_goal_pt(v::State, G::BallGoal, SS::StateSpace) = (norm(state2workspace(v, SS) - G.center) <= G.radius)
function sample_goal(G::BallGoal)
    while true
        v = G.center + 2*G.radius*(rand(typeof(G.center)) - .5)
        if norm(v - G.center) <= G.radius
            return v
        end
    end
end

## Convex Hull Goals
function is_goal_pt(v::State, G::ConvexHullWorkspaceGoal, SS::StateSpace)
    pt = state2workspace(v, SS)
    if G.n == 1
        pt == G.pts[1]
    elseif G.n == 2
        isapprox(norm(G.pts[1] - G.pts[2]), norm(G.pts[1] - pt) + norm(G.pts[2] - pt))    # not the most efficient way
    else
        G.b[1:G.m] = pt
        sol = SCS_solve(G.m+G.n+1, G.n, G.A, G.b, G.c, G.m+1, G.n, Int[], 0, Int[], 0, 0, 0, Float64[], 0, verbose=0)
        sol.status != :Infeasible
    end
end
function sample_goal{N,T}(G::ConvexHullWorkspaceGoal{N,T})
    G.n == 1 ? G.pts[1] : G.n == 2 ? G.pts[1] + rand(T)*(G.pts[2] - G.pts[1]) : sum(G.pts .* rand(Dirichlet(G.n, T(1))))
end

function is_goal_pt(v::State, G::ConvexHullStateSpaceGoal, SS::StateSpace)
    if G.n == 1
        v == G.sts[1]
    elseif G.n == 2
        isapprox(norm(G.sts[1] - G.sts[2]), norm(G.sts[1] - v) + norm(G.sts[2] - v))    # not the most efficient way
    else
        G.b[1:G.m] = v
        sol = SCS_solve(G.m+G.n+1, G.n, G.A, G.b, G.c, G.m+1, G.n, Int[], 0, Int[], 0, 0, 0, Float64[], 0, verbose=0)
        sol.status != :Infeasible
    end
end
function sample_goal{S}(G::ConvexHullStateSpaceGoal{S})
    T = eltype(S)
    G.n == 1 ? G.sts[1] : G.n == 2 ? G.sts[1] + rand(T)*(G.sts[2] - G.sts[1]) : sum(G.sts .* rand(Dirichlet(G.n, T(1))))
end