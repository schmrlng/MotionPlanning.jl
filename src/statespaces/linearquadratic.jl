export FinalTime, LinearQuadratic, LinearQuadraticOpt
export LinearQuadraticStateSpace, DoubleIntegrator
export waypoints

### Linear Quadratic Steering
type FinalTime{T<:FloatingPoint} <: ControlInfo
    t::T
end
abstract LinearQuadratic{T<:FloatingPoint} <: QuasiMetric
## Optimal Steering
type LinearQuadraticOpt{T} <: LinearQuadratic{T}
    steer::Function
    cmax::T         # for potential pruning
end
LinearQuadraticOpt(cost::Function) = LinearQuadraticOpt(cost::Function, 1.)
function LinearQuadraticOpt{T<:FloatingPoint}(cost::Function, cost_deriv::Function, cmax::T = 1.)
    function topt(x0, x1, tm)
        # Bisection
        b = tm
        cost_deriv(x0, x1, b) < 0 && return tm
        a = .01
        while cost_deriv(x0, x1, a) > 0; a /= 2.; end
        m = 0.
        for i in 1:10
            m = (a+b)/2
            cost_deriv(x0, x1, m) > 0 ? b = m : a = m
        end
        m
    end
    function steer(x0, x1, cm)
        t = topt(x0, x1, cm)
        cost(x0, x1, t), FinalTime(t)
    end
    LinearQuadraticOpt(steer, cmax)
end
## TODO: Approximate(ly Optimal) Steering
# type LinearQuadraticApprox{T<:FloatingPoint} <: LinearQuadratic
#
# end
controltype(d::LinearQuadratic) = FinalTime

### Linear Quadratic State Space
immutable LinearQuadraticStateSpace{T<:FloatingPoint} <: DifferentialStateSpace
    dim::Int
    lo::Vector{T}
    hi::Vector{T}
    dist::LinearQuadratic{T}
    x::Function

    A::Matrix{T}
    B::Matrix{T}
    c::Vector{T}
    R::Matrix{T}
end
## Optimal Steering
LinearQuadraticStateSpace(dim::Int, lo::Vector, hi::Vector,
                          A::Matrix, B::Matrix, c::Vector, R::Matrix,
                          cost::Function, cost_deriv::Function, x::Function) =
    LinearQuadraticStateSpace(dim, lo, hi, LinearQuadraticOpt(cost, cost_deriv), x, A, B, c, R)
## TODO: Approximate Steering
# LinearQuadraticStateSpace(dim::Int, lo::Vector{T}, hi::Vector{T},
#                           A::Matrix{T}, B::Matrix{T}, c::Vector{T}, R::Matrix{T},
#                           G::Function, Ginv::Function, expAt::Function, cdrift::Function) =
#     LinearQuadraticStateSpace{T}(dim, lo, hi, A, B, c, R, G, Ginv, expAt, cdrift) # so constructor works without {}

vector_to_state{T}(v::AbstractVector{T}, SS::LinearQuadraticStateSpace) = v
sample_space(SS::LinearQuadraticStateSpace) = vector_to_state(SS.lo + rand(SS.dim).*(SS.hi-SS.lo), SS)   # TODO: @devec
function volume(SS::LinearQuadraticStateSpace)
    warn("TODO: what is volume for a LinearQuadraticStateSpace?")
    prod(SS.hi-SS.lo)
end
function defaultNN(SS::LinearQuadraticStateSpace, init)
    V = typeof(init)[init]
    QuasiMetricNN_BruteForce(V, SS.dist)
end

function pairwise_distances{S<:State,T<:FloatingPoint}(dist::LinearQuadraticOpt{T}, V::Vector{S})
    N = length(V)
    VM = hcat(V...)
    DS = Array(T, N, N)
    US = Array(FinalTime, N, N)
    for j = 1 : N
        vj = view(VM,:,j)
        for i = 1 : N
            @inbounds DS[i,j], US[i,j] = dist.steer(view(VM,:,i), vj, dist.cmax)
        end
    end
    DS, US
end

waypoints(i, j, NN::QuasiMetricNN, SS::LinearQuadraticStateSpace, res=5) = [Vector2(SS.x(NN[i], NN[j], NN.US[i,j].t, s)) for s in linspace(0, NN.US[i,j].t, res)]
function inbounds(v, SS::LinearQuadraticStateSpace)
    for i in 1:length(v)
        (SS.lo[i] > v[i] || v[i] > SS.hi[i]) && return false
    end
    true
end
is_free_state(v, CC::PointRobot2D, SS::LinearQuadraticStateSpace) = inbounds(v, SS) && is_free_state(Vector2(v), CC)
function is_free_motion(v, w, CC::PointRobot2D, SS::LinearQuadraticStateSpace)   # TODO: inputs V, i, j instead of v, w
    t = (SS.dist.steer(v, w, SS.dist.cmax)[2]).t   # terrible, hmm
    for s in linspace(0, t, 5)
        y = SS.x(v, w, t, s)
        !inbounds(y, SS) && return false
        vy = Vector2(y)
        s > 0 && !is_free_motion(vx, vy, CC) && return false
        vx = vy
    end
    true
end
# TODO: is_free_path(path, CC::PointRobot2D, SS::LinearQuadraticStateSpace)

function plot_tree(SS::LinearQuadraticStateSpace, NN::QuasiMetricNN, A; kwargs...)
    pts = hcat(NN[find(A)]...)
    scatter(pts[1,:], pts[2,:], zorder=1; kwargs...)
    X = vcat([[hcat(waypoints(A[v], v, NN, SS, 20)...)[1,:]', nothing] for v in find(A)]...)
    Y = vcat([[hcat(waypoints(A[v], v, NN, SS, 20)...)[2,:]', nothing] for v in find(A)]...)
    plt.plot(X, Y, linewidth=.5, linestyle="-", zorder=1; kwargs...)
end

function plot_path(SS::LinearQuadraticStateSpace, NN::QuasiMetricNN, sol; kwargs...)
    wps = hcat([hcat(waypoints(sol[i], sol[i+1], NN, SS, 20)...) for i in 1:length(sol)-1]...)
    length(sol) > 1 && plot_path(wps; kwargs...)
    # plt.quiver([wps[row,1:3:end]' for row in 1:4]..., zorder=5, width=.003, headwidth=8)
end

### Double Integrator State Space

function DoubleIntegrator(d::Int, lo = zeros(d), hi = ones(d); vmax = 1.5, r = 1.)
    A = [zeros(d,d) eye(d); zeros(d,2d)]
    B = [zeros(d,d); eye(d)]
    c = zeros(2d)
    R = r*eye(d)
    G(t) = [(t^3/(3r))*eye(d) (t^2/(2r))*eye(d);
            (t^2/(2r))*eye(d) (t/r)*eye(d)]
    Ginv(t) = [(12r/t^3)*eye(d) (-6r/t^2)*eye(d);
               (-6r/t^2)*eye(d) (4r/t)*eye(d)]
    expAt(t) = [eye(d) t*eye(d); zeros(d,d) eye(d)]
    cdrift(t) = zeros(2d)

    if d == 2       # TODO: figure how to do this symbolically using Calculus.jl
        dx(x0, x1, t) = (x1[1] - x0[1] - t*x0[3], x1[2] - x0[2] - t*x0[4], x1[3] - x0[3], x1[4] - x0[4])
        function cost(x0, x1, t)
            dx1, dx2, dx3, dx4 = dx(x0, x1, t)
            t + 12r/t^3*(dx1^2 + dx2^2) - 12r/t^2*(dx1*dx3 + dx2*dx4) + 4r/t*(dx3^2 + dx4^2)
        end
        function cost_deriv(x0, x1, t)
            dx1, dx2, dx3, dx4 = dx(x0, x1, t)
            1. - 36r/t^4*(dx1^2 + dx2^2) - 12r/t^3*(2dx1*x0[3] + 2dx2*x0[4]) + 24r/t^3*(dx1*dx3 + dx2*dx4) + 12r/t^2*(x0[3]*dx3 + x0[4]*dx4) - 4r/t^2*(dx3^2 + dx4^2)
        end
        function x(x0, x1, t, s)
            dx1, dx2, dx3, dx4 = dx(x0, x1, t)
            m11, m12, m21, m22 = (-2s^3/t^3 + 3s^2/t^2, s^3/t^2 - s^2/t, -6s^2/t^3 + 6s/t^2, 3s^2/t^2 - 2s/t)
            [x0[1] + s*x0[3] + m11*dx1 + m12*dx3,
             x0[2] + s*x0[4] + m11*dx2 + m12*dx4,
             x0[3] + m21*dx1 + m22*dx3,
             x0[4] + m21*dx2 + m22*dx4]
        end
        return LinearQuadraticStateSpace(2d, [lo, -vmax*ones(d)], [hi, vmax*ones(d)], A, B, c, R, cost, cost_deriv, x)
    else
        error("TODO: Implement LinearQuadraticApprox before higher dimensions will work")
    end
end


# TODO: old code relevant to approx opt implementation
# function pairwise_distances_approx_opt{T}(V::Vector{Vector{T}}, SS::LinearQuadraticStateSpace, t_bound::Float64, res::Int64 = 10)
#     N = length(V)
#     V1 = hcat(V...)
#     t = t_bound/2.
#     V0bar = SS.expAt(t)*V1 .+ SS.cdrift(t)
#     return (t + pairwise(SqMahalanobis(SS.Ginv(t)), V0bar, V1)), fill(t, N, N)
# end