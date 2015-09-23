export FinalTime, LinearQuadratic, LQOptSteering
export LinearQuadraticQuasiMetricSpace, DoubleIntegrator

### Linear Quadratic Typedefs
## Quasimetrics
abstract LinearQuadratic{T<:FloatingPoint} <: QuasiMetric
type FinalTime{T<:FloatingPoint} <: ControlInfo
    t::T
end
controltype(d::LinearQuadratic) = FinalTime
# Optimal Steering
include("linearquadraticBVP.jl")
type LQOptSteering{T} <: LinearQuadratic{T}
    BVP::LinearQuadratic2BVP{T}
    cmax::T    # for potential pruning
end
function LQOptSteering(A::Matrix, B::Matrix, c::Vector, R::Matrix, cmax = 1.)
    LQOptSteering(LinearQuadratic2BVP(A, B, c, R), cmax)
end
setup_steering(dist::LQOptSteering, r) = (dist.cmax = r)

evaluate(M::LQOptSteering, v::AbstractVector, w::AbstractVector) = steer(M.BVP, v, w, M.cmax)[1]
defaultNN(M::LQOptSteering, init::AbstractVector) = QuasiMetricNN_BruteForce(typeof(init)[], M, init)
function pairwise_distances{S<:State,T<:FloatingPoint}(dist::LQOptSteering{T}, V::Vector{S}, W::Vector{S}, batchsize = 1001)
    M = length(V)
    N = length(W)
    Vmat = hcat(V...)
    Wmat = hcat(W...)
    t = dist.cmax
    Vbarmat = dist.BVP.expAt(t)*Vmat .+ dist.BVP.cdrift(t)
    Ginv = dist.BVP.Ginv(t)
    BRB = dist.BVP.B*inv(dist.BVP.R)*dist.BVP.B'

    cd = pairwise(SqMahalanobis(Ginv*BRB*Ginv), Vbarmat, Wmat)
    LHT = Ginv*(dist.BVP.A*Wmat .+ dist.BVP.c)

    T1 = Distances.dot_percol(Wmat, LHT)
    BLAS.gemm!('T', 'N', -2., Vbarmat, LHT, 1., cd)
    broadcast!(.+, cd, cd, 2*T1')
    broadcast!(.-, cd, 1, cd)

    IS, JS = findn(cd .> 0)
    VTS = [steer(dist.BVP, V[i], W[j], t) for (i,j) in zip(IS, JS)]
    II = find(T[v for (v,t) in VTS] .<= t)
    DS = sparse(IS[II], JS[II], T[v for (v,t) in VTS[II]], M, N)
    if V == W
        for i in 1:N
            DS[i,i] = 0.
        end
    end
    US = sparse(IS[II], JS[II], FinalTime{T}[FinalTime(t) for (v,t) in VTS[II]], M, N)
    DS, US
end
pairwise_distances{S<:State,T<:FloatingPoint}(dist::LQOptSteering{T}, V::Vector{S}, batchsize = 1001) = pairwise_distances(dist, V, V, batchsize)

## Quasimetric Space Instantiation
LinearQuadraticQuasiMetricSpace(d::Int, lo::Vector, hi::Vector, A::Matrix, B::Matrix, c::Vector, R::Matrix, C::Matrix) =
    RealVectorStateSpace(d, lo, hi, LQOptSteering(A, B, c, R), OutputMatrix(C))
function DoubleIntegrator(d::Int, lo = zeros(d), hi = ones(d); vmax = 1.5, r = 1.)
    A = [zeros(d,d) eye(d); zeros(d,2d)]
    B = [zeros(d,d); eye(d)]
    c = zeros(2d)
    R = r*eye(d)
    C = [eye(d) zeros(d,d)]
    LinearQuadraticQuasiMetricSpace(2d, vcat(lo, -vmax*ones(d)), vcat(hi, vmax*ones(d)), A, B, c, R, C)
end

### Waypoint Interpolation
waypoints{T}(v::AbstractVector{T}, w::AbstractVector{T}, M::LQOptSteering{T}, res = 20, t = steer(M.BVP, v, w, M.cmax)[2]) =
    [M.BVP.x(v, w, t, s) for s in linspace(0, t, res)]

collision_waypoints{T}(v::AbstractVector{T}, w::AbstractVector{T}, M::LQOptSteering{T}, s2w::State2Workspace, res = 5, t = steer(M.BVP, v, w, M.cmax)[2]) =
    [state2workspace(M.BVP.x(v, w, t, s), s2w) for s in linspace(0, t, res)]

time_waypoint{T}(v::AbstractVector{T}, w::AbstractVector{T}, M::LQOptSteering{T}, s::T, t = steer(M.BVP, v, w, M.cmax)[2]) =
    s < t ? M.BVP.x(v, w, t, s) : w