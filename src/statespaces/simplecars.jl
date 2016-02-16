export ReedsSheppMetricSpace, DubinsQuasiMetricSpace
export ReedsSheppExact, DubinsExact, SimpleCarMetric

### (Quasi)Metric Typedefs/Evaluation
immutable CarSegment{T<:AbstractFloat}
    t::Int          # segment type (L = 1, S = 0, R = -1)
    d::T            # segment length (radians for curved segments)
end
typealias CarPath{T} Vector{CarSegment{T}}
immutable ReedsSheppExact{T<:AbstractFloat} <: Metric
    r::T
    p::CarPath{T}    # scratch space to avoid extra memory allocation

    ReedsSheppExact(r::T) = new(r, Array(CarSegment{T}, 5))
end
ReedsSheppExact{T<:AbstractFloat}(r::T) = ReedsSheppExact{T}(r)
immutable DubinsExact{T<:AbstractFloat} <: QuasiMetric
    r::T
    p::CarPath{T}    # scratch space to avoid extra memory allocation

    DubinsExact(r::T) = new(r, Array(CarSegment{T}, 3))
end
DubinsExact{T<:AbstractFloat}(r::T) = DubinsExact{T}(r)
typealias SimpleCarMetric{T} Union{ReedsSheppExact{T}, DubinsExact{T}}

evaluate(RS::ReedsSheppExact, v::SE2State, w::SE2State) = reedsshepp(v, w, RS.r, RS.p)[1]
evaluate(D::DubinsExact, v::SE2State, w::SE2State) = dubins(v, w, D.r, D.p)[1]
changeprecision{T<:AbstractFloat}(::Type{T}, dist::ReedsSheppExact) = ReedsSheppExact(T(dist.r))
changeprecision{T<:AbstractFloat}(::Type{T}, dist::DubinsExact) = DubinsExact(T(dist.r))

## (Quasi)Metric Space Instantiation
ReedsSheppMetricSpace{T}(r::T, lo = zero(Vec{2,T}), hi = one(Vec{2,T})) =
    SE2StateSpace(lo, hi, ChoppedMetric(ReedsSheppExact(r), Euclidean(), T(Inf)), ExtractVector())
DubinsQuasiMetricSpace{T}(r::T, lo = zero(Vec{2,T}), hi = one(Vec{2,T})) =
    SE2StateSpace(lo, hi, ChoppedQuasiMetric(DubinsExact(r), Euclidean(), T(Inf)), ExtractVector())

function helper_data_structures{S<:SE2State,R<:ReedsSheppExact}(V::Vector{S}, M::ChoppedMetric{R})
    TreeDistanceDS(KDTree(statevec2mat(V)[1:2,:], M.lowerbound; reorder=false)), EmptyControlCache()
end
function helper_data_structures{S<:SE2State,R<:DubinsExact}(V::Vector{S}, M::ChoppedQuasiMetric{R})
    DS = TreeDistanceDS(KDTree(statevec2mat(V)[1:2,:], M.lowerbound; reorder=false))
    US = EmptyControlCache()   # TODO: non-empty ControlDS
    DS, US, DS, US
end
# Required for CLBM
evaluate(M::Euclidean, v::SE2State, w::SE2State) = norm(v.x - w.x)
inrange(tree::NearestNeighbors.NNTree, v::SE2State, args...) = inrange(tree, dense(v.x), args...)

### Steering
function propagate{T}(M::SimpleCarMetric{T}, v::SE2State{T}, u::StepControl{T,2})
    s = u.u[1]
    invr = u.u[2]
    if abs(u.t*s*invr) > 10*eps(T)
        SE2State(v.x[1] + (sin(v.θ + u.t*s*invr) - sin(v.θ))/invr,
                 v.x[2] + (cos(v.θ) - cos(v.θ + u.t*s*invr))/invr,
                 mod2pi(v.θ + u.t*s*invr))
    else
        SE2State(v.x[1] + u.t*s*cos(v.θ),
                 v.x[2] + u.t*s*sin(v.θ),
                 mod2pi(v.θ + u.t*s*invr))
    end
end
steering_control(RS::ReedsSheppExact, v::SE2State, w::SE2State) = reedsshepp(v, w, RS.r, RS.p)[2]
steering_control(D::DubinsExact, v::SE2State, w::SE2State) = dubins(v, w, D.r, D.p)[2]
function collision_waypoints{T}(M::SimpleCarMetric{T}, v::SE2State{T}, u::StepControl{T,2})
    s = u.u[1]
    invr = u.u[2]
    θres = T(pi)/12
    m = floor(Int, u.t*s*invr / θres)
    if m == 0
        [v]
    else
        unshift!([SE2State(v.x[1] + (sin(v.θ + i*θres) - sin(v.θ))/invr,
                           v.x[2] + (cos(v.θ) - cos(v.θ + i*θres))/invr,
                           mod2pi(v.θ + i*θres)) for i in 1:m], v)
    end
end

### Simple Car Steering Nuts and Bolts
carsegment2stepcontrol{T}(x::CarSegment{T}, r::T) = StepControl(abs(x.d)*r, Vec(T(sign(x.d)), x.t/r))

## Dubins Steering Nuts and Bolts
function dubinsLSL!{T}(d::T, a::T, b::T, c::T, path::CarPath{T})
    ca, sa, cb, sb = cos(a), sin(a), cos(b), sin(b)
    tmp = 2 + d*d - 2*(ca*cb + sa*sb - d*(sa - sb))
    tmp < 0 && return c
    th = atan2(cb - ca, d + sa - sb)
    t = mod2pi(-a + th)
    p = sqrt(max(tmp, T(0)))
    q = mod2pi(b - th)
    cnew = t + p + q
    c <= cnew && return c
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(0, p)
    path[3] = CarSegment{T}(1, q)
    cnew
end

function dubinsRSR!{T}(d::T, a::T, b::T, c::T, path::CarPath{T})
    ca, sa, cb, sb = cos(a), sin(a), cos(b), sin(b)
    tmp = 2 + d*d - 2*(ca*cb + sa*sb - d*(sb - sa))
    tmp < 0 && return c
    th = atan2(ca - cb, d - sa + sb)
    t = mod2pi(a - th)
    p = sqrt(max(tmp, T(0)))
    q = mod2pi(-b + th)
    cnew = t + p + q
    c <= cnew && return c
    path[1] = CarSegment{T}(-1, t)
    path[2] = CarSegment{T}(0, p)
    path[3] = CarSegment{T}(-1, q)
    cnew
end

function dubinsRSL!{T}(d::T, a::T, b::T, c::T, path::CarPath{T})
    ca, sa, cb, sb = cos(a), sin(a), cos(b), sin(b)
    tmp = d * d - 2 + 2 * (ca*cb + sa*sb - d * (sa + sb))
    tmp < 0 && return c
    p = sqrt(max(tmp, T(0)))
    th = atan2(ca + cb, d - sa - sb) - atan2(T(2), p)
    t = mod2pi(a - th)
    q = mod2pi(b - th)
    cnew = t + p + q
    c <= cnew && return c
    path[1] = CarSegment{T}(-1, t)
    path[2] = CarSegment{T}(0, p)
    path[3] = CarSegment{T}(1, q)
    cnew
end

function dubinsLSR!{T}(d::T, a::T, b::T, c::T, path::CarPath{T})
    ca, sa, cb, sb = cos(a), sin(a), cos(b), sin(b)
    tmp = -2 + d * d + 2 * (ca*cb + sa*sb + d * (sa + sb))
    tmp < 0 && return c
    p = sqrt(max(tmp, T(0)))
    th = atan2(-ca - cb, d + sa + sb) - atan2(-T(2), p)
    t = mod2pi(-a + th)
    q = mod2pi(-b + th)
    cnew = t + p + q
    c <= cnew && return c
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(0, p)
    path[3] = CarSegment{T}(-1, q)
    cnew
end

function dubinsRLR!{T}(d::T, a::T, b::T, c::T, path::CarPath{T})
    ca, sa, cb, sb = cos(a), sin(a), cos(b), sin(b)
    tmp = (6 - d * d  + 2 * (ca*cb + sa*sb + d * (sa - sb))) / 8
    abs(tmp) >= 1 && return c
    p = 2*T(pi) - acos(tmp)
    th = atan2(ca - cb, d - sa + sb)
    t = mod2pi(a - th + p/2)
    q = mod2pi(a - b - t + p)
    cnew = t + p + q
    c <= cnew && return c
    path[1] = CarSegment{T}(-1, t)
    path[2] = CarSegment{T}(1, p)
    path[3] = CarSegment{T}(-1, q)
    cnew
end

function dubinsLRL!{T}(d::T, a::T, b::T, c::T, path::CarPath{T})
    ca, sa, cb, sb = cos(a), sin(a), cos(b), sin(b)
    tmp = (6 - d * d  + 2 * (ca*cb + sa*sb - d * (sa - sb))) / 8
    abs(tmp) >= 1 && return c
    p = 2*T(pi) - acos(tmp)
    th = atan2(-ca + cb, d + sa - sb)
    t = mod2pi(-a + th + p/2)
    q = mod2pi(b - a - t + p)
    cnew = t + p + q
    c <= cnew && return c
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(-1, p)
    path[3] = CarSegment{T}(1, q)
    cnew
end

function dubins{T}(s1::SE2State{T}, s2::SE2State{T}, r::T = T(1), pmin = Array(CarSegment{T}, 3))
    const dubinsTypes = (dubinsLSL!, dubinsRSR!, dubinsRSL!, dubinsLSR!, dubinsRLR!, dubinsLRL!)
    v = (s2.x - s1.x) / r
    d = norm(v)
    th = atan2(v)
    a = mod2pi(s1.θ - th)
    b = mod2pi(s2.θ - th)

    cmin = T(Inf)
    for dubinsType! in dubinsTypes
        cmin = dubinsType!(d, a, b, cmin, pmin)
    end
    cmin*r, [carsegment2stepcontrol(pmin[i], r) for i in 1:3]
end

function save_dubins_cache(Xmax = 5., Ymax = 5., Nx = 101, Ny = 101, Nt = 101,
                           fname = joinpath(Pkg.dir("MotionPlanning"), "src", "statespaces",
                                            @sprintf("Dubins_%.2f_%.2f_%d_%d_%d.gz", Xmax, Ymax, Nx, Ny, Nt)))
    segment_to_string(s) = (s.t == 1 ? " L " : s.t == 0 ? " S " : " R ") * @sprintf("%.3f", s.d)
    v = SE2State(0., 0., 0.)
    fh = GZip.open(fname, "w")
    for t in [2pi*(i-1)/Nt for i in 1:Nt], y in linspace(-Ymax, Ymax, Ny), x in linspace(-Xmax, Xmax, Nx)
        c, p = dubins(v, SE2State(x, y, t))
        cacheline = @sprintf("%.6f ", c) * prod([segment_to_string(s) for s in p]) * "\n"
        write(fh, cacheline)
    end
    GZip.close(fh)
end

## Reeds-Shepp Steering Nuts and Bolts
# Utilities (pedantic about typing to guard against problems)
R{T}(x::T, y::T) = sqrt(x*x + y*y), atan2(y, x)
function M{T}(t::T)
    m = mod2pi(t)
    m > pi ? m - 2*T(pi) : m
end
function Tau{T}(u::T, v::T, E::T, N::T)
    delta = M(u - v)
    A = sin(u) - sin(delta)
    B = cos(u) - cos(delta) - 1
    r, θ = R(E*A + N*B, N*A - E*B)
    t = 2*cos(delta) - 2*cos(v) - 2*cos(u) + 3
    t < 0 ? M(θ + pi) : M(θ)
end
Omega{T}(u::T, v::T, E::T, N::T, t::T) = M(Tau(u, v, E, N) - u + v - t)
timeflip{T}(s::SE2State{T}) = SE2State(-s.x[1], s.x[2], -s.θ)
reflect{T}(s::SE2State{T}) = SE2State(s.x[1], -s.x[2], -s.θ)
backwards{T}(s::SE2State{T}) = SE2State(s.x[1]*cos(s.θ) + s.x[2]*sin(s.θ), s.x[1]*sin(s.θ) - s.x[2]*cos(s.θ), s.θ)
function timeflip!(u::ZeroOrderHoldControl)
    for i in 1:length(u)
        u[i] = StepControl(u[i].t, Vec(-u[i].u[1], u[i].u[2]))
    end
    u
end
function reflect!(u::ZeroOrderHoldControl)
    for i in 1:length(u)
        u[i] = StepControl(u[i].t, Vec(u[i].u[1], -u[i].u[2]))
    end
    u
end
backwards!(u::ZeroOrderHoldControl) = reverse!(u)

# TODO: something about an enum type in Julia v0.4?
const POST, POST_T, POST_R, POST_B, POST_R_T, POST_B_T, POST_B_R, POST_B_R_T = 0, 1, 2, 3, 4, 5, 6, 7

# Gotta check 'em all
function reedsshepp{T}(s1::SE2State{T}, s2::SE2State{T}, r::T = T(1), p = Array(CarSegment{T}, 5))
    dx, dy = (s2.x - s1.x) / r
    ct, st = cos(s1.θ), sin(s1.θ)
    target = SE2State(dx*ct + dy*st, -dx*st + dy*ct, mod2pi(s2.θ - s1.θ))

    tTarget = timeflip(target)
    rTarget = reflect(target)
    trTarget = reflect(tTarget)
    bTarget = backwards(target)
    btTarget = timeflip(bTarget)
    brTarget = reflect(bTarget)
    btrTarget = reflect(btTarget)
    
    c, l = T(Inf), 0
    # Only Tim Wheeler knows what the hell is happening down below
    # (8.1) C S C
    b, c, l = LpSpLp!(target, c, l, p); b && (post = POST)
    b, c, l = LpSpLp!(tTarget, c, l, p); b && (post = POST_T)
    b, c, l = LpSpLp!(rTarget, c, l, p); b && (post = POST_R)
    b, c, l = LpSpLp!(trTarget, c, l, p); b && (post = POST_R_T)

    # (8.2) C S C
    b, c, l = LpSpRp!(target, c, l, p); b && (post = POST)
    b, c, l = LpSpRp!(tTarget, c, l, p); b && (post = POST_T)
    b, c, l = LpSpRp!(rTarget, c, l, p); b && (post = POST_R)
    b, c, l = LpSpRp!(trTarget, c, l, p); b && (post = POST_R_T)

    # (8.3) C|C|C
    b, c, l = LpRmLp!(target, c, l, p); b && (post = POST)
    # b, c, l = LpRmLp!(tTarget, c, l, p); b && (post = POST_T) # (redundant)
    b, c, l = LpRmLp!(rTarget, c, l, p); b && (post = POST_R)
    # b, c, l = LpRmLp!(trTarget, c, l, p); b && (post = POST_R_T) # (redundant)

    # (8.4) C|C C
    b, c, l = LpRmLm!(target, c, l, p); b && (post = POST)
    b, c, l = LpRmLm!(tTarget, c, l, p); b && (post = POST_T)
    b, c, l = LpRmLm!(rTarget, c, l, p); b && (post = POST_R)
    b, c, l = LpRmLm!(trTarget, c, l, p); b && (post = POST_R_T)
    b, c, l = LpRmLm!(bTarget, c, l, p); b && (post = POST_B)
    b, c, l = LpRmLm!(btTarget, c, l, p); b && (post = POST_B_T)
    b, c, l = LpRmLm!(brTarget, c, l, p); b && (post = POST_B_R)
    b, c, l = LpRmLm!(btrTarget, c, l, p); b && (post = POST_B_R_T)

    # (8.7) C Cu|Cu C
    b, c, l = LpRpuLmuRm!(target, c, l, p); b && (post = POST)
    b, c, l = LpRpuLmuRm!(tTarget, c, l, p); b && (post = POST_T)
    b, c, l = LpRpuLmuRm!(rTarget, c, l, p); b && (post = POST_R)
    b, c, l = LpRpuLmuRm!(trTarget, c, l, p); b && (post = POST_R_T)

    # (8.8) C|Cu Cu|C
    b, c, l = LpRmuLmuRp!(target, c, l, p); b && (post = POST)
    b, c, l = LpRmuLmuRp!(tTarget, c, l, p); b && (post = POST_T)
    b, c, l = LpRmuLmuRp!(rTarget, c, l, p); b && (post = POST_R)
    b, c, l = LpRmuLmuRp!(trTarget, c, l, p); b && (post = POST_R_T)

    # (8.9)
    b, c, l = LpRmSmLm!(target, c, l, p); b && (post = POST)
    b, c, l = LpRmSmLm!(tTarget, c, l, p); b && (post = POST_T)
    b, c, l = LpRmSmLm!(rTarget, c, l, p); b && (post = POST_R)
    b, c, l = LpRmSmLm!(trTarget, c, l, p); b && (post = POST_R_T)
    b, c, l = LpRmSmLm!(bTarget, c, l, p); b && (post = POST_B)
    b, c, l = LpRmSmLm!(btTarget, c, l, p); b && (post = POST_B_T)
    b, c, l = LpRmSmLm!(brTarget, c, l, p); b && (post = POST_B_R)
    b, c, l = LpRmSmLm!(btrTarget, c, l, p); b && (post = POST_B_R_T)

    # (8.10)
    b, c, l = LpRmSmRm!(target, c, l, p); b && (post = POST)
    b, c, l = LpRmSmRm!(tTarget, c, l, p); b && (post = POST_T)
    b, c, l = LpRmSmRm!(rTarget, c, l, p); b && (post = POST_R)
    b, c, l = LpRmSmRm!(trTarget, c, l, p); b && (post = POST_R_T)
    b, c, l = LpRmSmRm!(bTarget, c, l, p); b && (post = POST_B)
    b, c, l = LpRmSmRm!(btTarget, c, l, p); b && (post = POST_B_T)
    b, c, l = LpRmSmRm!(brTarget, c, l, p); b && (post = POST_B_R)
    b, c, l = LpRmSmRm!(btrTarget, c, l, p); b && (post = POST_B_R_T)

    # (8.11) C|Cpi/2 S Cpi/2|C
    b, c, l = LpRmSmLmRp!(target, c, l, p); b && (post = POST)
    b, c, l = LpRmSmLmRp!(tTarget, c, l, p); b && (post = POST_T)
    b, c, l = LpRmSmLmRp!(rTarget, c, l, p); b && (post = POST_R)
    b, c, l = LpRmSmLmRp!(trTarget, c, l, p); b && (post = POST_R_T)

    u = [carsegment2stepcontrol(p[i], r) for i in 1:l]
    if post == POST_T
        timeflip!(u)
    elseif post == POST_R
        reflect!(u)
    elseif post == POST_B
        backwards!(u)
    elseif post == POST_R_T
        reflect!(timeflip!(u))
    elseif post == POST_B_T
        backwards!(timeflip!(u))
    elseif post == POST_B_R
        backwards!(reflect!(u))
    elseif post == POST_B_R_T
        backwards!(reflect!(timeflip!(u)))
    end
    c*r, u
end

## Where the magic happens
function LpSpLp!{T}(tx::T, ty::T, tt::T, c::T, l::Int, path::CarPath{T})
    r, θ = R(tx - sin(tt), ty - 1 + cos(tt))
    u = r
    t = mod2pi(θ)
    v = mod2pi(tt - t)
    cnew = t + u + v
    c <= cnew && return false, c, l
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(0, u)
    path[3] = CarSegment{T}(1, v)
    true, cnew, 3
end

function LpSpRp!{T}(tx::T, ty::T, tt::T, c::T, l::Int, path::CarPath{T})
    r, θ = R(tx + sin(tt), ty - 1 - cos(tt))
    r*r < 4 && return false, c, l
    u = sqrt(r*r - 4)
    r1, θ1 = R(u, T(2))
    t = mod2pi(θ + θ1)
    v = mod2pi(t - tt)
    cnew = t + u + v
    c <= cnew && return false, c, l
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(0, u)
    path[3] = CarSegment{T}(-1, v)
    true, cnew, 3
end

function LpRmLp!{T}(tx::T, ty::T, tt::T, c::T, l::Int, path::CarPath{T})
    E = tx - sin(tt)
    N = ty + cos(tt) - 1
    E*E + N*N > 16 && return false, c, l
    r, θ = R(E, N)
    u = acos(1 - r*r/8)
    t = mod2pi(θ - u/2 + pi)
    v = mod2pi(pi - u/2 - θ + tt)
    u = -u
    cnew = t - u + v
    c <= cnew && return false, c, l
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(-1, u)
    path[3] = CarSegment{T}(1, v)
    true, cnew, 3
end
    
function LpRmLm!{T}(tx::T, ty::T, tt::T, c::T, l::Int, path::CarPath{T})
    E = tx - sin(tt)
    N = ty + cos(tt) - 1
    E*E + N*N > 16 && return false, c, l
    r, θ = R(E, N)
    u = acos(1 - r*r/8)
    t = mod2pi(θ - u/2 + pi)
    v = mod2pi(pi - u/2 - θ + tt) - 2*T(pi)
    u = -u
    cnew = t - u - v
    c <= cnew && return false, c, l
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(-1, u)
    path[3] = CarSegment{T}(1, v)
    true, cnew, 3
end
    
function LpRpuLmuRm!{T}(tx::T, ty::T, tt::T, c::T, l::Int, path::CarPath{T})
    E = tx + sin(tt)
    N = ty - cos(tt) - 1
    p = (2 + sqrt(E*E + N*N)) / 4
    (p < 0 || p > 1) && return false, c, l
    u = acos(p)
    t = mod2pi(Tau(u, -u, E, N))
    v = mod2pi(Omega(u, -u, E, N, tt)) - 2*T(pi)
    cnew = t + 2*u - v
    c <= cnew && return false, c, l
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(-1, u)
    path[3] = CarSegment{T}(1, -u)
    path[4] = CarSegment{T}(-1, v)
    true, cnew, 4
end

function LpRmuLmuRp!{T}(tx::T, ty::T, tt::T, c::T, l::Int, path::CarPath{T})
    E = tx + sin(tt)
    N = ty - cos(tt) - 1
    p = (20 - E*E - N*N) / 16
    (p < 0 || p > 1) && return false, c, l
    u = -acos(p)
    t = mod2pi(Tau(u, u, E, N))
    v = mod2pi(Omega(u, u, E, N, tt))
    cnew = t - 2*u + v
    c <= cnew && return false, c, l
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(-1, u)
    path[3] = CarSegment{T}(1, u)
    path[4] = CarSegment{T}(-1, v)
    true, cnew, 4
end

function LpRmSmLm!{T}(tx::T, ty::T, tt::T, c::T, l::Int, path::CarPath{T})
    E = tx - sin(tt)
    N = ty + cos(tt) - 1
    D, β = R(E, N)
    D < 2 && return false, c, l
    γ = acos(2/D)
    F = sqrt(D*D/4 - 1)
    t = mod2pi(pi + β - γ)
    u = 2 - 2*F
    u > 0 && return false, c, l
    v = mod2pi(-3*T(pi)/2 + γ + tt - β) - 2*T(pi)
    cnew = t + T(pi)/2 - u - v
    c <= cnew && return false, c, l
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(-1, -T(pi)/2)
    path[3] = CarSegment{T}(0, u)
    path[4] = CarSegment{T}(1, v)
    true, cnew, 4
end

function LpRmSmRm!{T}(tx::T, ty::T, tt::T, c::T, l::Int, path::CarPath{T})
    E = tx + sin(tt)
    N = ty - cos(tt) - 1
    D, β = R(E, N)
    D < 2 && return false, c, l
    t = mod2pi(β + T(pi)/2)
    u = 2 - D
    u > 0 && return false, c, l
    v = mod2pi(-pi - tt + β) - 2*T(pi)
    cnew = t + T(pi)/2 - u - v
    c <= cnew && return false, c, l
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(-1, -T(pi)/2)
    path[3] = CarSegment{T}(0, u)
    path[4] = CarSegment{T}(-1, v)
    true, cnew, 4
end
    
function LpRmSmLmRp!{T}(tx::T, ty::T, tt::T, c::T, l::Int, path::CarPath{T})
    E = tx + sin(tt)
    N = ty - cos(tt) - 1
    D, β = R(E, N)
    D < 2 && return false, c, l
    γ = acos(2/D)
    F = sqrt(D*D/4 - 1)
    t = mod2pi(pi + β - γ)
    u = 4 - 2*F
    u > 0 && return false, c, l
    v = mod2pi(pi + β - tt - γ)
    cnew = t + pi - u + v
    c <= cnew && return false, c, l
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(-1, -T(pi)/2)
    path[3] = CarSegment{T}(0, u)
    path[4] = CarSegment{T}(1, -T(pi)/2)
    path[5] = CarSegment{T}(-1, v)
    true, cnew, 5
end

for f in (:LpSpLp!, :LpSpRp!, :LpRmLp!, :LpRmLm!, :LpRpuLmuRm!, :LpRmuLmuRp!, :LpRmSmLm!, :LpRmSmRm!, :LpRmSmLmRp!)
    @eval $f{T}(s::SE2State{T}, c::T, l::Int, p::CarPath{T}) = $f(s.x[1], s.x[2], s.θ, c, l, p)
end