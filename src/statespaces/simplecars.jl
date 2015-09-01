export ReedsSheppMetricSpace, DubinsQuasiMetricSpace

### Simple Car Typedefs
## Basics
immutable CarSegment{T<:FloatingPoint}
    t::Int          # segment type (L = 1, S = 0, R = -1)               # TODO: t for type here and t for theta in SE2State is terribly confusing
    d::T            # segment length (radians for curved segments)
end
typealias CarPath{T} Vector{CarSegment{T}}
immutable CirclePoints{T<:FloatingPoint}  # TODO: cache CW and CCW for faster waypoint interpolation
    r::T
    dt::T
    pts::Vector{Vector2{T}}
end
CirclePoints{T}(r::T, N::Int) = CirclePoints(r, 2pi/N, [r*Vector2(cos(x), sin(x)) for x in linspace(0,2pi,N+1)[1:end-1]])

## Metrics
immutable ReedsSheppExact{T<:FloatingPoint} <: Metric
    r::T
    CP::CirclePoints{T}
    p::CarPath{T}    # scratch space to avoid extra memory allocation

    ReedsSheppExact(r::T, CP::CirclePoints{T}) = new(r, CP, Array(CarSegment{T}, 5))
end
ReedsSheppExact{T<:FloatingPoint}(r::T, CP::CirclePoints{T}) = ReedsSheppExact{T}(r, CP)
immutable DubinsExact{T<:FloatingPoint} <: QuasiMetric
    r::T
    CP::CirclePoints{T}
    p::CarPath{T}    # scratch space to avoid extra memory allocation

    DubinsExact(r::T, CP::CirclePoints{T}) = new(r, CP, Array(CarSegment{T}, 3))
end
DubinsExact{T<:FloatingPoint}(r::T, CP::CirclePoints{T}) = DubinsExact{T}(r, CP)

evaluate(RS::ReedsSheppExact, v::SE2State, w::SE2State) = reedsshepp(v, w, RS.r, RS.p)[1]
evaluate(D::DubinsExact, v::SE2State, w::SE2State) = dubins(v, w, D.r, D.p)[1]
defaultNN(RS::ReedsSheppExact, init::SE2State) = ArcLength_Pruned(typeof(init)[init], RS)  # TODO: abstract ALB_Metric
defaultNN(D::DubinsExact, init::SE2State) = QMArcLength_Pruned(typeof(init)[init], D)

## Metric(-ish) Space Instantiation 
ReedsSheppMetricSpace(r, lo = Vector2(0.,0.), hi = Vector2(1.,1.); res = 16) =
    SE2StateSpace(3, lo, hi, ReedsSheppExact(r, CirclePoints(r, res)), ExtractVector())
DubinsQuasiMetricSpace(r, lo = Vector2(0.,0.), hi = Vector2(1.,1.); res = 16) =
    SE2StateSpace(3, lo, hi, DubinsExact(r, CirclePoints(r, res)), ExtractVector())

### Waypoint Interpolation
## Full Waypoints
function waypoints{T}(v::SE2State{T}, s::CarSegment{T}, r::T = 1., cp::CirclePoints{T} = CirclePoints(r, 50))
    s.t == 0 && return [v, SE2State{T}(v.x+s.d*Vector2(cos(v.t), sin(v.t)), v.t)]
    center = v.x + sign(s.t)*Vector2(-r*sin(v.t), r*cos(v.t))
    turnpts = push!(cp.pts[1:1+ifloor(abs(s.d)/cp.dt)], Vector2(r*cos(abs(s.d)), r*sin(abs(s.d))))
    if s.t*s.d < 0
        for i in 1:length(turnpts)      # manual @devec
            turnpts[i] = Vector2(turnpts[i][1], -turnpts[i][2])
        end
    end
    pts = Vector2{T}[(center + sign(s.t)*rotate(p, v.t-pi/2)) for p in turnpts]
    thetas = push!([v.t + i*s.t*CP.dt for i in 0:ifloor(abs(s.d)/cp.dt)], v.t + s.t*s.d)  # TODO: this whole method is hacky/suboptimal for now, see CirclePoints note
    [SE2State{T}(p,t) for (p,t) in zip(pts, thetas)]
end
function waypoints{T}(v::SE2State{T}, p::CarPath{T}, r::T = 1., cp::CirclePoints{T} = CirclePoints(r, 50))
    pts = Array(SE2State{T}, 0)
    for s in p
        s_pts = waypoints(v, s, r, cp)
        append!(pts, s_pts[1:end-1])
        v = s_pts[end]
    end
    push!(pts, v)
end
function waypoints{T}(v::SE2State{T}, w::SE2State{T}, p::CarPath{T}, r::T = 1., cp::CirclePoints{T} = CirclePoints(r, 50), fix_w::Bool = false)
    pts = waypoints(v, p, r, cp)
    if fix_w
        error("TODO: this code is all temporary anyway, so I won't implement this right now.")
    else
        pts[end] = w
    end
    pts
end

## Special Cases for collision_waypoints and workspace_waypoints
# ExtractVector (Planar Position)
function workspace_waypoints{T}(v::SE2State{T}, s::CarSegment{T}, s2w::ExtractVector, r::T = 1., cp::CirclePoints{T} = CirclePoints(r, 50))
    s.t == 0 && return Vector2{T}[v.x, v.x+s.d*Vector2(cos(v.t), sin(v.t))]  # TODO: remove when [a, b] no longer concatenates
    center = v.x + sign(s.t)*Vector2(-r*sin(v.t), r*cos(v.t))
    turnpts = push!(cp.pts[1:1+ifloor(abs(s.d)/cp.dt)], Vector2(r*cos(abs(s.d)), r*sin(abs(s.d))))
    if s.t*s.d < 0
        for i in 1:length(turnpts)      # manual @devec
            turnpts[i] = Vector2(turnpts[i][1], -turnpts[i][2])
        end
    end
    Vector2{T}[(center + sign(s.t)*rotate(p, v.t-pi/2)) for p in turnpts]
end
function workspace_waypoints{T}(v::SE2State{T}, p::CarPath{T}, s2w::ExtractVector, r::T = 1., cp::CirclePoints{T} = CirclePoints(r, 50))
    pts = Array(Vector2{T}, 0)
    for s in p
        s_pts = workspace_waypoints(v, s, s2w, r, cp)
        append!(pts, s_pts[1:end-1])
        v = SE2State(s_pts[end], v.t + s.t*s.d)
    end
    push!(pts, v.x)
end
function workspace_waypoints{T}(v::SE2State{T}, w::SE2State{T}, p::CarPath{T}, s2w::ExtractVector, r::T = 1., cp::CirclePoints{T} = CirclePoints(r, 50), fix_w::Bool = false)
    pts = workspace_waypoints(v, p, s2w, r, cp)
    if fix_w
        wx0 = pts[end]
        segment_lengths = map(norm, diff(pts))
        scale_factors = cumsum(segment_lengths) ./ sum(segment_lengths)
        for i in 1:length(scale_factors)
            pts[i+1] = pts[i+1] + scale_factors[i]*(w.x - wx0)
        end
    else
        pts[end] = w.x
    end
    pts
end
workspace_waypoints{T}(v::SE2State{T}, w::SE2State{T}, RS::ReedsSheppExact{T}, s2w::State2Workspace, cp::CirclePoints{T} = CirclePoints(RS.r, 50)) =
    workspace_waypoints(v, w, reedsshepp(v, w, RS.r, RS.p)[2], s2w, RS.r, cp)
workspace_waypoints{T}(v::SE2State{T}, w::SE2State{T}, D::DubinsExact{T}, s2w::State2Workspace, cp::CirclePoints{T} = CirclePoints(D.r, 50)) =
    workspace_waypoints(v, w, dubins(v, w, D.r, D.p)[2], s2w, D.r, cp)
collision_waypoints{T}(v::SE2State{T}, w::SE2State{T}, RS::ReedsSheppExact{T}, s2w::State2Workspace, cp::CirclePoints{T} = RS.CP) =
    workspace_waypoints(v, w, reedsshepp(v, w, RS.r, RS.p)[2], s2w, RS.r, cp)
collision_waypoints{T}(v::SE2State{T}, w::SE2State{T}, D::DubinsExact{T}, s2w::State2Workspace, cp::CirclePoints{T} = D.CP) =
    workspace_waypoints(v, w, dubins(v, w, D.r, D.p)[2], s2w, D.r, cp)

## Discretized by Time
function time_waypoint{T}(v::SE2State{T}, s::CarSegment{T}, r::T, t::T)
    s.t == 0 && return v.x + min(t, s.d)*Vector2(cos(v.t), sin(v.t))
    center = v.x + sign(s.t)*Vector2(-r*sin(v.t), r*cos(v.t))
    th = min(t/r, abs(s.d))
    turnpt = s.t*s.d > 0 ? Vector2(r*cos(th), r*sin(th)) : Vector2(r*cos(th), -r*sin(th))
    center + sign(s.t)*rotate(turnpt, v.t-pi/2)
end
function time_waypoint{T}(v::SE2State{T}, w::SE2State{T}, p::CarPath{T}, r::T, t::T)
    pt = w
    total_time = 0    # should be similar to the cost of p, up to endpoint-fixing stuff
    for s in p
        segment_time = s.t == 0 ? abs(s.d) : r*abs(s.d)
        if total_time <= t < total_time + segment_time
            pt = time_waypoint(v, s, r, t - total_time)
        end
        total_time = total_time + segment_time
        v = SE2State(time_waypoint(v, s, r, Inf), v.t + s.t*s.d)
    end
    t > total_time && return w.x
    pt + t/total_time*(w.x - v.x)
end

### Dubins Steering Nuts and Bolts
function dubinsLSL!{T}(d::T, a::T, b::T, c::T, path::CarPath{T})
    ca, sa, cb, sb = cos(a), sin(a), cos(b), sin(b)
    tmp = 2. + d*d - 2.*(ca*cb +sa*sb - d*(sa - sb))
    tmp < 0. && return c
    th = atan2(cb - ca, d + sa - sb)
    t = mod2pi(-a + th)
    p = sqrt(max(tmp, 0.))
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
    tmp = 2. + d*d - 2.*(ca*cb + sa*sb - d*(sb - sa))
    tmp < 0. && return c
    th = atan2(ca - cb, d - sa + sb)
    t = mod2pi(a - th)
    p = sqrt(max(tmp, 0.))
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
    tmp = d * d - 2. + 2. * (ca*cb + sa*sb - d * (sa + sb))
    tmp < 0. && return c
    p = sqrt(max(tmp, 0.))
    th = atan2(ca + cb, d - sa - sb) - atan2(2., p)
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
    tmp = -2. + d * d + 2. * (ca*cb + sa*sb + d * (sa + sb))
    tmp < 0. && return c
    p = sqrt(max(tmp, 0.))
    th = atan2(-ca - cb, d + sa + sb) - atan2(-2., p)
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
    tmp = .125 * (6. - d * d  + 2. * (ca*cb + sa*sb + d * (sa - sb)))
    abs(tmp) >= 1. && return c
    p = 2pi - acos(tmp)
    th = atan2(ca - cb, d - sa + sb)
    t = mod2pi(a - th + .5 * p)
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
    tmp = .125 * (6. - d * d  + 2. * (ca*cb + sa*sb - d * (sa - sb)))
    abs(tmp) >= 1. && return c
    p = 2pi - acos(tmp)
    th = atan2(-ca + cb, d + sa - sb)
    t = mod2pi(-a + th + .5 * p)
    q = mod2pi(b - a - t + p)
    cnew = t + p + q
    c <= cnew && return c
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(-1, p)
    path[3] = CarSegment{T}(1, q)
    cnew
end

function dubins{T}(s1::SE2State{T}, s2::SE2State{T}, r = one(T), pmin = Array(CarSegment{T}, 3))
    const dubinsTypes = (dubinsLSL!, dubinsRSR!, dubinsRSL!, dubinsLSR!, dubinsRLR!, dubinsLRL!)
    v = (s2.x - s1.x) / r
    d = norm(v)
    th = atan2(v)
    a = mod2pi(s1.t - th)
    b = mod2pi(s2.t - th)

    cmin = inf(T)
    for dubinsType! in dubinsTypes
        cmin = dubinsType!(d, a, b, cmin, pmin)
    end
    cmin*r, rescale!(pmin, r)
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

### Reeds-Shepp Steering Nuts and Bolts
## Utilities
R(x, y) = sqrt(x*x + y*y), atan2(y, x)
function M(t)
    m = mod2pi(t)
    m > pi ? m - 2pi : m
end
function Tau(u, v, E, N)
    delta = M(u - v)
    A = sin(u) - sin(delta)
    B = cos(u) - cos(delta) - 1.
    r, θ = R(E*A + N*B, N*A - E*B)
    t = 2*cos(delta) - 2*cos(v) - 2*cos(u) + 3.
    t < 0 ? M(θ + pi) : M(θ)
end
Omega(u, v, E, N, t) = M(Tau(u, v, E, N) - u + v - t)
timeflip{T}(s::SE2State{T}) = SE2State(-s.x[1], s.x[2], -s.t)
reflect{T}(s::SE2State{T}) = SE2State(s.x[1], -s.x[2], -s.t)
backwards{T}(s::SE2State{T}) = SE2State(s.x[1]*cos(s.t) + s.x[2]*sin(s.t), s.x[1]*sin(s.t) - s.x[2]*cos(s.t), s.t)
function timeflip!{T}(p::CarPath{T})
    for i in 1:length(p)
        p[i] = CarSegment{T}(p[i].t, -p[i].d)
    end
    p
end
function reflect!{T}(p::CarPath{T})
    for i in 1:length(p)
        p[i] = CarSegment{T}(-p[i].t, p[i].d)
    end
    p
end
backwards!{T}(p::CarPath{T}) = reverse!(p)
function rescale!{T}(p::CarPath{T}, r::T)
    for i in 1:length(p)
        if p[i].t == 0
            p[i] = CarSegment{T}(0, r*p[i].d)
        end
    end
    p
end

# TODO: something about an enum type in Julia v0.4?
const POST, POST_T, POST_R, POST_B, POST_R_T, POST_B_T, POST_B_R, POST_B_R_T = 0, 1, 2, 3, 4, 5, 6, 7

## Gotta check 'em all
function reedsshepp{T}(s1::SE2State{T}, s2::SE2State{T}, r = one(T), p = Array(CarSegment{T}, 5))
    dx, dy = (s2.x - s1.x) / r
    ct, st = cos(s1.t), sin(s1.t)
    target = SE2State(dx*ct + dy*st, -dx*st + dy*ct, mod2pi(s2.t - s1.t))

    tTarget = timeflip(target)
    rTarget = reflect(target)
    trTarget = reflect(tTarget)
    bTarget = backwards(target)
    btTarget = timeflip(bTarget)
    brTarget = reflect(bTarget)
    btrTarget = reflect(btTarget)
    
    c, l = inf(T), 0
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

    if post == POST
        p = p[1:l]
    elseif post == POST_T
        p = timeflip!(p[1:l])
    elseif post == POST_R
        p = reflect!(p[1:l])
    elseif post == POST_B
        p = backwards!(p[1:l])
    elseif post == POST_R_T
        p = reflect!(timeflip!(p[1:l]))
    elseif post == POST_B_T
        p = backwards!(timeflip!(p[1:l]))
    elseif post == POST_B_R
        p = backwards!(reflect!(p[1:l]))
    else
        p = backwards!(reflect!(timeflip!(p[1:l])))
    end
    c*r, rescale!(p, r)
end

## Where the magic happens
function LpSpLp!{T}(tx::T, ty::T, tt::T, c::T, l::Int, path::CarPath{T})
    r, θ = R(tx - sin(tt), ty - 1. + cos(tt))
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
    r, θ = R(tx + sin(tt), ty - 1. - cos(tt))
    r*r < 4. && return false, c, l
    u = sqrt(r*r - 4.)
    r1, θ1 = R(u, 2.)
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
    N = ty + cos(tt) - 1.
    E*E + N*N > 16. && return false, c, l
    r, θ = R(E, N)
    u = acos(1. - r*r/8)
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
    N = ty + cos(tt) - 1.
    E*E + N*N > 16. && return false, c, l
    r, θ = R(E, N)
    u = acos(1. - r*r/8)
    t = mod2pi(θ - u/2 + pi)
    v = mod2pi(pi - u/2 - θ + tt) - 2pi
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
    N = ty - cos(tt) - 1.
    p = (2. + sqrt(E*E + N*N)) / 4.
    (p < 0 || p > 1) && return false, c, l
    u = acos(p)
    t = mod2pi(Tau(u, -u, E, N))
    v = mod2pi(Omega(u, -u, E, N, tt)) - 2pi
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
    N = ty - cos(tt) - 1.
    p = (20. - E*E - N*N) / 16.
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
    N = ty + cos(tt) - 1.
    D, β = R(E, N)
    D < 2 && return false, c, l
    γ = acos(2/D)
    F = sqrt(D*D/4 - 1.)
    t = mod2pi(pi + β - γ)
    u = 2. - 2*F
    u > 0 && return false, c, l
    v = mod2pi(-3pi/2 + γ + tt - β) - 2pi
    cnew = t + pi/2 - u - v
    c <= cnew && return false, c, l
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(-1, -pi/2)
    path[3] = CarSegment{T}(0, u)
    path[4] = CarSegment{T}(1, v)
    true, cnew, 4
end

function LpRmSmRm!{T}(tx::T, ty::T, tt::T, c::T, l::Int, path::CarPath{T})
    E = tx + sin(tt)
    N = ty - cos(tt) - 1.
    D, β = R(E, N)
    D < 2 && return false, c, l
    t = mod2pi(β + pi/2)
    u = 2. - D
    u > 0 && return false, c, l
    v = mod2pi(-pi - tt + β) - 2pi
    cnew = t + pi/2 - u - v
    c <= cnew && return false, c, l
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(-1, -pi/2)
    path[3] = CarSegment{T}(0, u)
    path[4] = CarSegment{T}(-1, v)
    true, cnew, 4
end
    
function LpRmSmLmRp!{T}(tx::T, ty::T, tt::T, c::T, l::Int, path::CarPath{T})
    E = tx + sin(tt)
    N = ty - cos(tt) - 1.
    D, β = R(E, N)
    D < 2 && return false, c, l
    γ = acos(2/D)
    F = sqrt(D*D/4 - 1.)
    t = mod2pi(pi + β - γ)
    u = 4. - 2*F
    u > 0 && return false, c, l
    v = mod2pi(pi + β - tt - γ)
    cnew = t + pi - u + v
    c <= cnew && return false, c, l
    path[1] = CarSegment{T}(1, t)
    path[2] = CarSegment{T}(-1, -pi/2)
    path[3] = CarSegment{T}(0, u)
    path[4] = CarSegment{T}(1, -pi/2)
    path[5] = CarSegment{T}(-1, v)
    true, cnew, 5
end

for f in (:LpSpLp!, :LpSpRp!, :LpRmLp!, :LpRmLm!, :LpRpuLmuRm!, :LpRmuLmuRp!, :LpRmSmLm!, :LpRmSmRm!, :LpRmSmLmRp!)
    @eval $f{T}(s::SE2State{T}, c::T, l::Int, p::CarPath{T}) = $f(s.x[1], s.x[2], s.t, c, l, p)
end