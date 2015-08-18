### Simple Car State Spaces
abstract SimpleCarStateSpace <: StateSpace
immutable CarSegment{T<:FloatingPoint}
    t::Int          # segment type (L = 1, S = 0, R = -1)
    d::T            # segment length (radians for curved segments)
end
typealias CarPath{T} Vector{CarSegment{T}}

### Waypoint Interpolation
## Arbitrary Spacing (for CC)
immutable CirclePoints{T<:FloatingPoint}
    r::T
    dt::T
    pts::Vector{Vector2{T}}
end
CirclePoints{T}(r::T, N::Int) = CirclePoints(r, 2pi/N, [r*Vector2(cos(x), sin(x)) for x in linspace(0,2pi,N+1)[1:end-1]])

function waypoints{T}(v::SE2State{T}, s::CarSegment{T}, r::T = 1., cp::CirclePoints{T} = CirclePoints(r, 16))
    s.t == 0 && return Vector2{T}[v.x, v.x+s.d*Vector2(cos(v.t), sin(v.t))]
    center = v.x + sign(s.t)*Vector2(-r*sin(v.t), r*cos(v.t))
    turnpts = push!(cp.pts[1:1+ifloor(abs(s.d)/cp.dt)], Vector2(r*cos(abs(s.d)), r*sin(abs(s.d))))
    if s.t*s.d < 0
        for i in 1:length(turnpts)      # manual @devec
            turnpts[i] = Vector2(turnpts[i][1], -turnpts[i][2])
        end
    end
    [(center + sign(s.t)*rotate(p, v.t-pi/2)) for p in turnpts]
end
function waypoints{T}(v::SE2State{T}, p::CarPath{T}, r::T = 1., cp::CirclePoints{T} = CirclePoints(r, 16))
    pts = Array(Vector2{T}, 0)
    for s in p
        s_pts = waypoints(v, s, r, cp)
        append!(pts, s_pts[1:end-1])
        v = SE2State(s_pts[end], v.t + s.t*s.d)
    end
    push!(pts, v.x)
end
function waypoints{T}(v::SE2State{T}, w::SE2State{T}, p::CarPath{T}, r::T = 1., cp::CirclePoints{T} = CirclePoints(r, 16), fix_w::Bool = false)
    pts = waypoints(v, p, r, cp)
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
function dubinsLSL{T}(d::T, a::T, b::T)
    ca, sa, cb, sb = cos(a), sin(a), cos(b), sin(b)
    tmp = 2. + d*d - 2.*(ca*cb +sa*sb - d*(sa - sb))
    tmp < 0. && return (Inf, CarSegment{T}[])
    th = atan2(cb - ca, d + sa - sb)
    t = mod2pi(-a + th)
    p = sqrt(max(tmp, 0.))
    q = mod2pi(b - th)
    (t + p + q, [CarSegment{T}(1, t), CarSegment{T}(0, p), CarSegment{T}(1, q)])
end

function dubinsRSR{T}(d::T, a::T, b::T)
    ca, sa, cb, sb = cos(a), sin(a), cos(b), sin(b)
    tmp = 2. + d*d - 2.*(ca*cb + sa*sb - d*(sb - sa))
    tmp < 0. && return (Inf, CarSegment{T}[])
    th = atan2(ca - cb, d - sa + sb)
    t = mod2pi(a - th)
    p = sqrt(max(tmp, 0.))
    q = mod2pi(-b + th)
    (t + p + q, [CarSegment{T}(-1, t), CarSegment{T}(0, p), CarSegment{T}(-1, q)])
end

function dubinsRSL{T}(d::T, a::T, b::T)
    ca, sa, cb, sb = cos(a), sin(a), cos(b), sin(b)
    tmp = d * d - 2. + 2. * (ca*cb + sa*sb - d * (sa + sb))
    tmp < 0. && return (Inf, CarSegment{T}[])
    p = sqrt(max(tmp, 0.))
    th = atan2(ca + cb, d - sa - sb) - atan2(2., p)
    t = mod2pi(a - th)
    q = mod2pi(b - th)
    (t + p + q, [CarSegment{T}(-1, t), CarSegment{T}(0, p), CarSegment{T}(1, q)])
end

function dubinsLSR{T}(d::T, a::T, b::T)
    ca, sa, cb, sb = cos(a), sin(a), cos(b), sin(b)
    tmp = -2. + d * d + 2. * (ca*cb + sa*sb + d * (sa + sb))
    tmp < 0. && return (Inf, CarSegment{T}[])
    p = sqrt(max(tmp, 0.))
    th = atan2(-ca - cb, d + sa + sb) - atan2(-2., p)
    t = mod2pi(-a + th)
    q = mod2pi(-b + th)
    (t + p + q, [CarSegment{T}(1, t), CarSegment{T}(0, p), CarSegment{T}(-1, q)])
end

function dubinsRLR{T}(d::T, a::T, b::T)
    ca, sa, cb, sb = cos(a), sin(a), cos(b), sin(b)
    tmp = .125 * (6. - d * d  + 2. * (ca*cb + sa*sb + d * (sa - sb)))
    abs(tmp) >= 1. && return (Inf, CarSegment{T}[])
    p = 2pi - acos(tmp)
    th = atan2(ca - cb, d - sa + sb)
    t = mod2pi(a - th + .5 * p)
    q = mod2pi(a - b - t + p)
    (t + p + q, [CarSegment{T}(-1, t), CarSegment{T}(1, p), CarSegment{T}(-1, q)])
end

function dubinsLRL{T}(d::T, a::T, b::T)
    ca, sa, cb, sb = cos(a), sin(a), cos(b), sin(b)
    tmp = .125 * (6. - d * d  + 2. * (ca*cb + sa*sb - d * (sa - sb)))
    abs(tmp) >= 1. && return (Inf, CarSegment{T}[])
    p = 2pi - acos(tmp)
    th = atan2(-ca + cb, d + sa - sb)
    t = mod2pi(-a + th + .5 * p)
    q = mod2pi(b - a - t + p)
    (t + p + q, [CarSegment{T}(1, t), CarSegment{T}(-1, p), CarSegment{T}(1, q)])
end

function dubins{T}(s1::SE2State{T}, s2::SE2State{T})
    const dubinsTypes = (dubinsLSL, dubinsRSR, dubinsRSL, dubinsLSR, dubinsRLR, dubinsLRL)
    v = s2.x - s1.x
    d = norm(v)
    th = atan2(v)
    a = mod2pi(s1.t - th)
    b = mod2pi(s2.t - th)

    cmin, pmin = Inf, CarSegment{T}[]
    for dubinsType in dubinsTypes
        c, p = dubinsType(d, a, b)
        if c < cmin
            cmin, pmin = c, p
        end
    end
    cmin, pmin
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

### Reeds-Shepp Steering Nuts and Bolts (TODO: if I'm really worried about speed, everything should get '!'s)
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
timeflip{T}(p::CarPath{T}) = [CarSegment{T}(s.t, -s.d) for s in p]
reflect{T}(p::CarPath{T}) = [CarSegment{T}(-s.t, s.d) for s in p]
backwards{T}(p::CarPath{T}) = reverse(p)

## Gotta check 'em all
function reedsshepp{T}(s1::SE2State{T}, s2::SE2State{T})
    dx, dy = s2.x - s1.x
    c, s = cos(s1.t), sin(s1.t)
    target = SE2State(dx*c + dy*s, -dx*s + dy*c, s2.t - s1.t)

    tTarget = timeflip(target)
    rTarget = reflect(target)
    trTarget = reflect(tTarget)
    bTarget = backwards(target)
    btTarget = timeflip(bTarget)
    brTarget = reflect(bTarget)
    btrTarget = reflect(btTarget)

    cmin, pmin = Inf, CarSegment{T}[]
    # Only Tim Wheeler knows what the hell is happening down below
    # (8.1) C S C
    c, p = LpSpLp(target); if c < cmin; cmin, pmin = c, p; end
    c, p = LpSpLp(tTarget); if c < cmin; cmin, pmin = c, timeflip(p); end
    c, p = LpSpLp(rTarget); if c < cmin; cmin, pmin = c, reflect(p); end
    c, p = LpSpLp(trTarget); if c < cmin; cmin, pmin = c, reflect(timeflip(p)); end

    # (8.2) C S C
    c, p = LpSpRp(target); if c < cmin; cmin, pmin = c, p; end
    c, p = LpSpRp(tTarget); if c < cmin; cmin, pmin = c, timeflip(p); end
    c, p = LpSpRp(rTarget); if c < cmin; cmin, pmin = c, reflect(p); end
    c, p = LpSpRp(trTarget); if c < cmin; cmin, pmin = c, reflect(timeflip(p)); end

    # (8.3) C|C|C
    c, p = LpRmLp(target); if c < cmin; cmin, pmin = c, p; end                                   # L+R-L+
    # c, p = LpRmLp(tTarget); if c < cmin; cmin, pmin = c, timeflip(p); end                      # L-R+L- (redundant)
    c, p = LpRmLp(rTarget); if c < cmin; cmin, pmin = c, reflect(p); end                         # R+L-R+
    # c, p = LpRmLp(trTarget); if c < cmin; cmin, pmin = c, reflect(timeflip(p)); end            # R-L+R- (redundant)

    # (8.4) C|C C
    c, p = LpRmLm(target); if c < cmin; cmin, pmin = c, p; end                                   # L+R-L-
    c, p = LpRmLm(tTarget); if c < cmin; cmin, pmin = c, timeflip(p); end                        # L-R+L+
    c, p = LpRmLm(rTarget); if c < cmin; cmin, pmin = c, reflect(p); end                         # R+L-R-
    c, p = LpRmLm(trTarget); if c < cmin; cmin, pmin = c, reflect(timeflip(p)); end              # R-L+R+
    c, p = LpRmLm(bTarget); if c < cmin; cmin, pmin = c, backwards(p); end                       # L-R-L+
    c, p = LpRmLm(btTarget); if c < cmin; cmin, pmin = c, backwards(timeflip(p)); end            # L+R+L-
    c, p = LpRmLm(brTarget); if c < cmin; cmin, pmin = c, backwards(reflect(p)); end             # R-L-R+
    c, p = LpRmLm(btrTarget); if c < cmin; cmin, pmin = c, backwards(reflect(timeflip(p))); end  # R+L+R-

    # (8.7) C Cu|Cu C
    c, p = LpRpuLmuRm(target); if c < cmin; cmin, pmin = c, p; end
    c, p = LpRpuLmuRm(tTarget); if c < cmin; cmin, pmin = c, timeflip(p); end
    c, p = LpRpuLmuRm(rTarget); if c < cmin; cmin, pmin = c, reflect(p); end
    c, p = LpRpuLmuRm(trTarget); if c < cmin; cmin, pmin = c, reflect(timeflip(p)); end

    # (8.8) C|Cu Cu|C
    c, p = LpRmuLmuRp(target); if c < cmin; cmin, pmin = c, p; end
    c, p = LpRmuLmuRp(tTarget); if c < cmin; cmin, pmin = c, timeflip(p); end
    c, p = LpRmuLmuRp(rTarget); if c < cmin; cmin, pmin = c, reflect(p); end
    c, p = LpRmuLmuRp(trTarget); if c < cmin; cmin, pmin = c, reflect(timeflip(p)); end

    # (8.9)
    c, p = LpRmSmLm(target); if c < cmin; cmin, pmin = c, p; end
    c, p = LpRmSmLm(tTarget); if c < cmin; cmin, pmin = c, timeflip(p); end
    c, p = LpRmSmLm(rTarget); if c < cmin; cmin, pmin = c, reflect(p); end
    c, p = LpRmSmLm(trTarget); if c < cmin; cmin, pmin = c, reflect(timeflip(p)); end
    c, p = LpRmSmLm(bTarget); if c < cmin; cmin, pmin = c, backwards(p); end
    c, p = LpRmSmLm(btTarget); if c < cmin; cmin, pmin = c, backwards(timeflip(p)); end
    c, p = LpRmSmLm(brTarget); if c < cmin; cmin, pmin = c, backwards(reflect(p)); end
    c, p = LpRmSmLm(btrTarget); if c < cmin; cmin, pmin = c, backwards(reflect(timeflip(p))); end

    # (8.10)
    c, p = LpRmSmRm(target); if c < cmin; cmin, pmin = c, p; end
    c, p = LpRmSmRm(tTarget); if c < cmin; cmin, pmin = c, timeflip(p); end
    c, p = LpRmSmRm(rTarget); if c < cmin; cmin, pmin = c, reflect(p); end
    c, p = LpRmSmRm(trTarget); if c < cmin; cmin, pmin = c, reflect(timeflip(p)); end
    c, p = LpRmSmRm(bTarget); if c < cmin; cmin, pmin = c, backwards(p); end
    c, p = LpRmSmRm(btTarget); if c < cmin; cmin, pmin = c, backwards(timeflip(p)); end
    c, p = LpRmSmRm(brTarget); if c < cmin; cmin, pmin = c, backwards(reflect(p)); end
    c, p = LpRmSmRm(btrTarget); if c < cmin; cmin, pmin = c, backwards(reflect(timeflip(p))); end

    # (8.11) C|Cpi/2 S Cpi/2|C
    c, p = LpRmSmLmRp(target); if c < cmin; cmin, pmin = c, p; end
    c, p = LpRmSmLmRp(tTarget); if c < cmin; cmin, pmin = c, timeflip(p); end
    c, p = LpRmSmLmRp(rTarget); if c < cmin; cmin, pmin = c, reflect(p); end
    c, p = LpRmSmLmRp(trTarget); if c < cmin; cmin, pmin = c, reflect(timeflip(p)); end

    cmin, pmin
end

## Where the magic happens
function LpSpLp{T}(tx::T, ty::T, tt::T)
    r, θ = R(tx - sin(tt), ty - 1. + cos(tt))
    u = r
    t = mod2pi(θ)
    v = mod2pi(tt - t)
    (t + u + v, [CarSegment{T}(1, t), CarSegment{T}(0, u), CarSegment{T}(1, v)])
end

function LpSpRp{T}(tx::T, ty::T, tt::T)
    r, θ = R(tx + sin(tt), ty - 1. - cos(tt))
    r*r < 4. && return (Inf, CarSegment{T}[])
    u = sqrt(r*r - 4.)
    r1, θ1 = R(u, 2.)
    t = mod2pi(θ + θ1)
    v = mod2pi(t - tt)
    (t + u + v, [CarSegment{T}(1, t), CarSegment{T}(0, u), CarSegment{T}(-1, v)])
end

function LpRmLp{T}(tx::T, ty::T, tt::T)
    E = tx - sin(tt)
    N = ty + cos(tt) - 1.
    E*E + N*N > 16. && return (Inf, CarSegment{T}[])
    r, θ = R(E, N)
    u = acos(1. - r*r/8)
    t = mod2pi(θ - u/2 + pi)
    v = mod2pi(pi - u/2 - θ + tt)
    u = -u
    (t - u + v, [CarSegment{T}(1, t), CarSegment{T}(-1, u), CarSegment{T}(1, v)])
end

function LpRmLm{T}(tx::T, ty::T, tt::T)
    E = tx - sin(tt)
    N = ty + cos(tt) - 1.
    E*E + N*N > 16. && return (Inf, CarSegment{T}[])
    r, θ = R(E, N)
    u = acos(1. - r*r/8)
    t = mod2pi(θ - u/2 + pi)
    v = mod2pi(pi - u/2 - θ + tt) - 2pi
    u = -u
    (t - u - v, [CarSegment{T}(1, t), CarSegment{T}(-1, u), CarSegment{T}(1, v)])
end

function LpRpuLmuRm{T}(tx::T, ty::T, tt::T)
    E = tx + sin(tt)
    N = ty - cos(tt) - 1.
    p = (2. + sqrt(E*E + N*N)) / 4.
    (p < 0 || p > 1) && return (Inf, CarSegment{T}[])
    u = acos(p)
    t = mod2pi(Tau(u, -u, E, N))
    v = mod2pi(Omega(u, -u, E, N, tt)) - 2pi
    (t + 2*u - v, [CarSegment{T}(1, t), CarSegment{T}(-1, u), CarSegment{T}(1, -u), CarSegment{T}(-1, v)])
end

function LpRmuLmuRp{T}(tx::T, ty::T, tt::T)
    E = tx + sin(tt)
    N = ty - cos(tt) - 1.
    p = (20. - E*E - N*N) / 16.
    (p < 0 || p > 1) && return (Inf, CarSegment{T}[])
    u = -acos(p)
    t = mod2pi(Tau(u, u, E, N))
    v = mod2pi(Omega(u, u, E, N, tt))
    (t - 2*u + v, [CarSegment{T}(1, t), CarSegment{T}(-1, u), CarSegment{T}(1, u), CarSegment{T}(-1, v)])
end

function LpRmSmLm{T}(tx::T, ty::T, tt::T)
    E = tx - sin(tt)
    N = ty + cos(tt) - 1.
    D, β = R(E, N)
    D < 2 && return (Inf, CarSegment{T}[])
    γ = acos(2/D)
    F = sqrt(D*D/4 - 1.)
    t = mod2pi(pi + β - γ)
    u = 2. - 2*F
    u > 0 && return (Inf, CarSegment{T}[])
    v = mod2pi(-3pi/2 + γ + tt - β) - 2pi
    (t + pi/2 - u - v, [CarSegment{T}(1, t), CarSegment{T}(-1, -pi/2), CarSegment{T}(0, u), CarSegment{T}(1, v)])
end

function LpRmSmRm{T}(tx::T, ty::T, tt::T)
    E = tx + sin(tt)
    N = ty - cos(tt) - 1.
    D, β = R(E, N)
    D < 2 && return (Inf, CarSegment{T}[])
    t = mod2pi(β + pi/2)
    u = 2. - D
    u > 0 && return (Inf, CarSegment{T}[])
    v = mod2pi(-pi - tt + β) - 2pi
    (t + pi/2 - u - v, [CarSegment{T}(1, t), CarSegment{T}(-1, -pi/2), CarSegment{T}(0, u), CarSegment{T}(-1, v)])
end

function LpRmSmLmRp{T}(tx::T, ty::T, tt::T)
    E = tx + sin(tt)
    N = ty - cos(tt) - 1.
    D, β = R(E, N)
    D < 2 && return (Inf, CarSegment{T}[])
    γ = acos(2/D)
    F = sqrt(D*D/4 - 1.)
    t = mod2pi(pi + β - γ)
    u = 4. - 2*F
    u > 0 && return (Inf, CarSegment{T}[])
    v = mod2pi(pi + β - tt - γ)
    (t + pi - u + v, [CarSegment{T}(1, t), CarSegment{T}(-1, -pi/2), CarSegment{T}(0, u), CarSegment{T}(1, -pi/2), CarSegment{T}(-1, v)])
end

for f in (:LpSpLp, :LpSpRp, :LpRmLp, :LpRmLm, :LpRpuLmuRm, :LpRmuLmuRp, :LpRmSmLm, :LpRmSmRm, :LpRmSmLmRp)
    @eval $f(s::SE2State) = $f(s.x[1], s.x[2], s.t)
end