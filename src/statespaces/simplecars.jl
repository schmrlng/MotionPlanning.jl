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