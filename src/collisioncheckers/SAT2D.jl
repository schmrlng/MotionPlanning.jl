export Shape2D, Circle, Polygon, Box2D, Line, Compound2D, colliding, colliding_ends_free

import Base.atan2

# ---------- Vector Utilities ----------

norm2{T<:FloatingPoint}(v::Vector2{T}) = dot(v,v)
perp{T<:FloatingPoint}(v::Vector2{T}) = Vector2{T}(v[2],-v[1])
project{T<:FloatingPoint}(v1::Vector2{T}, v2::Vector2{T}) = (dot(v1,v2) / norm2(v2)) * v2
projectN{T<:FloatingPoint}(v1::Vector2{T}, v2::Vector2{T}) = dot(v1,v2) * v2
reflect{T<:FloatingPoint}(v1::Vector2{T}, ax::Vector2{T}) = v1 - 2*project(v1,ax)
reflectN{T<:FloatingPoint}(v1::Vector2{T}, ax::Vector2{T}) = v1 - 2*projectN(v1,ax)
function rotate{T<:FloatingPoint}(v::Vector2{T}, theta::T)
    c = cos(theta)
    s = sin(theta)
    Vector2{T}(v[1]*c - v[2]*s, v[1]*s + v[2]*c)
end
atan2{T<:FloatingPoint}(v::Vector2{T}) = atan2(v[2], v[1])
function projectNextrema{T<:FloatingPoint}(pts::AbstractVector{Vector2{T}}, n::Vector2{T})
    dmin = inf(T)
    dmax = -inf(T)
    for p in pts
        d = dot(p,n)
        d < dmin && (dmin = d)
        d > dmax && (dmax = d)
    end
    Vector2{T}(dmin, dmax)
end
function voronoi_region{T<:FloatingPoint}(p::Vector2{T}, v::Vector2{T})
    d = dot(p,v)
    d < 0 && return -1
    d > norm2(v) && return 1
    return 0
end
overlapping{T<:FloatingPoint}(I1::Vector2{T}, I2::Vector2{T}) = I1[1] <= I2[2] && I2[1] <= I1[2]
inrange{T<:FloatingPoint}(x::T, I::Vector2{T}) = I[1] <= x <= I[2]
minmaxV{T<:FloatingPoint}(x::T, y::T) = x < y ? Vector2{T}(x,y) : Vector2{T}(y,x)  # Base.minmax returns tuple 
wrap1(i,N) = i < 1 ? N : i > N ? 1 : i      # marginally faster than mod1

# ---------- Shape Definitions ----------

abstract Shape2D{T<:FloatingPoint}

type Circle{T} <: Shape2D{T}
    c::Vector2{T}
    r::T
    xrange::Vector2{T}
    yrange::Vector2{T}

    function Circle(c::Vector2{T},r::T,xrange::Vector2{T},yrange::Vector2{T})
        r <= 0 && error("Radius must be positive")
        new(c,r,xrange,yrange)
    end
end
Circle{T}(c::AbstractVector{T}, r::T) = Circle{T}(Vector2{T}(c...), r,
                                                  Vector2{T}(c[1]-r, c[1]+r),
                                                  Vector2{T}(c[2]-r, c[2]+r))
projectNextrema{T}(C::Circle{T}, n::Vector2{T}) = (d = dot(C.c,n); Vector2{T}(d-C.r, d+C.r))

type Polygon{T} <: Shape2D{T}
    points::Vector{Vector2{T}}
    edges::Vector{Vector2{T}}
    normals::Vector{Vector2{T}}
    xrange::Vector2{T}
    yrange::Vector2{T}
    nextrema::Vector{Vector2{T}}

    function Polygon(points::Vector{Vector2{T}})
        N = length(points)
        N < 3 && error("Polygons need at least 3 points! Try Line?")
        sum([(points[wrap1(i+1,N)][1] - points[i][1])*(points[wrap1(i+1,N)][2] + points[i][2]) for i in 1:N]) > 0 && reverse!(points)
        edges = Vector2{T}[points[wrap1(i+1,N)] - points[i] for i in 1:N]
        normals = Vector2{T}[unit(perp(g)) for g in edges]
        any(-pi .<= diff([T[atan2(n) for n in normals], atan2(normals[1])]) .<= 0) && error("Polygon needs to be convex")
        xrange = Vector2{T}(extrema([points[i][1] for i in 1:N])...)
        yrange = Vector2{T}(extrema([points[i][2] for i in 1:N])...)
        nextrema = Vector2{T}[projectNextrema(points, normals[i]) for i in 1:N]
        new(points, edges, normals, xrange, yrange, nextrema)
    end
end
Polygon{T}(points::Vector{Vector{T}}) = Polygon{T}([Vector2(p) for p in points])
Polygon{T}(points::Vector{Vector2{T}}) = Polygon{T}(points)   # so Polygon(pts) works without {T}
Polygon{T}(points::Vector{(T,T)}) = Polygon{T}([Vector2(p...) for p in points])
Box2D{T}(xr::AbstractVector{T}, yr::AbstractVector{T}) = Polygon{T}(Vector2{T}[Vector2(xr[1], yr[1]),
                                                                               Vector2(xr[2], yr[1]),
                                                                               Vector2(xr[2], yr[2]),
                                                                               Vector2(xr[1], yr[2])])
projectNextrema{T}(P::Polygon{T}, n::Vector2{T}) = projectNextrema(P.points, n)

type Line{T} <: Shape2D{T}   # special case of Polygon; only 1 edge/normal
    v::Vector2{T}
    w::Vector2{T}
    edge::Vector2{T}
    normal::Vector2{T}
    xrange::Vector2{T}
    yrange::Vector2{T}
    ndotv::T

    function Line(v::Vector2{T}, w::Vector2{T})
        edge = w - v
        normal = perp(edge)     # notably not normalized, for speed; can't be used for Circle SAT
        xrange = minmaxV(v[1],w[1])
        yrange = minmaxV(v[2],w[2])
        ndotv = dot(v, normal)
        new(v, w, edge, normal, xrange, yrange, ndotv)
    end
end
Line{T}(v::Vector2{T}, w::Vector2{T}) = Line{T}(v,w)   # so Line(v,w) works without {T}
projectNextrema{T}(L::Line{T}, n::Vector2{T}) = minmaxV(dot(L.v,n), dot(L.w,n))

typealias PolyOrLine{T} Union(Polygon{T}, Line{T})

type Compound2D{T} <: Shape2D{T}
    parts::Vector{Shape2D{T}}
    xrange::Vector2{T}
    yrange::Vector2{T}
end
function Compound2D{T}(parts::Vector{Shape2D{T}})
    xrange = Vector2{T}(extrema(vcat([S.xrange for S in parts]...))...)
    yrange = Vector2{T}(extrema(vcat([S.yrange for S in parts]...))...)
    Compound2D{T}(parts, xrange, yrange)
end

# type AABB{T} <: Shape2D{T}  # a bit late to the party; could refactor other Shapes
#     xrange::Vector2{T}
#     yrange::Vector2{T}
# end
# AABB{T}(x1::T, x2::T, y1::T, y2::T) = AABB{T}(Vector2(x1,x2), Vector2(y1,y2))

# ---------- Separating Axis ----------

is_separating_axis(S1::Shape2D, S2::Shape2D, ax::Vector2) = !overlapping(projectNextrema(S1, ax),
                                                                         projectNextrema(S2, ax))
is_separating_axis(P::Polygon, S::Shape2D, i::Integer) = !overlapping(P.nextrema[i],
                                                                      projectNextrema(S, P.normals[i]))
is_separating_axis(L::Line, P::Polygon) = !inrange(L.ndotv, projectNextrema(P, L.normal))

# ---------- Collision Checking ----------

## Broadphase
AABBseparated(S1::Shape2D, S2::Shape2D) = !(overlapping(S1.xrange, S2.xrange) && overlapping(S1.yrange, S2.yrange))
pointinAABB(p::Vector2, S::Shape2D) = inrange(p[1], S.xrange) && inrange(p[2], S.yrange)
## Point collision
colliding(p::Vector2, C::Circle) = norm2(p - C.c) <= C.r^2
colliding(C::Circle, p::Vector2) = norm2(p - C.c) <= C.r^2
function colliding(p::Vector2, P::Polygon)
    !pointinAABB(p,P) && return false
    for i in 1:length(P.normals)
        !inrange(dot(p, P.normals[i]), P.nextrema[i]) && return false
    end
    true
end
colliding(P::Polygon, p::Vector2) = colliding(p,P)
function colliding(p::Vector2, C::Compound2D)
    !pointinAABB(p,C) && return false
    for P in C.parts
        colliding(p,P) && return true
    end
    false
end
colliding(C::Compound2D, p::Vector2) = colliding(p,C)
## Shape collision
colliding(C1::Circle, C2::Circle) = norm2(C1.c - C2.c) <= (C1.r + C2.r)^2
function colliding(C::Circle, P::Polygon)
    AABBseparated(C,P) && return false
    N = length(P.edges)
    for i in 1:N
        pi2c = C.c - P.points[i]
        VR = voronoi_region(pi2c, P.edges[i])
        if VR == 0
            dot(pi2c, P.normals[i]) > C.r && return false
        else
            j = wrap1(i+VR,N)
            pj2c = C.c - P.points[j]
            VR*voronoi_region(pj2c, P.edges[j]) < 0 && norm2(VR > 0 ? pj2c : pi2c) > C.r^2 && return false
        end
    end
    true
end
colliding(P::Polygon, C::Circle) = colliding(C,P)
function colliding(P1::Polygon, P2::Polygon)
    AABBseparated(P1,P2) && return false
    for i in 1:length(P1.normals)
        is_separating_axis(P1,P2,i) && return false
    end
    for i in 1:length(P2.normals)
        is_separating_axis(P2,P1,i) && return false
    end
    true
end
function colliding(C::Compound2D, S::Shape2D)
    AABBseparated(C,S) && return false
    for P in C.parts
        colliding(P,S) && return true
    end
    false
end
colliding(S::Circle, C::Compound2D) = colliding(C,S)   # colliding(S::Shape2D, C::Compound2D) is type-ambiguous?
colliding(S::Polygon, C::Compound2D) = colliding(C,S)
colliding(S::Line, C::Compound2D) = colliding(C,S)
## Swept collisions
function colliding_ends_free(L::Line, C::Circle)
    AABBseparated(L,C) && return false
    cv = L.v-C.c
    cw = L.w-C.c
    d2 = norm2(L.edge)
    d2*C.r^2 < cross(cv, cw)^2 && return false
    0 <= dot(cv,L.edge) <= d2
end
colliding(L::Line, C::Circle) = colliding_ends_free(L,C) || colliding(L.v,C) || colliding(L.w,C)
colliding_ends_free(C::Circle, L::Line) = colliding_ends_free(L,C) 
colliding(C::Circle, L::Line) = colliding(L,C)
function colliding_ends_free(L::Line, P::Polygon)
    AABBseparated(L,P) && return false
    is_separating_axis(L,P) && return false
    for i in 1:length(P.normals)
        is_separating_axis(P,L,i) && return false
    end
    true
end
colliding(L::Line, P::Polygon) = colliding_ends_free(L, P)
colliding_ends_free(P::Polygon, L::Line) = colliding_ends_free(L,P) 
colliding(P::Polygon, L::Line) = colliding(L,P)

# ---------- Plotting ----------

plot(C::Circle; kwargs...) = plot_circle(C.c, C.r; kwargs...)
plot(P::Polygon; kwargs...) = plot_polygon(P.points; kwargs...)
plot(C::Compound2D; kwargs...) = for P in C.parts; plot(P; kwargs...); end
plot(L::Line; kwargs...) = plot_path([L.v L.w]; kwargs...)

# function composable(P::Polygon, opts = (fill("red"), fillopacity(1.)))
#     compose(context(), polygon([(p[1],p[2]) for p in P.points]), opts...)
# end

# function composable(C::Compound2D, opts = (fill("red"), fillopacity(1.)))
#     compose(context(), [composable(P, opts) for P in C.parts]...)
# end

# function composable(L::Line, opts = (stroke("black"),))
#     compose(context(), line([(L.v[1], L.v[2]), (L.w[1], L.w[2])]), opts...)
# end