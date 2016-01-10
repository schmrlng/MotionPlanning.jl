import Base.eltype
export Shape2D, Circle, Polygon, Box2D, Line, Compound2D, colliding, colliding_ends_free, closest, close

# ---------- Shape Definitions ----------

abstract Shape2D{T<:AbstractFloat}
typealias AnyVec2{T} Union{AbstractVector{T}, Vec{2,T}, Tuple{T,T}}
eltype{T<:AbstractFloat}(::Type{Shape2D{T}}) = T
eltype{T<:Shape2D}(::Type{T}) = eltype(super(T))

type Circle{T} <: Shape2D{T}
    c::Vec{2,T}
    r::T
    xrange::Vec{2,T}
    yrange::Vec{2,T}

    function Circle(c::Vec{2,T}, r::T, xrange::Vec{2,T}, yrange::Vec{2,T})
        r <= 0 && error("Radius must be positive")
        new(c, r, xrange, yrange)
    end
end
Circle(c::AnyVec2, r) = Circle{eltype(c)}(Vec(c[1], c[2]), r,
                                          Vec(c[1]-r, c[1]+r),
                                          Vec(c[2]-r, c[2]+r))
projectNextrema(C::Circle, n::Vec{2}) = (d = dot(C.c,n); Vec(d-C.r, d+C.r))

type Polygon{T} <: Shape2D{T}
    points::Vector{Vec{2,T}}
    edges::Vector{Vec{2,T}}
    normals::Vector{Vec{2,T}}
    xrange::Vec{2,T}
    yrange::Vec{2,T}
    nextrema::Vector{Vec{2,T}}

    function Polygon(points::Vector{Vec{2,T}})
        N = length(points)
        N < 3 && error("Polygons need at least 3 points! Try Line?")
        sum([(points[wrap1(i+1,N)][1] - points[i][1])*(points[wrap1(i+1,N)][2] + points[i][2]) for i in 1:N]) > 0 && reverse!(points)
        edges = [points[wrap1(i+1,N)] - points[i] for i in 1:N]
        normals = [unit(perp(g)) for g in edges]
        any(-pi .<= diff([T[atan2(n) for n in normals]; atan2(normals[1])]) .<= 0) && error("Polygon needs to be convex")
        xrange = Vec(extrema([points[i][1] for i in 1:N]))
        yrange = Vec(extrema([points[i][2] for i in 1:N]))
        nextrema = [projectNextrema(points, normals[i]) for i in 1:N]
        new(points, edges, normals, xrange, yrange, nextrema)
    end
end
Polygon{T<:AnyVec2}(points::Vector{T}) = Polygon{eltype(T)}(Vec{2,eltype(T)}[Vec(p) for p in points])
Box2D(xr::AnyVec2, yr::AnyVec2) = Polygon([Vec(xr[1], yr[1]),
                                           Vec(xr[2], yr[1]),
                                           Vec(xr[2], yr[2]),
                                           Vec(xr[1], yr[2])])
projectNextrema(P::Polygon, n::Vec{2}) = projectNextrema(P.points, n)

type Line{T} <: Shape2D{T}   # special case of Polygon; only 1 edge/normal
    v::Vec{2,T}
    w::Vec{2,T}
    edge::Vec{2,T}
    normal::Vec{2,T}
    xrange::Vec{2,T}
    yrange::Vec{2,T}
    ndotv::T

    function Line(v::Vec{2,T}, w::Vec{2,T})
        edge = w - v
        normal = perp(edge)     # notably not normalized, for speed; can't be used for Circle SAT
        xrange = minmaxV(v[1],w[1])
        yrange = minmaxV(v[2],w[2])
        ndotv = dot(v, normal)
        new(v, w, edge, normal, xrange, yrange, ndotv)
    end
end
Line(v::AnyVec2, w::AnyVec2) = Line{eltype(v)}(Vec(v), Vec(w))
projectNextrema(L::Line, n::Vec{2}) = minmaxV(dot(L.v,n), dot(L.w,n))

type Compound2D{T} <: Shape2D{T}
    parts::Vector{Shape2D{T}}
    xrange::Vec{2,T}
    yrange::Vec{2,T}
end
function Compound2D{S<:Shape2D}(parts::Vector{S})
    T = eltype(S)
    isempty(parts) && return Compound2D{T}(Shape2D{T}[], zero(Vec{2,T}), zero(Vec{2,T}))
    Compound2D(parts...)
end
function Compound2D{T}(parts::Shape2D{T}...)
    xrange = Vec(minimum([P.xrange[1] for P in parts]), maximum([P.xrange[2] for P in parts]))
    yrange = Vec(minimum([P.yrange[1] for P in parts]), maximum([P.yrange[2] for P in parts]))
    Compound2D{T}(Shape2D{T}[parts...], xrange, yrange)
end

typealias PolyOrLine{T} Union{Polygon{T}, Line{T}}
typealias Basic2D{T} Union{Circle{T}, Polygon{T}}    # TODO: Figure out why I have Shape2D, PolyOrLine, and Basic2D

# type AABB{T} <: Shape2D{T}  # a bit late to the party; could refactor other Shapes
#     xrange::Vector2{T}
#     yrange::Vector2{T}
# end
# AABB{T}(x1::T, x2::T, y1::T, y2::T) = AABB{T}(Vector2(x1,x2), Vector2(y1,y2))

# ---------- Separating Axis ----------

is_separating_axis(S1::Shape2D, S2::Shape2D, ax::Vec{2}) = !overlapping(projectNextrema(S1, ax),
                                                                         projectNextrema(S2, ax))
is_separating_axis(P::Polygon, S::Shape2D, i::Integer) = !overlapping(P.nextrema[i],
                                                                      projectNextrema(S, P.normals[i]))
is_separating_axis(L::Line, P::Polygon) = !inrange(L.ndotv, projectNextrema(P, L.normal))

# ---------- Collision Checking ----------

## Broadphase
AABBseparated(S1::Shape2D, S2::Shape2D) = !(overlapping(S1.xrange, S2.xrange) && overlapping(S1.yrange, S2.yrange))
pointinAABB(p::Vec{2}, S::Shape2D) = inrange(p[1], S.xrange) && inrange(p[2], S.yrange)
## Point collision
colliding(p::Vec{2}, C::Circle) = norm2(p - C.c) <= C.r^2
colliding(C::Circle, p::Vec{2}) = norm2(p - C.c) <= C.r^2
function colliding(p::Vec{2}, P::Polygon)
    !pointinAABB(p,P) && return false
    for i in 1:length(P.normals)
        !inrange(dot(p, P.normals[i]), P.nextrema[i]) && return false
    end
    true
end
colliding(P::Polygon, p::Vec{2}) = colliding(p,P)
function colliding(p::Vec{2}, C::Compound2D)
    !pointinAABB(p,C) && return false
    for P in C.parts
        colliding(p,P) && return true
    end
    false
end
colliding(C::Compound2D, p::Vec{2}) = colliding(p,C)
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
## Swept collisions
function colliding_ends_free(L::Line, C::Circle)
    AABBseparated(L,C) && return false
    vc = C.c-L.v
    d2 = norm2(L.edge)
    d2*C.r^2 < cross(L.edge, vc)^2 && return false
    0 <= dot(vc, L.edge) <= d2
end
function colliding_ends_free(L::Line, P::Polygon)
    AABBseparated(L,P) && return false
    is_separating_axis(L,P) && return false
    for i in 1:length(P.normals)
        is_separating_axis(P,L,i) && return false
    end
    true
end
colliding_ends_free(B::Basic2D, L::Line) = colliding_ends_free(L,B) 
colliding(L::Line, B::Basic2D) = colliding_ends_free(L,B) || colliding(L.v,B) || colliding(L.w,B)
colliding(B::Basic2D, L::Line) = colliding(L,B)
colliding(S::Line, C::Compound2D) = colliding(C,S)
colliding_ends_free(S::Line, C::Compound2D) = colliding(C,S)
colliding_ends_free(C::Compound2D, S::Line) = colliding(C,S)

# ---------- Transformations ----------

inflate(C::Circle, eps; roundcorners = true) = Circle(C.c, C.r + eps)
function inflate(P::Polygon, eps; roundcorners = true)
    push_out_corner_vector(n0, n1) = (abs(cross(n0, n1)) < 1e-6 ? n0 : (perp(n1) - perp(n0)) / cross(n0, n1))
    N = length(P.points)
    S = eltype(P.points)
    !roundcorners && return Polygon(S[P.points[i] + eps*push_out_corner_vector(P.normals[wrap1(i-1,N)], P.normals[i]) for i in 1:N])
    Compound2D(Polygon(vcat([S[P.points[i] + eps*P.normals[wrap1(i-1,N)], P.points[i] + eps*P.normals[i]] for i in 1:N]...)), [Circle(p, eps) for p in P.points]...)
end
inflate{T}(C::Compound2D{T}, eps; roundcorners = true) = Compound2D(Shape2D{T}[inflate(P, eps; roundcorners = roundcorners) for P in C.parts])

# ---------- Closest Point ---------

function closest(p::Vec{2}, C::Circle)
    xmin = C.c + C.r*unit(p - C.c)
    norm2(p-xmin), xmin
end
closest(p::Vec{2}, C::Circle, W::AbstractMatrix) = closest(p, C, eigfact(full(W)))
function closest(p::Vec{2}, C::Circle, EF::Base.Eigen)
    ctop = p - C.c
    v1 = Vec(EF.vectors[1:2,1])
    v2 = Vec(EF.vectors[1:2,2])
    p1 = dot(v1, ctop)
    p2 = dot(v2, ctop)
    s1, s2 = EF.values
    lambda = 1.
    f = (p1*s1/(lambda+s1))^2 + (p2*s2/(lambda+s2))^2 - C.r^2
    while abs(f) > 1e-8
        fp = -2/(lambda+s1)*(p1*s1/(lambda+s1))^2 + -2/(lambda+s2)*(p2*s2/(lambda+s2))^2
        alpha = 1.
        lambdanew = 1.
        fnew = 1.
        while true      # crappy linesearch
            lambdanew = lambda - alpha*f/fp
            fnew = (p1*s1/(lambdanew+s1))^2 + (p2*s2/(lambdanew+s2))^2 - C.r^2
            abs(fnew) < abs(f) && break
            alpha /= 2.
        end
        f = fnew
        lambda = lambdanew
    end
    xmin = C.c + v1*p1*s1/(lambda+s1) + v2*p2*s2/(lambda+s2)
    s1*(p1-p1*s1/(lambda+s1))^2 + s2*(p2-p2*s2/(lambda+s2))^2, xmin
end
closest(p::Vec{2}, P::Polygon) = closest(p, P.points)
function closest(p::Vec{2}, points::Vector)
    N = length(points)
    d2min, vmin = Inf, points[1]
    for i in 1:N
        edge = points[wrap1(i+1,N)] - points[i]
        x = dot(edge, p - points[i]) / norm2(edge)
        v = (x < 0 ? points[i] : (x < 1 ? points[i] + x*edge : points[wrap1(i+1,N)]))
        d2 = norm2(p - v)
        if d2 < d2min
            d2min, vmin = d2, v
        end
    end
    d2min, vmin
end
function closest(p::Vec{2}, P::Polygon, W::AbstractMatrix)
    L = Matrix2x2(chol(W))
    xmin = inv(L)*closest(L*p, eltype(P.points)[L*pt for pt in P.points])[2]
    dot(xmin-p, Matrix2x2(W)*(xmin-p)), xmin
end
function closest(p::Vec{2}, C::Compound2D)
    d2min, vmin = Inf, p
    for P in C.parts
        d2, v = closest(p, P)
        if d2 < d2min
            d2min, vmin = d2, v
        end
    end
    d2min, vmin
end
function closest(p::Vec{2}, C::Compound2D, W::AbstractMatrix)  # refactor/combine with above... someday
    d2min, vmin = Inf, p
    for P in C.parts
        d2, v = closest(p, P, W)
        if d2 < d2min
            d2min, vmin = d2, v
        end
    end
    d2min, vmin
end

function close(p::Vec{2}, P::Basic2D, W::AbstractMatrix, r2)  # pretty janky
    cplist = [closest(p, P, W)]
    cplist[1][1] < r2 ? cplist[1:1] : cplist[1:0]
end
close(p::Vec{2}, C::Compound2D, W::AbstractMatrix, r2) = sort!(vcat([close(p, P, W, r2) for P in C.parts]...), by=first)

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