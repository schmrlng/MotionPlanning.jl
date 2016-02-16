import Base.eltype
export Shape2D, Circle, Polygon, Box2D, Line, Compound2D
export colliding, colliding_ends_free, closest, closeR

# ---------- Shape Definitions ----------

abstract Shape2D{T<:AbstractFloat}
typealias AnyVec2{T} Union{AbstractVector{T}, Vec{2,T}, Tuple{T,T}}
eltype{T<:AbstractFloat}(::Type{Shape2D{T}}) = T
eltype{T<:Shape2D}(::Type{T}) = eltype(super(T))

immutable Circle{T} <: Shape2D{T}
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
changeprecision{T<:AbstractFloat}(::Type{T}, C::Circle) = Circle(changeprecision(T, C.c), T(C.r))

immutable Polygon{T} <: Shape2D{T}
    points::Vector{Vec{2,T}}
    edges::Vector{Vec{2,T}}
    normals::Vector{Vec{2,T}}
    xrange::Vec{2,T}
    yrange::Vec{2,T}
    nextrema::Vector{Vec{2,T}}

    function Polygon(points::Vector{Vec{2,T}})
        N = length(points)
        N < 3 && error("Polygons need at least 3 points! Try Line?")
        if sum([(points[wrap1(i+1,N)][1] - points[i][1])*(points[wrap1(i+1,N)][2] + points[i][2]) for i in 1:N]) > 0
            reverse!(points)
        end
        edges = [points[wrap1(i+1,N)] - points[i] for i in 1:N]
        normals = [normalize(perp(g)) for g in edges]
        any(-pi .<= diff([T[atan2(n) for n in normals]; atan2(normals[1])]) .<= 0) && error("Polygon must be convex")
        xrange = Vec(extrema([points[i][1] for i in 1:N]))
        yrange = Vec(extrema([points[i][2] for i in 1:N]))
        nextrema = [projectNextrema(points, normals[i]) for i in 1:N]
        new(points, edges, normals, xrange, yrange, nextrema)
    end
end
Polygon{T<:AnyVec2}(points::Vector{T}) = Polygon{eltype(T)}(Vec{2,eltype(T)}[Vec(p[1],p[2]) for p in points])
Box2D(xr::AnyVec2, yr::AnyVec2) = Polygon([Vec(xr[1], yr[1]),
                                           Vec(xr[2], yr[1]),
                                           Vec(xr[2], yr[2]),
                                           Vec(xr[1], yr[2])])
projectNextrema(P::Polygon, n::Vec{2}) = projectNextrema(P.points, n)
changeprecision{T<:AbstractFloat}(::Type{T}, P::Polygon) = Polygon(changeprecision(T, P.points))

immutable Line{T} <: Shape2D{T}   # special case of Polygon; only 1 edge/normal
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
Line(v::AnyVec2, w::AnyVec2) = Line{eltype(v)}(Vec(v[1],v[2]), Vec(w[1],w[2]))
projectNextrema(L::Line, n::Vec{2}) = minmaxV(dot(L.v,n), dot(L.w,n))
changeprecision{T<:AbstractFloat}(::Type{T}, L::Line) = Line(changeprecision(T, L.v), changeprecision(T, L.w))

immutable Compound2D{T} <: Shape2D{T}
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
changeprecision{T<:AbstractFloat}(::Type{T}, C::Compound2D) =
    Compound2D(Shape2D{T}[changeprecision(T,P) for P in C.parts])

typealias PolyOrLine{T} Union{Polygon{T}, Line{T}}
typealias Basic2D{T} Union{Circle{T}, Polygon{T}}    # TODO: Figure out why I have Shape2D, PolyOrLine, and Basic2D

# type AABB{T} <: Shape2D{T}  # a bit late to the party; could refactor other Shapes
#     xrange::Vector2{T}
#     yrange::Vector2{T}
# end
# AABB{T}(x1::T, x2::T, y1::T, y2::T) = AABB{T}(Vector2(x1,x2), Vector2(y1,y2))

# ---------- Separating Axis ----------

is_separating_axis(S1::Shape2D, S2::Shape2D, ax::Vec{2}) = !overlapping(projectNextrema(S1,ax), projectNextrema(S2,ax))
is_separating_axis(P::Polygon, S::Shape2D, i::Integer) = !overlapping(P.nextrema[i], projectNextrema(S, P.normals[i]))
is_separating_axis(L::Line, P::Polygon) = !ininterval(L.ndotv, projectNextrema(P, L.normal))

# ---------- Collision Checking ----------

## Broadphase
AABBseparated(S1::Shape2D, S2::Shape2D) = !(overlapping(S1.xrange, S2.xrange) && overlapping(S1.yrange, S2.yrange))
pointinAABB(p::Vec{2}, S::Shape2D) = ininterval(p[1], S.xrange) && ininterval(p[2], S.yrange)
## Point collision
colliding(p::Vec{2}, C::Circle) = norm2(p - C.c) <= C.r^2
colliding(C::Circle, p::Vec{2}) = norm2(p - C.c) <= C.r^2
function colliding(p::Vec{2}, P::Polygon)
    !pointinAABB(p,P) && return false
    @all [!ininterval(dot(p, P.normals[i]), P.nextrema[i]) for i in 1:length(P.normals)]
end
colliding(P::Polygon, p::Vec{2}) = colliding(p,P)
function colliding(p::Vec{2}, C::Compound2D)
    !pointinAABB(p,C) && return false
    @any [colliding(p,P) for P in C.parts]
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
    !(@any [is_separating_axis(P1,P2,i) for i in 1:length(P1.normals)]) &&
    !(@any [is_separating_axis(P2,P1,i) for i in 1:length(P2.normals)])
end
function colliding(C::Compound2D, S::Shape2D)
    AABBseparated(C,S) && return false
    @any [colliding(P,S) for P in C.parts]
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
    !(@any [is_separating_axis(P,L,i) for i in 1:length(P.normals)])
end
colliding_ends_free(B::Basic2D, L::Line) = colliding_ends_free(L,B) 
colliding(L::Line, B::Basic2D) = colliding_ends_free(L,B) || colliding(L.v,B) || colliding(L.w,B)
colliding(B::Basic2D, L::Line) = colliding(L,B)
colliding(L::Line, C::Compound2D) = colliding(C,L)
function colliding_ends_free(L::Line, C::Compound2D)
    AABBseparated(L,C) && return false
    @any [colliding_ends_free(L,P) for P in C.parts]
end
colliding_ends_free(C::Compound2D, L::Line) = colliding(C,L)

# ---------- Transformations ----------

inflate{T}(C::Circle{T}, eps; roundcorners = true) = Circle(C.c, C.r + T(eps))
function inflate{T}(P::Polygon{T}, eps; roundcorners = true)
    push_out_corner_vector(n0, n1) = (abs(cross(n0, n1)) < T(1e-6) ? n0 : (perp(n1) - perp(n0)) / cross(n0, n1))
    eps = T(eps)
    N = length(P.points)
    S = eltype(P.points)
    if !roundcorners
        return Polygon(S[P.points[i] + eps*push_out_corner_vector(P.normals[wrap1(i-1,N)], P.normals[i]) for i in 1:N])
    end
    Compound2D(
        Polygon(vcat([S[P.points[i] + eps*P.normals[wrap1(i-1,N)], P.points[i] + eps*P.normals[i]] for i in 1:N]...)),
        [Circle(p, eps) for p in P.points]...
    )
end
inflate{T}(C::Compound2D{T}, eps; roundcorners = true) =
    Compound2D(Shape2D{T}[inflate(P, T(eps); roundcorners = roundcorners) for P in C.parts])

# ---------- Closest Point ---------

@unfix function closest(p::Vec{2}, C::Circle)
    xmin = C.c + C.r*normalize(p - C.c)
    norm2(p-xmin), xmin
end
@unfix closest(p::Vec{2}, C::Circle, W::Mat{2,2}) = closest(p, C, eigfact(dense(W)))
@unfix function closest{T}(p::Vec{2,T}, C::Circle{T}, EF::Base.Eigen)
    ctop = p - C.c
    v1 = Vec(EF[:vectors][1:2,1])
    v2 = Vec(EF[:vectors][1:2,2])
    p1 = dot(v1, ctop)
    p2 = dot(v2, ctop)
    s1, s2 = EF[:values]
    lambda = one(T)
    f = (p1*s1/(lambda+s1))^2 + (p2*s2/(lambda+s2))^2 - C.r^2
    while abs(f) > 1e-8
        fp = -2/(lambda+s1)*(p1*s1/(lambda+s1))^2 + -2/(lambda+s2)*(p2*s2/(lambda+s2))^2
        alpha = one(T)
        lambdanew = one(T)
        fnew = one(T)
        while true      # crappy linesearch
            lambdanew = lambda - alpha*f/fp
            fnew = (p1*s1/(lambdanew+s1))^2 + (p2*s2/(lambdanew+s2))^2 - C.r^2
            abs(fnew) < abs(f) && break
            alpha /= 2
        end
        f = fnew
        lambda = lambdanew
    end
    xmin = C.c + v1*p1*s1/(lambda+s1) + v2*p2*s2/(lambda+s2)
    s1*(p1-p1*s1/(lambda+s1))^2 + s2*(p2-p2*s2/(lambda+s2))^2, xmin
end
@unfix closest(p::Vec{2}, P::Polygon) = closest_polypts(p, P.points)
@unfix function closest_polypts{T}(p::Vec{2,T}, points::Vector{Vec{2,T}})
    N = length(points)
    d2min, vmin = T(Inf), points[1]
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
@unfix function closest(p::Vec{2}, P::Polygon, W::Mat{2,2})
    L = Mat(dense(chol(dense(W))))
    xmin = inv(L)*closest_polypts(L*p, eltype(P.points)[L*pt for pt in P.points])[2]
    dot(xmin-p, W*(xmin-p)), xmin
end
@unfix function closest{T}(p::Vec{2,T}, C::Compound2D{T})
    d2min, vmin = T(Inf), p
    for P in C.parts
        d2, v = closest(p, P)
        if d2 < d2min
            d2min, vmin = d2, v
        end
    end
    d2min, vmin
end
@unfix function closest{T}(p::Vec{2,T}, C::Compound2D{T}, W::Mat{2,2,T})  # TODO: macro, iterator, or v0.5 generator
    d2min, vmin = T(Inf), p
    for P in C.parts
        d2, v = closest(p, P, W)
        if d2 < d2min
            d2min, vmin = d2, v
        end
    end
    d2min, vmin
end

@unfix function closeR(p::Vec{2}, P::Basic2D, W::Mat{2,2}, r2)  # pretty janky
    cplist = [closest(p, P, W)]
    cplist[1][1] < r2 ? cplist[1:1] : cplist[1:0]
end
@unfix closeR(p::Vec{2}, C::Compound2D, W::Mat{2,2}, r2) =
    sort!(vcat([closeR(p, P, W, r2) for P in C.parts]...), by=first)

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