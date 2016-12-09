import Base: atan2, cross, dot, normalize
export norm2, perp, project, projectN, reflect, reflectN, rotate
export projectNextrema, voronoi_region, overlapping, ininterval, minmaxV, wrap1

norm2(v::SVector{2}) = dot(v,v)
perp(v::SVector{2}) = SVector(v[2],-v[1])
cross(v1::SVector{2}, v2::SVector{2}) = v1[1]*v2[2] - v1[2]*v2[1]
project(v1::SVector{2}, v2::SVector{2}) = (dot(v1,v2) / norm2(v2)) * v2
projectN(v1::SVector{2}, v2::SVector{2}) = dot(v1,v2) * v2
reflect(v1::SVector{2}, ax::SVector{2}) = v1 - 2*project(v1,ax)
reflectN(v1::SVector{2}, ax::SVector{2}) = v1 - 2*projectN(v1,ax)
function rotate(v::SVector{2}, θ)
    c = cos(θ)
    s = sin(θ)
    SVector(v[1]*c - v[2]*s, v[1]*s + v[2]*c)
end
atan2(v::SVector{2}) = atan2(v[2], v[1])
function projectNextrema{T}(pts::AbstractVector{SVector{2,T}}, n::SVector{2,T})
    dmin = T(Inf)
    dmax = -T(Inf)
    for p in pts
        d = dot(p,n)
        d < dmin && (dmin = d)
        d > dmax && (dmax = d)
    end
    SVector(dmin, dmax)
end
function voronoi_region(p::SVector{2}, v::SVector{2})
    d = dot(p,v)
    d < 0 && return -1
    d > norm2(v) && return 1
    return 0
end
overlapping(I1::SVector{2}, I2::SVector{2}) = I1[1] <= I2[2] && I2[1] <= I1[2]
ininterval(x, I::SVector{2}) = I[1] <= x <= I[2]
minmaxV(x, y) = x < y ? SVector(x,y) : SVector(y,x)  # Base.minmax returns tuple 
wrap1(i, N) = i < 1 ? N : i > N ? 1 : i      # marginally faster than mod1