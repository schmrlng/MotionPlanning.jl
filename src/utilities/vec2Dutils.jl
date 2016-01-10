import Base: atan2, full, cross
import FixedSizeArrays: unit
export norm2, perp, project, projectN, reflect, reflectN, rotate
export projectNextrema, voronoi_region, overlapping, inrange, minmaxV, wrap1

full{T}(v::Vec{2,T}) = convert(Vector{T}, v)
unit(v::Vec{2}) = v/norm(v)
norm2(v::Vec{2}) = dot(v,v)
perp(v::Vec{2}) = Vec(v[2],-v[1])
cross(v1::Vec{2}, v2::Vec{2}) = v1[1]*v2[2]-v1[2]*v2[1]
project(v1::Vec{2}, v2::Vec{2}) = (dot(v1,v2) / norm2(v2)) * v2
projectN(v1::Vec{2}, v2::Vec{2}) = dot(v1,v2) * v2
reflect(v1::Vec{2}, ax::Vec{2}) = v1 - 2*project(v1,ax)
reflectN(v1::Vec{2}, ax::Vec{2}) = v1 - 2*projectN(v1,ax)
function rotate(v::Vec{2}, theta)
    c = cos(theta)
    s = sin(theta)
    Vec(v[1]*c - v[2]*s, v[1]*s + v[2]*c)
end
atan2(v::Vec{2}) = atan2(v[2], v[1])
function projectNextrema{T<:AbstractFloat}(pts::AbstractVector{Vec{2,T}}, n::Vec{2,T})
    dmin = T(Inf)
    dmax = -T(Inf)
    for p in pts
        d = dot(p,n)
        d < dmin && (dmin = d)
        d > dmax && (dmax = d)
    end
    Vec(dmin, dmax)
end
function voronoi_region(p::Vec{2}, v::Vec{2})
    d = dot(p,v)
    d < 0 && return -1
    d > norm2(v) && return 1
    return 0
end
overlapping(I1::Vec{2}, I2::Vec{2}) = I1[1] <= I2[2] && I2[1] <= I1[2]
inrange(x, I::Vec{2}) = I[1] <= x <= I[2]
minmaxV(x, y) = x < y ? Vec(x,y) : Vec(y,x)  # Base.minmax returns tuple 
wrap1(i, N) = i < 1 ? N : i > N ? 1 : i      # marginally faster than mod1