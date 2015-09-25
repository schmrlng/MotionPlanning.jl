import Base.atan2, Base.full, Base.inv
export norm2, perp, project, projectN, reflect, reflectN, rotate
export projectNextrema, voronoi_region, overlapping, inrange, minmaxV, wrap1

full(v::Vector2) = convert(Vector{eltype(v)}, v)
function inv{T}(M::Matrix2x2{T})
    d = M[1,1]*M[2,2] - M[1,2]*M[2,1]
    Matrix2x2(Vector2(M[2,2]/d, -M[2,1]/d), Vector2(-M[1,2]/d, M[1,1]/d))
end
norm2{T<:AbstractFloat}(v::Vector2{T}) = dot(v,v)
perp{T<:AbstractFloat}(v::Vector2{T}) = Vector2{T}(v[2],-v[1])
project{T<:AbstractFloat}(v1::Vector2{T}, v2::Vector2{T}) = (dot(v1,v2) / norm2(v2)) * v2
projectN{T<:AbstractFloat}(v1::Vector2{T}, v2::Vector2{T}) = dot(v1,v2) * v2
reflect{T<:AbstractFloat}(v1::Vector2{T}, ax::Vector2{T}) = v1 - 2*project(v1,ax)
reflectN{T<:AbstractFloat}(v1::Vector2{T}, ax::Vector2{T}) = v1 - 2*projectN(v1,ax)
function rotate{T<:AbstractFloat}(v::Vector2{T}, theta::T)
    c = cos(theta)
    s = sin(theta)
    Vector2{T}(v[1]*c - v[2]*s, v[1]*s + v[2]*c)
end
atan2{T<:AbstractFloat}(v::Vector2{T}) = atan2(v[2], v[1])
function projectNextrema{T<:AbstractFloat}(pts::AbstractVector{Vector2{T}}, n::Vector2{T})
    dmin = T(Inf)
    dmax = -T(Inf)
    for p in pts
        d = dot(p,n)
        d < dmin && (dmin = d)
        d > dmax && (dmax = d)
    end
    Vector2{T}(dmin, dmax)
end
function voronoi_region{T<:AbstractFloat}(p::Vector2{T}, v::Vector2{T})
    d = dot(p,v)
    d < 0 && return -1
    d > norm2(v) && return 1
    return 0
end
overlapping{T<:AbstractFloat}(I1::Vector2{T}, I2::Vector2{T}) = I1[1] <= I2[2] && I2[1] <= I1[2]
inrange{T<:AbstractFloat}(x::T, I::Vector2{T}) = I[1] <= x <= I[2]
minmaxV{T<:AbstractFloat}(x::T, y::T) = x < y ? Vector2{T}(x,y) : Vector2{T}(y,x)  # Base.minmax returns tuple 
wrap1(i,N) = i < 1 ? N : i > N ? 1 : i      # marginally faster than mod1