export SimpleConvexSet, AxisAlignedBox, Ball, SetComplement
export inflate, boundingbox

abstract type SimpleConvexSet end

struct AxisAlignedBox{N,T} <: SimpleConvexSet
    lo::SVector{N,T}
    hi::SVector{N,T}
end
AxisAlignedBox(lo, hi) = AxisAlignedBox{length(lo),eltype(lo)}(lo, hi)
Base.in(x, B::AxisAlignedBox{N}) where {N} = all(B.lo[i] <= x[i] <= B.hi[i] for i in 1:N)
Base.rand(B::AxisAlignedBox) = B.lo + (B.hi - B.lo).*rand(typeof(B.lo))
inflate(B::AxisAlignedBox, ε) = AxisAlignedBox(B.lo .- ε, B.hi .+ ε)
boundingbox(B::AxisAlignedBox) = B
@recipe function f(B::AxisAlignedBox; dims=(1, 2))
    seriestype :=  :shape
    fillcolor  --> :match
    linecolor  --> :match
    label      --> ""
    x, y = dims
    SVector(B.lo[x], B.hi[x], B.hi[x], B.lo[x]), SVector(B.lo[y], B.lo[y], B.hi[y], B.hi[y])
end

struct Ball{N,T} <: SimpleConvexSet
    c::SVector{N,T}
    r::T
end
Ball(c, r) = Ball{length(c),eltype(c)}(c, r)
Base.in(x, B::Ball{N}) where {N} = norm(x - B.c) <= B.r
Base.rand(B::Ball) = (x = B.c + 2*B.r*(rand(typeof(B.c)) .- 1//2); x in B ? x : rand(B))
inflate(B::Ball, ε) = Ball(B.c, B.r + ε)
boundingbox(B::Ball) = AxisAlignedBox(B.c .- B.r, B.c .+ B.r)
@recipe function f(B::Ball; dims=(1, 2), n=64)
    seriestype :=  :shape
    fillcolor  --> :match
    linecolor  --> :match
    label      --> ""
    x, y = dims
    θ = range(-π, stop=π, length=n)
    B.c[x] .+ B.r*cos.(θ), B.c[y] .+ B.r*sin.(θ)
end

struct SetComplement{S}
    set::S
end
SetComplement(S::SetComplement) = S.set
Base.in(x, S::SetComplement) = !(x in S.set)
@recipe f(S::SetComplement; dims=(1,2)) = nothing
