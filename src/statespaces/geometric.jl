export RealVectorMetricSpace, BoundedEuclideanStateSpace, UnitHypercube

immutable RealVectorMetricSpace{T<:FloatingPoint} <: RealVectorStateSpace
    dim::Int
    lo::Vector{T}
    hi::Vector{T}
    dist::Metric
end

vector_to_state{T}(v::AbstractVector{T}, SS::RealVectorMetricSpace) = SS.dim == 2 ? convert(Vector2{T}, v) : v
sample_space(SS::RealVectorMetricSpace) = vector_to_state(SS.lo + rand(SS.dim).*(SS.hi-SS.lo), SS)   # TODO: @devec
function volume(SS::RealVectorMetricSpace)
    SS.dist == Euclidean() && return prod(SS.hi-SS.lo)
    error("Volume not yet implemented for non-Euclidean metrics!")
end
function defaultNN{T}(SS::RealVectorMetricSpace{T}, init)
    V = typeof(init)[init]
    SS.dist == Euclidean() && return EuclideanNN_KDTree(V)
    MetricNN_BruteForce(V, SS.dist)
end

### Bounded Euclidean State Space

BoundedEuclideanStateSpace(d::Int, lo::Vector, hi::Vector) = RealVectorMetricSpace(d, lo, hi, Euclidean())
UnitHypercube(d::Int) = BoundedEuclideanStateSpace(d, zeros(d), ones(d))