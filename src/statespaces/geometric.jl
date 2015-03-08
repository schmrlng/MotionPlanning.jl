export RealVectorMetricSpace, BoundedEuclideanStateSpace, UnitHypercube
export shortcut, adaptive_shortcut

immutable RealVectorMetricSpace{T<:FloatingPoint} <: RealVectorStateSpace
    dim::Int
    lo::Vector{T}
    hi::Vector{T}
    dist::Metric
end

vector_to_state{T}(v::AbstractVector{T}, SS::RealVectorMetricSpace) = SS.dim == 2 ? convert(Vector2{T}, v) : v
sample_space(SS::RealVectorMetricSpace) = vector_to_state(SS.lo + rand(SS.dim).*(SS.hi-SS.lo), SS)
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

### ADAPTIVE-SHORTCUT (Hsu 2000)

function shortcut(path::Path, CC::CollisionChecker)
    N = length(path)
    if N == 2
        return path
    end
    if is_free_motion(path[1], path[end], CC)
        return path[[1,end]]
    end
    mid = iceil(N/2)
    return [shortcut(path[1:mid], CC)[1:end-1], shortcut(path[mid:end], CC)]
end

function cut_corner(v1::AbstractVector, v2::AbstractVector, v3::AbstractVector, CC::CollisionChecker)
    m1 = (v1 + v2)/2
    m2 = (v3 + v2)/2
    while ~is_free_motion(m1, m2, CC)
        m1 = (m1 + v2)/2
        m2 = (m2 + v2)/2
    end
    return typeof(v1)[v1, m1, m2, v3]
end

function adaptive_shortcut(path::Path, CC::CollisionChecker, iterations::Int = 10)
    while (short_path = shortcut(path, CC)) != path
        path = short_path
    end
    for i in 1:iterations
        path = [path[1:1], vcat([cut_corner(path[j-1:j+1]..., CC)[2:3] for j in 2:length(path)-1]...), path[end:end]]
        while (short_path = shortcut(path, CC)) != path
            path = short_path
        end
    end
    return path, sum(mapslices(norm, diff(hcat(path...), 2), 1))
end