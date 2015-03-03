export RealVectorMetricSpace, BoundedEuclideanStateSpace, UnitHypercube
export shortcut, adaptive_shortcut

immutable RealVectorMetricSpace <: RealVectorStateSpace
    dim::Int
    lo::Vector
    hi::Vector
    dist::Metric
end

vector_to_state(v::Vector, SS::RealVectorMetricSpace) = SS.dim == 2 ? Vector2(v) : v
sample_space(SS::RealVectorMetricSpace) = (v = SS.lo + rand(SS.dim).*(SS.hi-SS.lo); SS.dim == 2 ? Vector2(v) : v)

### Bounded Euclidean State Space

BoundedEuclideanStateSpace(d::Int, lo::Vector, hi::Vector) = RealVectorMetricSpace(d, lo, hi, Euclidean())
UnitHypercube(d::Int) = BoundedEuclideanStateSpace(d, zeros(d), ones(d))

volume(SS::RealVectorMetricSpace) = prod(SS.hi-SS.lo)
steer(SS::RealVectorMetricSpace, v::Vector, w::Vector, eps::Float64, distvw = norm(w - v)) = v + (w - v) * min(eps/distvw, 1)

pairwise_distances{T}(V::Vector{Vector{T}}, SS::RealVectorMetricSpace, r_bound::Float64) = pairwise(SS.dist, hcat(V...))

### ADAPTIVE-SHORTCUT (Hsu 2000)

function shortcut{T}(path::Vector{Vector{T}}, obs::ObstacleSet)
    N = length(path)
    if N == 2
        return path
    end
    if is_free_motion(path[1], path[end], obs)
        return path[[1,end]]
    end
    mid = iceil(N/2)
    return [shortcut(path[1:mid], obs)[1:end-1], shortcut(path[mid:end], obs)]
end

function cut_corner(v1::Vector, v2::Vector, v3::Vector, obs::ObstacleSet)
    m1 = (v1 + v2)/2
    m2 = (v3 + v2)/2
    while ~is_free_motion(m1, m2, obs)
        m1 = (m1 + v2)/2
        m2 = (m2 + v2)/2
    end
    return typeof(v1)[v1, m1, m2, v3]
end

function adaptive_shortcut{T}(path::Vector{Vector{T}}, obs::ObstacleSet, iterations::Int = 10)
    while (short_path = shortcut(path, obs)) != path
        path = short_path
    end
    for i in 1:iterations
        path = [path[1:1], vcat([cut_corner(path[j-1:j+1]..., obs)[2:3] for j in 2:length(path)-1]...), path[end:end]]
        while (short_path = shortcut(path, obs)) != path
            path = short_path
        end
    end
    return path, sum(mapslices(norm, diff(hcat(path...), 2), 1))
end