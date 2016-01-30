export BoundedEuclideanStateSpace, UnitHypercube

## Metric Space Instantiation 
BoundedEuclideanStateSpace{d}(lo::Vec{d}, hi::Vec{d}) = RealVectorStateSpace(lo, hi, Euclidean(), Identity())
UnitHypercube{T<:AbstractFloat}(d::Int, ::Type{T} = Float64) = BoundedEuclideanStateSpace(zero(Vec{d,T}), one(Vec{d,T}))

colwise{S<:Vec}(d::Euclidean, v::S, W::Vector{S}) = colwise(d, full(v), hcat(W...))
colwise(d::Euclidean, v::Vec, W::Matrix) = colwise(d, full(v), W)
helper_data_structures{S}(V::Vector{S}, M::Euclidean) = (TreeDistanceDS(KDTree(hcat(V...), M; reorder=false)), EmptyControlDS())

waypoints{S}(::Euclidean, v::S, w::S) = S[v, w]
function time_waypoint(M::Euclidean, v, w, t)
    n = evaluate(M, v, w)
    v + min(t/n, 1)*(w - v)
end
cost_waypoint(M::Euclidean, v, w, c) = time_waypoint(M, v, w, c)