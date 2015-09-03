export BoundedEuclideanStateSpace, UnitHypercube

## Metric Space Instantiation 
BoundedEuclideanStateSpace(d::Int, lo::Vector, hi::Vector) = RealVectorStateSpace(d, lo, hi, Euclidean(), Identity())
UnitHypercube(d::Int) = BoundedEuclideanStateSpace(d, zeros(d), ones(d))

waypoints{S}(v::S, w::S, ::Euclidean) = S[v, w]
function time_waypoint{S}(v::S, w::S, M::Euclidean, t)
    n = evaluate(M, v, w)
    v + min(t/n, 1)*(w - v)
end
cost_waypoint{S}(v::S, w::S, M::Euclidean, c) = time_waypoint(v, w, M, c)
defaultNN(::Euclidean, init) = EuclideanNN_KDTree(typeof(init)[init])