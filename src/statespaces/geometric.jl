export BoundedEuclideanStateSpace, UnitHypercube

### Metric Typedefs/Evaluation
evaluate(M::Euclidean, v::Vec, w::Vec) = norm(w - v)
colwise{S<:State}(d::Euclidean, v::Vec, W::Vector{S}) = colwise(d, dense(v), statevec2mat(W))
colwise{S<:State}(d::Euclidean, W::Vector{S}, v::Vec) = colwise(d, dense(v), statevec2mat(W))
colwise(d::Euclidean, v::Vec, W::Matrix) = colwise(d, dense(v), W)
colwise(d::Euclidean, W::Matrix, v::Vec) = colwise(d, dense(v), W)

### Metric Space Instantiation
# TODO: implement show method with dynamics and cost function
BoundedEuclideanStateSpace{d}(lo::Vec{d}, hi::Vec{d}) = RealVectorStateSpace(lo, hi, Euclidean(), Identity())
UnitHypercube{T<:AbstractFloat}(d::Int, ::Type{T} = Float64) = BoundedEuclideanStateSpace(zero(Vec{d,T}), one(Vec{d,T}))

helper_data_structures{S}(V::Vector{S}, M::Euclidean) = (TreeDistanceDS(KDTree(statevec2mat(V), M; reorder=false)),
                                                         EmptyControlCache())
inrange(tree::NearestNeighbors.NNTree, v::Vec, args...) = inrange(tree, dense(v), args...)

### Steering
propagate(M::Euclidean, v::Vec, u::StepControl) = v + u.t*u.u
steering_control(M::Euclidean, v::Vec, w::Vec) = StepControl(evaluate(M, v, w), normalize(w - v))
collision_waypoints(::Euclidean, v::Vec, w::Vec) = (v, w)
# TODO: specialized methods that don't check ends

# time_waypoint(M::Euclidean, v, w, t) = v + min(t/evaluate(M, v, w), 1)*(w - v)
# cost_waypoint(M::Euclidean, v, w, c) = time_waypoint(M, v, w, c)