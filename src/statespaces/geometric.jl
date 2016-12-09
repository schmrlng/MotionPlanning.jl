export BoundedEuclideanStateSpace, UnitHypercube

### Metric Typedefs/Evaluation
evaluate(M::Euclidean, v::SVector, w::SVector) = norm(w - v)
colwise{S<:AbstractVector}(d::Euclidean, v::AbstractVector, W::Vector{S}) = colwise(d, v, statevec2mat(W))
colwise{S<:AbstractVector}(d::Euclidean, W::Vector{S}, v::AbstractVector) = colwise(d, v, statevec2mat(W))

### Metric Space Instantiation
# TODO: implement show method with dynamics and cost function
BoundedEuclideanStateSpace(lo::AbstractVector, hi::AbstractVector) =
    BoundedStateSpace(SVector(lo), SVector(hi), Euclidean(), Identity())
UnitHypercube{T<:AbstractFloat}(d::Int, ::Type{T} = Float64) = BoundedEuclideanStateSpace(zeros(T,d), ones(T,d))

helper_data_structures{S}(V::Vector{S}, M::Euclidean) = (TreeDistanceDS(KDTree(V, M; reorder=false)),
                                                         EmptyControlCache())

### Steering
propagate(M::Euclidean, v::AbstractVector, u::StepControl) = v + u.t*typeof(v)(u.u)
steering_control(M::Euclidean, v::AbstractVector, w::AbstractVector) = StepControl(evaluate(M, v, w), normalize(w - v))
collision_waypoints(::Euclidean, v::AbstractVector, w::AbstractVector) = (v, w)
# TODO: specialized methods that don't check ends
