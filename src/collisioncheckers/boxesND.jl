### Obstacle Types

abstract ObstacleSet

immutable ObstacleList <: ObstacleSet
    list::Vector{ObstacleSet}
end

typealias Matrix3D{T} Array{T,3}
immutable AABoxes <: ObstacleSet
    boxes::Matrix3D
end

AABoxes(box_hcat::Matrix) = AABoxes(reshape(box_hcat, size(box_hcat, 1), 2, div(size(box_hcat, 2), 2)))
AABoxes{T}(box_list::Vector{Matrix{T}}) = AABoxes(hcat(box_list...))
volume(obs::ObstacleSet) = volume_naive(obs)

### ObstacleList

is_free_state(v::Vector, obs_list::ObstacleList) = all([is_free_state(v, obs) for obs in obs_list.list])
is_free_motion(v::Vector, w::Vector, obs_list::ObstacleList) = all([is_free_motion(v, w, obs) for obs in obs_list.list])
is_free_path{T}(path::Vector{Vector{T}}, obs_list::ObstacleList) = all([is_free_path(path, obs) for obs in obs_list.list])
closest_obs_pt(v::Vector, obs_list::ObstacleList) = mapreduce(obs -> closest_obs_pt(v, obs), (x,y) -> (x[1] < y[1] ? x : y), obs_list.list)
volume_naive(obs_list::ObstacleList) = sum([volume_naive(obs) for obs in obs_list.list])
plot_obstacles(obs_list::ObstacleList) = [plot_obstacles(obs) for obs in obs_list.list]

### AABoxes

is_free_state(v::Vector, obs::AABoxes) = is_free_state(v, obs.boxes)
is_free_motion(v::Vector, w::Vector, obs::AABoxes) = is_free_motion(v, w, obs.boxes)
is_free_path{T}(path::Vector{Vector{T}}, obs::AABoxes) = is_free_path(path, obs.boxes)
closest_obs_pt(v::Vector, obs::AABoxes) = closest_obs_pt(v, obs.boxes)
volume_naive(obs::AABoxes) = sum(mapslices(prod, obs.boxes[:,2,:] - obs.boxes[:,1,:], [1,2]))
plot_obstacles(obs::AABoxes) = mapslices(o -> plot_rectangle(o, color = "red"), obs.boxes, [1,2])

### Point/box obstacle checking in R^n

function is_free_state(v::Vector, o::AbstractMatrix)
    for i = 1:size(o,1)
        if !(o[i,1] <= v[i] <= o[i,2])
            return true
        end
    end
    return false
end

function is_free_state(v::Vector, obs::Matrix3D)
    for k = 1:size(obs,3)
        !is_free_state(v, view(obs,:,:,k)) && return false
    end
    return true
end

function is_free_motion_broadphase(bb_min::Vector, bb_max::Vector, obs::Matrix3D, k::Int)
    for i = 1:size(obs,1)
        if obs[i,2,k] < bb_min[i] || obs[i,1,k] > bb_max[i]
            return true
        end
    end
    return false
end

function face_contains_projection(v::Vector, v_to_w::Vector, lambda::Number, j::Int, o::AbstractMatrix)
    for i = 1:size(o,1)
        if i != j && !(o[i,1] <= v[i] + v_to_w[i]*lambda <= o[i,2])
            return false
        end
    end
    return true
end

function is_free_motion(v::Vector, w::Vector, o::AbstractMatrix)
    v_to_w = w - v
    @devec lambdas = (blend(v .< o[:,1], o[:,1], o[:,2]) .- v) ./ v_to_w
    for i in 1:size(o,1)
        face_contains_projection(v, v_to_w, lambdas[i], i, o) && return false
    end
    return true
end

function is_free_motion(v::Vector, w::Vector, obs::Matrix3D)
    bb_min = min(v,w)
    bb_max = max(v,w)
    for k = 1:size(obs,3)
        if !is_free_motion_broadphase(bb_min, bb_max, obs, k)
            if !is_free_motion(v, w, view(obs,:,:,k))
                return false
            end
        end
    end
    return true
end


function is_free_path_naive{T}(path::Vector{Vector{T}}, obs::Matrix3D)
    for i in 1:length(path)-1
        if !is_free_motion(path[i], path[i+1], obs)
            return false
        end
    end
    return true
end
is_free_path{T}(path::Vector{Vector{T}}, obs::Matrix3D) = is_free_path_naive(path, obs)

function closest_obs_pt(v::Vector, o::AbstractMatrix)
    x = [clamp(v[i], o[i,1], o[i,2]) for i in 1:length(v)]
    return (norm(x-v), x)
end

function closest_obs_pt(v::Vector, obs::Matrix3D)
    min_d = Inf
    min_x = v
    for k = 1:size(obs,3)
        (d, x) = closest_obs_pt(v, view(obs,:,:,k))
        if d < min_d
            min_d = d
            min_x = x
        end
    end
    return (min_d, min_x)
end