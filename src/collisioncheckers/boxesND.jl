export PointRobotNDBoxes

# ---------- Point Robot (amongst N-d boxes) ----------

typealias Matrix3D{T} Array{T,3}
type PointRobotNDBoxes <: SweptCollisionChecker
    boxes::Matrix3D
    count::Int
end
PointRobotNDBoxes(boxes::Matrix3D) = PointRobotNDBoxes(boxes, 0)
PointRobotNDBoxes(box_hcat::Matrix) = PointRobotNDBoxes(reshape(box_hcat, size(box_hcat, 1), 2, div(size(box_hcat, 2), 2)))
PointRobotNDBoxes{T}(box_list::Vector{Matrix{T}}) = PointRobotNDBoxes(hcat(box_list...))

is_free_state(v::AbstractVector, CC::PointRobotNDBoxes) = is_free_state(v, CC.boxes)
is_free_motion(v::AbstractVector, w::AbstractVector, CC::PointRobotNDBoxes) = is_free_motion(v, w, CC.boxes)
is_free_path(path::Path, CC::PointRobotNDBoxes) = is_free_path(path, CC.boxes)

inflate(CC::PointRobotNDBoxes, eps) = eps > 0 ? PointRobotNDBoxes(CC.boxes .+ [-eps eps]) : CC  # TODO: rounded corners/edges
closest(p::AbstractVector, CC::PointRobotNDBoxes, W::AbstractMatrix) = closest(p, CC.boxes, W)
close(p::AbstractVector, CC::PointRobotNDBoxes, W::AbstractMatrix, r2) = close(p, CC.boxes, W, r2)

plot(CC::PointRobotNDBoxes, lo = zeros(2), hi = ones(2); kwargs...) =
    mapslices(o -> plot_rectangle(o, color = "red", edgecolor = "none"; kwargs...), CC.boxes, [1,2])

### Point/box obstacle checking in R^n (side note: SAT-style might be faster)
# TODO: devectorized @all/@any or something of the sort

function is_free_state(v::AbstractVector, o::AbstractMatrix)
    for i = 1:size(o,1)
        !(o[i,1] <= v[i] <= o[i,2]) && return true
    end
    return false
end

function is_free_state(v::AbstractVector, obs::Matrix3D)
    for k = 1:size(obs,3)
        !is_free_state(v, view(obs,:,:,k)) && return false
    end
    return true
end

function is_free_motion_broadphase(bb_min::AbstractVector, bb_max::AbstractVector, obs::Matrix3D, k::Int)
    for i = 1:size(obs,1)
        (obs[i,2,k] < bb_min[i] || obs[i,1,k] > bb_max[i]) && return true
    end
    return false
end

function face_contains_projection(v::AbstractVector, v_to_w::AbstractVector, lambda::Number, j::Int, o::AbstractMatrix)
    for i = 1:size(o,1)
        (i != j && !(o[i,1] <= v[i] + v_to_w[i]*lambda <= o[i,2])) && return false
    end
    return true
end

function is_free_motion(v::AbstractVector, w::AbstractVector, o::AbstractMatrix)
    v_to_w = w - v
    @devec lambdas = (blend(v .< o[:,1], o[:,1], o[:,2]) .- v) ./ v_to_w
    for i in 1:size(o,1)
        face_contains_projection(v, v_to_w, lambdas[i], i, o) && return false
    end
    return true
end

function is_free_motion(v::AbstractVector, w::AbstractVector, obs::Matrix3D)
    bb_min = min(v,w)
    bb_max = max(v,w)
    for k = 1:size(obs,3)
        !is_free_motion_broadphase(bb_min, bb_max, obs, k) && !is_free_motion(v, w, view(obs,:,:,k)) && return false
    end
    return true
end

function is_free_path(path::Path, obs::Matrix3D)
    for i in 1:length(path)-1
        !is_free_motion(path[i], path[i+1], obs) && return false
    end
    return true
end

function closest(p::AbstractVector, o::AbstractMatrix, W::AbstractMatrix)
    # v = Variable(length(p))
    # problem = minimize(quad_form(v-p, W), view(o,:,1) <= v, v <= view(o,:,2))
    # solve!(problem, ECOS.ECOSSolver(verbose=false))
    # return (problem.optval, vec(v.value))
    L = chol(W)
    vmin = bvls(L, L*p, view(o,:,1), view(o,:,2))
    d2min = dot(vmin - p, W*(vmin - p))
    d2min, vmin
end

function closest(p::AbstractVector, obs::Matrix3D, W::AbstractMatrix)
    d2min, vmin = Inf, p
    for k = 1:size(obs,3)
        (d2, v) = closest(p, view(obs,:,:,k), W)
        if d2 < d2min
            d2min, vmin = d2, v
        end
    end
    d2min, vmin
end

function close(p::AbstractVector, obs::Matrix3D, W::AbstractMatrix, r2)
    cps = [closest(p, view(obs,:,:,k), W) for k in 1:size(obs,3)]
    sort!(cps[[c[1] for c in cps] .< r2], by=first)
end