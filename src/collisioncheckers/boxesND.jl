export PointRobotNDBoxes

# ---------- Point Robot (amongst N-d boxes) ----------

immutable BoxBounds{N,T<:AbstractFloat}
    lo::SVector{N,T}
    hi::SVector{N,T}
end
BoxBounds(lo::AbstractVector, hi::AbstractVector) = BoxBounds(SVector(lo), SVector(hi))
BoxBounds(lohi::Matrix) = BoxBounds(lohi[:,1], lohi[:,2])
inflate{N,T}(BB::BoxBounds{N,T}, ε) = ε > 0 ? BoxBounds(BB.lo - T(ε), BB.hi + T(ε)) : BB
changeprecision{T<:AbstractFloat}(::Type{T}, BB::BoxBounds) =
    BoxBounds(changeprecision(T, BB.lo), changeprecision(T, BB.hi))

type PointRobotNDBoxes{N,T} <: SweptCollisionChecker
    boxes::Vector{BoxBounds{N,T}}
    count::Int
end
PointRobotNDBoxes{N,T}(boxes::Vector{BoxBounds{N,T}}) = PointRobotNDBoxes(boxes, 0)
PointRobotNDBoxes{T}(box_list::Vector{Matrix{T}}) = PointRobotNDBoxes(map(BoxBounds, box_list))
PointRobotNDBoxes{T}(boxes::Array{T,3}) = PointRobotNDBoxes(vec(mapslices(BoxBounds, boxes, [1,2])))
PointRobotNDBoxes(boxhcat::Matrix) = PointRobotNDBoxes(reshape(boxhcat, size(boxhcat, 1), 2, div(size(boxhcat, 2), 2)))
changeprecision{T<:AbstractFloat}(::Type{T}, CC::PointRobotNDBoxes) = PointRobotNDBoxes(changeprecision(T, CC.boxes))

is_free_state(v::AbstractVector, CC::PointRobotNDBoxes) = is_free_state(v, CC.boxes)
is_free_motion(v::AbstractVector, w::AbstractVector, CC::PointRobotNDBoxes) = (CC.count += 1;
                                                                               is_free_motion(v, w, CC.boxes))
is_free_path(P::Path, CC::PointRobotNDBoxes) = is_free_path(P, CC.boxes)

inflate{N,T}(CC::PointRobotNDBoxes{N,T}, ε) = ε > 0 ? PointRobotNDBoxes([inflate(B, T(ε)) for B in CC.boxes]) : CC
addobstacle(CC::PointRobotNDBoxes, o) = PointRobotNDBoxes(vcat(CC.boxes, BoxBounds(o)))
addblocker(CC::PointRobotNDBoxes, v::AbstractVector, r) = addobstacle(CC, BoxBounds(v-r, v+r))
closest(p::AbstractVector, CC::PointRobotNDBoxes, W::AbstractMatrix) = closest(p, CC.boxes, W)
closeR(p::AbstractVector, CC::PointRobotNDBoxes, W::AbstractMatrix, r2) = closeR(p, CC.boxes, W, r2)

plot(CC::PointRobotNDBoxes, lo = zeros(2), hi = ones(2); kwargs...) =
    map(o -> plot_rectangle(o.lo, o.hi, color = "red", edgecolor = "none",
                            xmin = lo[1], xmax = hi[1], ymin = lo[2], ymax = hi[2]; kwargs...), CC.boxes)

### Point/box obstacle checking in R^n (side note: SAT-style might be faster)

is_free_state{N}(v::AbstractVector, BB::BoxBounds{N}) = @any [!(BB.lo[i] <= v[i] <= BB.hi[i]) for i in 1:N]
is_free_state{N,T}(v::AbstractVector, BL::Vector{BoxBounds{N,T}}) = @all [is_free_state(v, BL[k]) for k in 1:length(BL)]
is_free_motion_broadphase{N}(l::AbstractVector, h::AbstractVector, BB::BoxBounds{N}) =
    @any [BB.hi[i] < l[i] || BB.lo[i] > h[i] for i in 1:N]
function is_free_motion{N}(v::AbstractVector, w::AbstractVector, BB::BoxBounds{N})
    v_to_w = w - v
    corner = blend(map(<, v, BB.lo), BB.lo, BB.hi)
    lambdas = (corner - v) ./ v_to_w
    !(@any [(@all [i == j || BB.lo[j] <= v[j] + v_to_w[j]*lambdas[i] <= BB.hi[j] for j in 1:N]) for i in 1:N])
end
function is_free_motion{N,T}(v::AbstractVector, w::AbstractVector, BL::Vector{BoxBounds{N,T}})
    bb_min = map(min, v, w)    # TODO: check back when StaticArrays has implemented min(v,w)
    bb_max = map(max, v, w)
    @all [is_free_motion_broadphase(bb_min, bb_max, BL[k]) || is_free_motion(v, w, BL[k]) for k in 1:length(BL)]
end
is_free_path{B<:BoxBounds}(P::Path, BL::Vector{B}) = @all [is_free_motion(P[i], P[i+1], BL) for i in 1:length(BL)-1]

### Closest Point

function closest(p::AbstractVector, BB::BoxBounds, W::AbstractMatrix)
    # v = Variable(length(p))
    # problem = minimize(quad_form(v-p, W), view(o,:,1) <= v, v <= view(o,:,2))
    # solve!(problem, ECOS.ECOSSolver(verbose=false))
    # return (problem.optval, vec(v.value))
    L = chol(W)
    vmin = bvls(L, L*p, BB.lo, BB.hi)
    d2min = dot(vmin - p, W*(vmin - p))
    d2min, vmin
end

function closest{N,T}(p::AbstractVector, BL::Vector{BoxBounds{N,T}}, W::AbstractMatrix)
    d2min, vmin = T(Inf), p
    for k = 1:length(BL)
        (d2, v) = closest(p, BL[k], W)
        if d2 < d2min
            d2min, vmin = d2, v
        end
    end
    d2min, vmin
end

function closeR{N,T}(p::AbstractVector, BL::Vector{BoxBounds{N,T}}, W::AbstractMatrix, r2)
    cps = [closest(p, BL[k], W) for k in 1:length(BL)]
    sort!(cps[[c[1] for c in cps] .< r2], by=first)
end