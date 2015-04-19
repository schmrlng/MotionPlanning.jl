import Base.eltype, Base.getindex, Base.convert, Base.hcat
export RSState, TurningPoints, RSMetric, RSSegment, RSDiscrete, ReedsSheppStateSpace
export waypoints

### Types and Utilties
## State
immutable RSState{T<:FloatingPoint} <: AbstractState
    x::Vector2{T}
    t::T
end
RSState{T<:FloatingPoint}(x::T, y::T, t::T) = RSState(Vector2{T}(x,y), t)
eltype{T}(::RSState{T}) = T
eltype{T}(::Type{RSState{T}}) = T
convert{T}(::Type{Vector{T}}, s::RSState{T}) = convert(Vector{T}, s.x)
hcat(X::RSState...) = hcat([X[i].x for i in 1:length(X)]...)
function reorient(o::RSState, v::RSState)
    RSState(rotate(v.x - o.x, -o.t), mod2pi(v.t-o.t))   # TODO: bust out SAT2D Vector2D utilities into their own file
end

immutable TurningPoints{T<:FloatingPoint}
    r::T
    dt::T
    pts::Vector{Vector2{T}}
end
TurningPoints{T}(r::T, N::Int) = TurningPoints(r, 2pi/N, [r*Vector2(cos(x), sin(x)) for x in linspace(0,2pi,N+1)[1:end-1]])

## Metric
abstract RSMetric <: Metric     # TODO: RSInterpolated
immutable RSSegment{T<:FloatingPoint}
    t::Int          # segment type
    d::T            # segment length
end
immutable RSDiscrete{T<:FloatingPoint} <: RSMetric
    r::T
    Xmax::T
    Ymax::T
    Nx::Int
    Ny::Int
    Nt::Int
    costs::InterpGrid{T,3,BCnil,InterpQuadratic}
    paths::Array{Vector{RSSegment},3}
end
function RSDiscrete{T<:FloatingPoint}(r::T = 1., Xmax::T = 5., Ymax::T = 5., Nx::Int = 101, Ny::Int = 101, Nt::Int = 101)
    RS_costs = Array(Float64, Nx, Ny, Nt + 1)  # + 1 for cost interpolation; slice 1 and Nt + 1 are identical
    RS_paths = Array(Vector{RSSegment}, Nx, Ny, Nt)

    fname = joinpath(Pkg.dir("MotionPlanning"),"src","statespaces",
                     @sprintf("ReedsShepp_%.2f_%.2f_%d_%d_%d.gz", Xmax, Ymax, Nx, Ny, Nt))
    RS_lines = readlines(GZip.open(fname))
    for i in 1:length(RS_lines)
        sl = split(RS_lines[i])
        RS_costs[i] = r*parsefloat(T, sl[1])
        RS_paths[i] = [RSSegment(sl[2j] == "L" ? 1 : sl[2j] == "R" ? -1 : 0,
                                 (sl[2j] == "S" ? r : 1)*parsefloat(T, sl[2j+1])) for j in 1:div(length(sl) - 1, 2)]
    end
    RS_costs[:,:,Nt+1] = RS_costs[:,:,1]
    RSDiscrete(r, r*Xmax, r*Ymax, Nx, Ny, Nt, InterpGrid(RS_costs, BCnil, InterpQuadratic), RS_paths)
end

## State Space
immutable ReedsSheppStateSpace{T<:FloatingPoint} <: StateSpace
    dim::Int
    lo::Vector2{T}           # workspace only; 0 <= theta <= 2pi
    hi::Vector2{T}
    dist::RSDiscrete{T}
    r::T
    TP::TurningPoints{T}
end
ReedsSheppStateSpace(r, lo = Vector2(0.,0.), hi = Vector2(1.,1.); res = 16) = ReedsSheppStateSpace(3, lo, hi, RSDiscrete(r), r, TurningPoints(r, res))

vector_to_state{T}(v::AbstractVector{T}, SS::ReedsSheppStateSpace{T}) = RSState(Vector2(v), zero(T))   # Required for goal sampling TODO: RSGoal?
sample_space{T}(SS::ReedsSheppStateSpace{T}) = RSState(SS.lo + Vector2{T}(rand(T), rand(T)).*(SS.hi-SS.lo), convert(T, 2pi*rand(T)))   # TODO: @devec
function volume(SS::ReedsSheppStateSpace)
    # warn("TODO: what is volume for a ReedsSheppStateSpace?")
    2pi*prod(SS.hi-SS.lo)
end
function defaultNN(SS::ReedsSheppStateSpace, init)
    V = typeof(init)[init]
    ArcLength_Pruned(V, SS.dist)
end

### Steering Nuts and Bolts
function RSvec2grid(v::RSState, w::RSState, RS::RSDiscrete)
    p = reorient(v, w)
    xGrid = p.x[1]*(RS.Nx - 1)/(2*RS.Xmax) + RS.Nx/2  + .5
    yGrid = p.x[2]*(RS.Ny - 1)/(2*RS.Ymax) + RS.Ny/2  + .5
    tGrid = p.t*RS.Nt/(2pi) + 1
    xGrid, yGrid, tGrid
end
function RSvec2sub(v::RSState, w::RSState, RS::RSDiscrete)
    x, y, t = RSvec2grid(v, w, RS)
    int(x), int(y), wrap1(int(t), RS.Nt)
end
evaluate(RS::RSDiscrete, v::RSState, w::RSState) = RS.costs[RSvec2grid(v,w,RS)...]

function waypoints{T}(v::RSState{T}, s::RSSegment{T}, r::T, TP::TurningPoints)
    s.t == 0 && return Vector2{T}[v.x, v.x+s.d*Vector2(cos(v.t), sin(v.t))]
    center = v.x + sign(s.t)*Vector2(-r*sin(v.t), r*cos(v.t))
    turnpts = push!(TP.pts[1:1+ifloor(abs(s.d)/TP.dt)], Vector2(r*cos(abs(s.d)), r*sin(abs(s.d))))
    if s.t*s.d < 0
        for i in 1:length(turnpts)      # manual @devec
            turnpts[i] = Vector2(turnpts[i][1], -turnpts[i][2])
        end
    end
    [(center + sign(s.t)*rotate(p, v.t-pi/2)) for p in turnpts]
end
function waypoints{T}(v::RSState{T}, w::RSState{T}, SS::ReedsSheppStateSpace, TP::TurningPoints = SS.TP)
    pts = Array(Vector2{T}, 0)
    for s in SS.dist.paths[RSvec2sub(v, w, SS.dist)...]
        s_pts = waypoints(v, s, SS.r, TP)
        append!(pts, s_pts[1:end-1])
        v = RSState(s_pts[end], v.t + s.t*s.d)
    end
    push!(pts, w.x)
end
waypoints(i::Int, j::Int, NN::NearNeighborCache, SS::ReedsSheppStateSpace, TP::TurningPoints = SS.TP) = waypoints(NN[i], NN[j], SS, TP)

inbounds(v::RSState, SS::ReedsSheppStateSpace) = (SS.lo[1] < v.x[1] < SS.hi[1] && SS.lo[2] < v.x[2] < SS.hi[2])
inbounds(v::AbstractVector, SS::ReedsSheppStateSpace) = (SS.lo[1] < v[1] < SS.hi[1] && SS.lo[2] < v[2] < SS.hi[2])
is_free_state(v::RSState, CC::PointRobot2D, SS::ReedsSheppStateSpace) = inbounds(v, SS) && is_free_state(v.x, CC)
function is_free_motion(v::RSState, w::RSState, CC::PointRobot2D, SS::ReedsSheppStateSpace)   # TODO: inputs V, i, j instead of v, w
    wps = waypoints(v, w, SS)
    for i in 1:length(wps)-1
        (!inbounds(wps[i], SS) || !is_free_motion(wps[i], wps[i+1], CC)) && return false
    end
    true
end
# TODO: is_free_path(path, CC::PointRobot2D, SS::ReedsSheppStateSpace)

##### StateSpace

function plot_tree(SS::ReedsSheppStateSpace, NN::NearNeighborCache, A; kwargs...)
    pts = hcat(NN[find(A)]...)
    scatter(pts[1,:], pts[2,:], zorder=1; kwargs...)
    TP = TurningPoints(SS.r, 50)
    X = vcat([[hcat(waypoints(A[v], v, NN, SS, TP)...)[1,:]', nothing] for v in find(A)]...)
    Y = vcat([[hcat(waypoints(A[v], v, NN, SS, TP)...)[2,:]', nothing] for v in find(A)]...)
    plt.plot(X, Y, linewidth=.5, linestyle="-", zorder=1; kwargs...)
end

function plot_path(SS::ReedsSheppStateSpace, NN::NearNeighborCache, sol; kwargs...)
    TP = TurningPoints(SS.r, 50)
    wps = hcat([hcat(waypoints(sol[i], sol[i+1], NN, SS, TP)...) for i in 1:length(sol)-1]...)
    length(sol) > 1 && plot_path(wps; kwargs...)
    # plt.quiver([wps[row,1:3:end]' for row in 1:4]..., zorder=5, width=.003, headwidth=8)
end