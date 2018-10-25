export Pairwise, pairwise, collision_waypoints_itr
export SteeringEdge

# Pairwise Iterator
struct Pairwise{I}
    itr::I
end
pairwise(itr) = Pairwise(itr)
Base.length(p::Pairwise) = length(p.itr) - 1
Base.eltype(::Type{Pairwise{I}}) where {I} = Tuple{eltype(I),eltype(I)}
@inline function Base.iterate(p::Pairwise, state=iterate(p.itr))
    state === nothing && return nothing
    first_item, itr_state = state
    state = iterate(p.itr, itr_state)
    state === nothing && return nothing
    (first_item, state[1]), state
end

# SteeringBVP Collision Checking
collision_waypoints_itr(f::DifferentialDynamics, x::State, controls) = waypoints_itr(f, x, controls, 8)
function is_free_edge(CC::CollisionChecker, bvp::SteeringBVP, x0::State, xf::State, controls::Nothing=nothing)
    is_free_edge(CC, bvp.dynamics, x0, bvp(x0, xf).controls)
end
function is_free_edge(CC::CollisionChecker, bvp::SteeringBVP, x0::State, xf::State, controls)
    is_free_edge(CC, bvp.dynamics, x0, controls)
end
function is_free_edge(CC::CollisionChecker, f::DifferentialDynamics, x0::State, controls)
    CC.edge_count[] += 1
    all(is_free_motion(CC, x1, x2) for (x1, x2) in pairwise(collision_waypoints_itr(f, x0, controls)))
end
## Geometric Steering
function is_free_edge(CC::CollisionChecker, bvp::GeometricSteering, x0::State, xf::State, controls::Nothing=nothing)
    CC.edge_count[] += 1
    is_free_motion(CC, x0, xf)
end
function is_free_edge(CC::CollisionChecker, bvp::GeometricSteering, x0::State, xf::State, controls)
    CC.edge_count[] += 1
    is_free_motion(CC, x0, xf)
end

# Partial Steering (for RRT and variants)
function steer_towards(bvp::SteeringBVP{<:Any,Time}, x0, xf, r) # includes GeometricSteering, TODO: interpolate_cost
    c, u = bvp(x0, xf)
    r < c ? propagate(bvp.dynamics, x0, u, r) : xf
end

# Plotting
struct SteeringEdge{BVP,S,S2C,U}
    state2config::S2C
    bvp::BVP
    x0::S
    xf::S
    controls::U
end
SteeringEdge(state2config, bvp, x0, xf) = SteeringEdge(state2config, bvp, x0, xf, bvp(x0, xf).controls)
@recipe function f(E::SteeringEdge; dims=(1, 2), num_waypoints=10, plot_endpoints=true) # , edge_color=:grey, edge_alpha=1, edge_markershape=:none
    # color :=  edge_color
    # alpha :=  edge_alpha
    # label --> ""
    x, y = dims
    plot_endpoints && @series begin
        seriestype  := :scatter
        # markershape := edge_markershape
        q0, qf = E.state2config(E.x0), E.state2config(E.xf)
        SVector(q0[x], qf[x]), SVector(q0[y], qf[y])
    end
    @series begin
        pts = waypoints(E.bvp.dynamics, E.x0, E.controls, num_waypoints)
        [p[x] for p in pts], [p[y] for p in pts]
    end
end
